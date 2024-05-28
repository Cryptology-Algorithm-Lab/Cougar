use blstrs::*;

use crate::Leopard::utils::*;
use crate::Leopard::setup::LeopardParams;
use merlin::Transcript;
use crate::Leopard::transcript::*;
use crate::Leopard::msm_pippenger::Pippenger;
use group::ff::Field;
use std::time::Instant;
use group::Group;
use std::mem::size_of;

#[derive(Debug)]
pub struct InnerProductProof {
    pub(crate) P: Gt,
    pub(crate) L_vec: Vec<Gt>,
    pub(crate) R_vec: Vec<Gt>,
    pub(crate) a: Scalar,
    pub(crate) b: Scalar,
}

impl InnerProductProof {
    pub fn proof(
        transcript: &mut Transcript,
        LP: &LeopardParams,
        a_vec: Vec<Scalar>,
        b_vec: Vec<Scalar>,
    ) -> InnerProductProof{

        let mut m = LP.m; let mut n = LP.n;
        let mut g = &mut LP.g.clone()[..];
        let mut h = &mut LP.h.clone()[..];
        let mut hH = &mut LP.H.clone()[..];
        let mut a = &mut a_vec.clone()[..];
        let mut b = &mut b_vec.clone()[..];

        let gamma: Scalar = transcript.challenge_scalar(b"gamma");
        let U: Gt = (LP.U * gamma).into();
        let c = inner_prod(a, b);
        let P = calc_LR1(g, h, hH, a, b) + U * c;
        let lg_m = m.next_power_of_two().trailing_zeros() as usize;
        let lg_n = n.next_power_of_two().trailing_zeros() as usize;
        
        
        // Setup
        let mut L_vec: Vec<Gt> = Vec::with_capacity((lg_m+lg_n));
        let mut R_vec: Vec<Gt> = Vec::with_capacity((lg_m+lg_n));
        
        // Column Reduction
        while m!= 1 {
            m >>= 1;

            // Step 1. Divide
            let (a_p, a_m) = a.split_at_mut(m*n);
            let (b_p, b_m) = b.split_at_mut(m*n);
            let (g_p, g_m) = g.split_at_mut(m);
            let (h_p, h_m) = h.split_at_mut(m);
            
            // Step 2. Calculate L & R
            let c_L: Scalar = inner_prod(a_p, b_m);
            let c_R: Scalar = inner_prod(a_m, b_p);
            
            let L: Gt = calc_LR1(g_m, h_p, hH, a_p, b_m) + U * c_L;
            let R: Gt = calc_LR1(g_p, h_m, hH, a_m, b_p) + U * c_R;
            
            // Step 3. Push L, R & Update Transcript
            L_vec.push(L);
            R_vec.push(R);
            
            transcript.append_point(b"L", &L);
            transcript.append_point(b"R", &R);
            
            
            // Step 4. Get challenge
            let x = transcript.challenge_scalar(b"x");
            let x_inv = x.invert().unwrap();

            // Step 5. Update Parameters
            for i in 0..(m*n) {
                a_p[i] = a_p[i] * x + a_m[i] * x_inv;
                b_p[i] = b_p[i] * x_inv + b_m[i] * x;
            }
            
            for i in 0..m {
                g_p[i] = (g_p[i] * x_inv + g_m[i] * x).into();
                h_p[i] = (h_p[i] * x + h_m[i] * x_inv).into();
            }

            a = a_p;
            b = b_p;
            g = g_p;
            h = h_p;
                        
        }
        
        // Row Reduction
        // Before that, we first need to copy H as two parts
        let mut gH = &mut LP.H.clone()[..];
        
        while n!=1 {
            n >>= 1;
            
            // Step 1. Divide
            let (a_L, a_R) = a.split_at_mut(m*n);
            let (b_L, b_R) = b.split_at_mut(m*n);
            let (gH_L, gH_R) = gH.split_at_mut(n);
            let (hH_L, hH_R) = hH.split_at_mut(n);

            // Step 2. Calculate L & R
            let c_L: Scalar = inner_prod(&a_L, &b_R);
            let c_R: Scalar = inner_prod(&a_R, &b_L);
            
            let L: Gt = calc_LR2(&g, &h, &gH_R, &hH_L, &a_L, &b_R) + U * c_L;
            let R: Gt = calc_LR2(&g, &h, &gH_L, &hH_R, &a_R, &b_L) + U * c_R;
            
            // Step 3. Push L, R & Update Transcript
            L_vec.push(L);
            R_vec.push(R);
            
            transcript.append_point(b"L", &L);
            transcript.append_point(b"R", &R);

            // Step 4. Get challenge
            let x = transcript.challenge_scalar(b"x");
            let x_inv = x.invert().unwrap();

            // Step 5. Update Parameters
            for i in 0..m*n {
                a_L[i] = a_L[i] * x + a_R[i] * x_inv;
                b_L[i] = b_L[i] * x_inv + b_R[i] * x;
            }

            for i in 0..n {
                gH_L[i] = (gH_L[i] * x_inv + gH_R[i] * x).into();
                hH_L[i] = (hH_L[i] * x + hH_R[i] * x_inv).into();
            }
            
            a = a_L;
            b = b_L;
            gH = gH_L;
            hH = hH_L;
        }
        
        // Done!
        InnerProductProof {
            P: P,
            L_vec: L_vec,
            R_vec: R_vec, 
            a: a[0], 
            b: b[0]}
    }

    pub(crate) fn verification_scalar(&self, m:usize, n:usize, transcript: &mut Transcript)
    -> (Vec<Scalar>, Vec<Scalar>, Vec<Scalar>, Vec<Scalar>){
        let lg_mn = self.L_vec.len();
        let lg_m = m.next_power_of_two().trailing_zeros() as usize;
        let lg_n = n.next_power_of_two().trailing_zeros() as usize;
        assert_eq!(lg_mn, lg_m + lg_n);
        
        
        // Step 1. Check challenges: pass 
        let mut challenges = Vec::with_capacity(lg_mn);
        for (L,R) in self.L_vec.iter().zip(self.R_vec.iter()) {
            transcript.validate_and_append_point(b"L", L);
            transcript.validate_and_append_point(b"R", R);
            challenges.push(transcript.challenge_scalar(b"x"));
        }
        
        // Step 2. Basic Setting
        let mut challenges_inv = challenges.clone();
        let (inv_left, inv_right) = batch_invert(lg_m, lg_n, &mut challenges_inv);

        // Step 4. Compute x_sq & x_inv_sq
        challenges = challenges.iter().map(|x| x.square()).collect();
        challenges_inv = challenges_inv.iter().map(|x| x.square()).collect();
        
        let mut s_L = Vec::with_capacity(m);
        let mut s_R = Vec::with_capacity(n);
        
        s_L.push(inv_left);
        s_R.push(inv_right);
        

        for i in 1..m {
            let lg_i = (32 - 1 - (i as u32).leading_zeros()) as usize;
            let k = 1 << lg_i;
            let u_lg_i_sq = challenges[(lg_mn - lg_n - 1) - lg_i];
            s_L.push(s_L[i-k] * u_lg_i_sq);
        }

        for i in 1..n {
            let lg_i = (32 - 1 - (i as u32).leading_zeros()) as usize;
            let k = 1 << lg_i;
            let u_lg_i_sq = challenges[(lg_mn - 1) - lg_i];
            s_R.push(s_R[i-k] * u_lg_i_sq);
        }
        
        (challenges, challenges_inv, s_L, s_R)
    }
    
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        LP: &LeopardParams,
    ) {
        // Step 1. Scalar Preprocessing
        let (m, n) = (LP.m, LP.n);
        let lgm = m.next_power_of_two().trailing_zeros() as usize;
        let lgn= n.next_power_of_two().trailing_zeros() as usize; 

        let gamma = transcript.challenge_scalar(b"gamma");
        let U: Gt = (LP.U * gamma).into();
        
        let (x_sq, x_inv_sq, s_L, s_R) = self.verification_scalar(m, n, transcript);
        
        // Sanity Check
        assert_eq!(x_sq.len(), lgm + lgn);
        assert_eq!(x_inv_sq.len(), lgm + lgn);
        assert_eq!(s_L.len(), m);
        assert_eq!(s_R.len(), n);
        
        
        // We need to split them for column reduction & row reduction
        let (x_sq_L, x_sq_R) = x_sq.split_at(lgm);
        let (x_inv_sq_L, x_inv_sq_R) = x_inv_sq.split_at(lgm);
        let (s_inv_L, s_inv_R) = (s_L.iter().rev(), s_R.clone().into_iter().rev());
        
        // Step 2. g, h, H precomputing
        let gas: Vec<Scalar> = s_L.iter().map(|s| self.a * s).collect();
        let hbsi: Vec<Scalar> = s_inv_L.into_iter().map(|s| self.b * s).collect();
        
        let g: Vec<G1Projective> = LP.g.iter().map(|x| x.into()).collect();
        let h: Vec<G1Projective> = LP.h.iter().map(|x| x.into()).collect();
        let H: Vec<G2Projective> = LP.H.iter().map(|x| x.into()).collect();
        let s_inv_R: Vec<Scalar> = s_inv_R.into_iter().map(|x| x).collect();
        
        // Step 3. Calulate Left & Right
        // Left: Optimization with MEXP
        let pair_ll: G1Affine = G1Projective::multi_exp(&g,  &gas).into();
        let pair_lr: G2Affine = G2Projective::multi_exp(&H, &s_R).into();
        let pair_rl: G1Affine = G1Projective::multi_exp(&h, &hbsi).into();
        let pair_rr: G2Affine = G2Projective::multi_exp(&H, &s_inv_R).into();
        
        
        let left = pairing(&pair_ll, &pair_lr) + 
                    pairing(&pair_rl, &pair_rr) + 
                    U * (self.a * self.b);

        let right = self.P + Gt::pippenger(
            self.L_vec.clone().into_iter().chain(self.R_vec.clone().into_iter()),
            x_sq.clone().into_iter().chain(x_inv_sq.clone().into_iter())
            );
      
        assert_eq!(left, right);
    }
}