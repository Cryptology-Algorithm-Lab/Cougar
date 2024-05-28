use halo2curves::grumpkin::{
    G1Affine as GG1Affine, Fr as GScalar, G1 as GG1
};

use halo2_proofs::poly::{
    EvaluationDomain,
    ExtendedLagrangeCoeff,
    LagrangeCoeff,
    Coeff,
    Polynomial,
};

use halo2curves::bn256::{
    G1Affine as BG1Affine, G2Affine as BG2Affine, G1 as BG1, G2 as BG2,
    Fr as BScalar, Gt, G2Prepared as BG2Prepared
};

use halo2curves::bn256::{pairing, multi_miller_loop};
use halo2curves::pairing::MillerLoopResult;

use crate::Cougar::setup::{
    CougarParams,
    LeopardParams,
};

use crate::Cougar::utils::{
    inner_prod, Leopard_commit, inner_commit, eval_as_poly,
};

use crate::Cougar::execution_trace::{
    op_PTROW,
    op_PTCOL,
};

use crate::Cougar::custom_gates::{
    constructor_T,
    constructor_Z,
};

use crate::Cougar::transcript::{
    Transcript,
    TranscriptWrite,
    Blake2bWrite,
    Challenge255,
};

use crate::Cougar::LeopardPCS::leopard_prover::{
    leopard_merged_prover,
};


use crate::Cougar::permutation::{
    get_id_and_permpolys,
    get_z_poly,
};


use ff::Field;

pub fn divide_to_coord(points: &[GG1Affine]) -> (Vec<BScalar>, Vec<BScalar>)
{
    let n = points.len();
    let mut x_vec: Vec<BScalar> = Vec::with_capacity(n);
    let mut y_vec: Vec<BScalar> = Vec::with_capacity(n);
    
    for i in 0..n {
        x_vec.push(points[i].x);
        y_vec.push(points[i].y);
    }
    
    (x_vec, y_vec)
}

pub fn cougar_encode(
    E: &EvaluationDomain<BScalar>,
    coeffs: &[BScalar],
    D: usize,
    offset: usize,
    stepsize: usize,
) -> Polynomial<BScalar, Coeff>
{
    let n = coeffs.len();
    let mut ret: Vec<BScalar> = vec![BScalar::ZERO; D];
    
    let mut cursor = offset;
    
    for i in 0..n {
        ret[cursor] = coeffs[i];
        cursor += stepsize;
    }
    let mut tmp = E.lagrange_from_vec(ret);
    E.lagrange_to_coeff(tmp)
}

// Leopard + Cougar Encoding
pub fn Leopard_with_encoding (
    E: &EvaluationDomain<BScalar>,
    ck_pc: &LeopardParams,
    points: &[GG1Affine],
    D: usize,
    offset: usize,
    stepsize: usize,
) -> (Gt, Gt)
{
    // Step 1. Divide
    let (x_vec, y_vec) = divide_to_coord(points);
    
    // Step 2. Encode
    let f_x = cougar_encode(E, &x_vec, D, offset, stepsize);
    let f_y = cougar_encode(E, &y_vec, D, offset, stepsize);
    
    // Step 3. Apply Leopard
    (Leopard_commit(ck_pc, &f_x[..]), Leopard_commit(ck_pc, &f_y[..]))
}


pub fn manual_recode(ext: &mut [Vec<BScalar>], points: &[GG1Affine], stepsize: usize, offset: usize) 
-> usize
{
    let n = points.len();
    let mut cursor = offset;
    
    for i in 0..n {
        ext[5][cursor] = points[i].x;
        ext[6][cursor] = points[i].y;
        cursor += stepsize;
    }
    let diff = cursor - offset;
    diff
}

pub struct A_poly {
    // Coeffs
    pub coeff_x: Vec<BScalar>,
    pub coeff_y: Vec<BScalar>,
    // Records: Start Index, Stepsize, # of Records
    pub record: Vec<(usize, usize, usize)>,
    // Values
    pub stepsize: usize,
    pub idx: usize,
}

impl A_poly {
    pub fn new(D: usize, stepsize: usize) -> Self {
         Self {
             coeff_x: vec![BScalar::ZERO; D],
             coeff_y: vec![BScalar::ZERO; D],
             record: Vec::new(),
             stepsize: stepsize,
             idx: 0,
        }
    }
    
    pub fn register(&mut self, P_vec: &[GG1Affine]) {
        let n = P_vec.len();
        let (x_vec, y_vec) = divide_to_coord(P_vec);
        self.record.push((self.idx, self.stepsize, n));
        
        for i in 0..n {
            self.coeff_x[self.idx] = x_vec[i];
            self.coeff_y[self.idx] = y_vec[i];
            self.idx += self.stepsize;
        }
        
    }
    
    pub fn register_three(&mut self, L_vec: &[GG1Affine], R_vec: &[GG1Affine], P_vec: &[GG1Affine]) {
        let n = L_vec.len();
        let (l_x, l_y) = divide_to_coord(L_vec);
        let (r_x, r_y) = divide_to_coord(R_vec);
        let (p_x, p_y) = divide_to_coord(P_vec);
        self.record.push((self.idx, 3 * self.stepsize, n));
        self.record.push((self.idx + self.stepsize, 3 * self.stepsize, n));
        self.record.push((self.idx + 2 * self.stepsize, 3 * self.stepsize, n));
        
        for i in 0..n {
            self.coeff_x[self.idx] = l_x[i];
            self.coeff_y[self.idx] = l_y[i];
            self.idx += self.stepsize;
            
            self.coeff_x[self.idx] = r_x[i];
            self.coeff_y[self.idx] = r_y[i];
            self.idx += self.stepsize;
            
            self.coeff_x[self.idx] = p_x[i];
            self.coeff_y[self.idx] = p_y[i];
            self.idx += self.stepsize;            
        }
    }
    
    pub fn count(&self) -> usize {
        self.record.len()
    }
    
}


pub fn cougar_prover(
    // Commitment Parameters
    CP: &CougarParams,
    // Transcript for Fiat-Shamir
    transcript: &mut Blake2bWrite<Vec<u8>, BG1Affine, Challenge255<BG1Affine>>,
    // Witnesses
    a_vec: Vec<GScalar>,
    b_vec: Vec<GScalar>,
)
{
    // Setting Up
    let mut m = CP.m; let mut n = CP.n; let d = CP.d; let D = CP.D;
    let ck_pc = &CP.ck_pc;
    
    // Evaluation Domain for xi
    let eval_xi: &EvaluationDomain<BScalar> = &CP.eval_xi;    
    let eval_zeta: &EvaluationDomain<BScalar> = &CP.eval_zeta;    
    
    // Execution Trace (of size D)
    let mut ext: Vec<Vec<BScalar>> = vec![vec![BScalar::ZERO ; D]; 7];
    
    // Permutation Table (position-wise encoding)
    // In our application, transpositions suffice!
    let mut perm: Vec<Vec<(usize, usize)>> = vec![Vec::new(); 7];
    
    // Parameter PT4.Row & PT4.Col
    let lg_m = m.next_power_of_two().trailing_zeros() as usize;
    let lg_n = n.next_power_of_two().trailing_zeros() as usize;
    
    
    // Precomputation on A_x, A_y
    let mut A_xy = A_poly::new(D, d);
    
    // Slicing
    let mut a = &mut a_vec.clone()[..];
    let mut b = &mut b_vec.clone()[..];
    let mut G = &mut CP.G.clone()[..];
    let mut H = &mut CP.H.clone()[..];
    
    // Calculate Initial P_vec
    let mut P_vec: Vec<GG1Affine> = inner_commit(&G, &a).into_iter()
                                    .chain(inner_commit(&H, &b).into_iter()).collect();
    
    
    // Compute Inner Product Value
    let mut c = inner_prod(a, b);
    
    // Global Offset on the Execution Trace
    let mut offset = 0;
    
    // Protocol.Row
    for i in 0..lg_m {
        m >>=1;
        // Step 1. Split
        let (a_L, a_R) = a.split_at_mut(m * n);
        let (b_L, b_R) = b.split_at_mut(m * n);
        let (G_L, G_R) = G.split_at_mut(m);
        let (H_L, H_R) = H.split_at_mut(m);
        
        // Step 2. Calculate L, R, c_L, c_R
        let L_vec: Vec<GG1Affine> = inner_commit(&G_R, &a_L).into_iter()
                                    .chain(inner_commit(&H_L, &b_R).into_iter()).collect();
        let R_vec: Vec<GG1Affine> = inner_commit(&G_L, &a_R).into_iter()
                                    .chain(inner_commit(&H_R, &b_L).into_iter()).collect();
        let c_L = inner_prod(&a_L, &b_R);
        let c_R = inner_prod(&a_R, &b_L);
                        
        let (L_x, L_y) = Leopard_with_encoding(
            &eval_xi, &ck_pc, &L_vec, D, A_xy.idx, 3 * d
        );
        let (R_x, R_y) = Leopard_with_encoding(
            &eval_xi, &ck_pc, &R_vec, D, A_xy.idx + d, 3 * d
        );
        let (P_x, P_y) = Leopard_with_encoding(
            &eval_xi, &ck_pc, &P_vec, D, A_xy.idx + 2 * d, 3 * d
        );
        
        // Step 3, Send L, R and Get Challenge
        transcript.write_base(c_L);
        transcript.write_base(c_R);
        transcript.write_gt(L_x);
        transcript.write_gt(L_y);
        transcript.write_gt(R_x);
        transcript.write_gt(R_y);
        transcript.write_gt(P_x);
        transcript.write_gt(P_y);
        
        A_xy.register_three(&L_vec, &R_vec, &P_vec);
        
        let x: GScalar = *transcript.squeeze_challenge_base::<()>();
        let x_inv = x.invert().unwrap();
        
        // Step 4. Compute P
        let (P_vec, diff) = op_PTROW(&L_vec, &R_vec, &P_vec, x, x_inv, &mut ext, &mut perm, offset);
        offset += diff;
        
        // Step 5. Update
        for i in 0..m {
            G_L[i] = (G_L[i] + G_R[i] * x_inv).into();
            H_L[i] = (H_L[i] + H_R[i] * x).into();
        }
        
        for i in 0..(m*n) {
            a_L[i] = a_L[i] + a_R[i] * x;
            b_L[i] = b_L[i] + b_R[i] * x_inv;
        }
        c = c_L * x_inv + c + c_R * x_inv;
        a = a_L;
        b = b_L;
        G = G_L;
        H = H_L;
    }
    
    
    ///////////////////////////////
    // HERE COMES THE INTERLUDE! //
    ///////////////////////////////
    
    // Compute the commitment of the last P_vec
    let (P_x, P_y) = Leopard_with_encoding(&eval_xi, &ck_pc, &P_vec, D, A_xy.idx, d);
    transcript.write_gt(P_x);
    transcript.write_gt(P_y);
    
    // Manually Store P_vec to Execution Trace & A_xy!
    // A_xy? Ok...!
    A_xy.register(&P_vec);
    
    // Execution Trace; With a Special Care :)
    let n_slots = P_vec.len() - 1;
    let zero = BScalar::ZERO; let one = BScalar::ONE;
    
    for i in 0..n_slots {
        for j in 0..d {
            ext[0][offset + j] = zero; ext[1][offset + j] = zero; 
            ext[2][offset + j] = P_vec[i].x; ext[3][offset + j] = P_vec[i].y; 
            ext[4][offset + j] = one; ext[5][offset + j] = P_vec[i].x;
            ext[6][offset + j] = P_vec[i].y;
        }
        offset += d;
    }
    // Last: Only a single element!
    ext[0][offset] = zero; ext[1][offset] = zero;
    ext[2][offset] = P_vec[n_slots].x; ext[3][offset] = P_vec[n_slots].y; 
    ext[4][offset] = one; ext[5][offset] = P_vec[n_slots].x;
    ext[6][offset] = P_vec[n_slots].y;
    offset += 1;
    
    // Protocol.Col
    for j in 0..lg_n {
        n>>=1;
        
        // Step 1. Split
        let (a_L, a_R) = a.split_at_mut(n);
        let (b_L, b_R) = b.split_at_mut(n);
        let (P_L, P_R) = P_vec.split_at(2 * n);
        let (P1, P2) = P_L.split_at(n);
        let (P3, P4) = P_R.split_at(n);
        
        // Step 2, Calculate C_L, C_R
        let c_L = inner_prod(&a_L, &b_R);
        let c_R = inner_prod(&a_R, &b_L);
        
        // Step 3. Send them and Get challenge;
        transcript.write_base(c_L);
        transcript.write_base(c_R);
        
        let x: GScalar = *transcript.squeeze_challenge_base::<()>();
        let x_inv = x.invert().unwrap();
        
        // Step 4. Compute P_hat and Send it.
        let (P_vec, diff) = op_PTCOL(P1, P2, P3, P4, x, x_inv, &mut ext, offset);
        offset += diff;        
        let (P_x, P_y) = Leopard_with_encoding(&eval_xi, &ck_pc, &P_vec, D, A_xy.idx, d);
        A_xy.register(&P_vec);
        
        transcript.write_gt(P_x);
        transcript.write_gt(P_y);
 
        // Step 5. Update
        for i in 0..a_L.len() {
            a_L[i] = a_L[i] + a_R[i] * x;
            b_L[i] = b_L[i] + b_R[i] * x_inv;
        }
        a = a_L;
        b = b_L;
        c = c_L * x_inv + c + c_R * x;
    }
    
    
    // Last Phase
    transcript.write_base(a[0]);
    transcript.write_base(b[0]);
    let P_pub_L = G[0] * a[0];
    let P_pub_R = H[0] * b[0];
    transcript.write_scalar(P_pub_L.x);
    transcript.write_scalar(P_pub_L.y);
    transcript.write_scalar(P_pub_R.x);
    transcript.write_scalar(P_pub_R.y);
    
    ///////////////////
    // AggMEC STARTS //
    ///////////////////
    // Setups
    // Heres for AggMEC
    // Step 1. Compute a1x, a2x, w1, w2, q1, q2
    // A1, A2 => Precomputed!
    let (A_x_vec, A_y_vec) = (A_xy.coeff_x.clone(), A_xy.coeff_y.clone());
    let (W_x_vec, W_y_vec) = (ext[5].clone(), ext[6].clone());
    
    // Compute Q1, Q2
    let A_x = eval_xi.lagrange_from_vec(A_x_vec);
    let A_y = eval_xi.lagrange_from_vec(A_y_vec);
    let W_x = eval_xi.lagrange_from_vec(W_x_vec);
    let W_y = eval_xi.lagrange_from_vec(W_y_vec);
    
    // Technique? 
    let Q_x = W_x.clone() + &(A_x.clone() * (-BScalar::ONE));
    let Q_y = W_y.clone() + &(A_y.clone() * (-BScalar::ONE));
    
    // Strategy: Embed Q_x
    let mut Q_x = eval_xi.lagrange_to_coeff(Q_x);
    let mut Q_y = eval_xi.lagrange_to_coeff(Q_y);
    eval_xi.distribute_powers_zeta(&mut Q_x[..], true);
    eval_xi.distribute_powers_zeta(&mut Q_y[..], true);
    let Q_x = eval_xi.coeff_to_lagrange(Q_x);
    let Q_y = eval_xi.coeff_to_lagrange(Q_y);
    let Q_x = eval_zeta.extended_from_vec(Q_x[..].to_vec());
    let Q_y = eval_zeta.extended_from_vec(Q_y[..].to_vec());
    let Q_x = eval_zeta.divide_by_vanishing_poly(Q_x);
    let Q_y = eval_zeta.divide_by_vanishing_poly(Q_y);
    let coeff_Q_x = eval_zeta.extended_to_coeff_wo_trunc(Q_x);
    let coeff_Q_y = eval_zeta.extended_to_coeff_wo_trunc(Q_y);
    let coeff_Q_x = eval_xi.coeff_from_vec(coeff_Q_x);
    let coeff_Q_y = eval_xi.coeff_from_vec(coeff_Q_y);
    let Q_x = eval_xi.coeff_to_lagrange(coeff_Q_x.clone());
    let Q_y = eval_xi.coeff_to_lagrange(coeff_Q_y.clone());    
    
    // Do iFFT 
    let coeff_W_x = eval_xi.lagrange_to_coeff(W_x);
    let coeff_W_y = eval_xi.lagrange_to_coeff(W_y);    

    transcript.write_gt(Leopard_commit(&ck_pc, &coeff_W_x[..]));
    transcript.write_gt(Leopard_commit(&ck_pc, &coeff_W_y[..]));
    transcript.write_gt(Leopard_commit(&ck_pc, &coeff_Q_x[..]));
    transcript.write_gt(Leopard_commit(&ck_pc, &coeff_Q_y[..]));
    
    // Step 2. Get challenges
    let z: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    let rho: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    
    // Step 3. Compute Fv 
    // TODO: MultiThread
    let F_V = {
        let n_stores = A_xy.count();
        let rho_big = rho.pow_vartime([n_stores as u64, 0, 0, 0]);
        let mut acc = BScalar::ONE;
        
        let (mut A_x_vec, mut A_y_vec) = (A_xy.coeff_x, A_xy.coeff_y);
        for i in 0..n_stores {
            let (start, stepsize, count) = A_xy.record[i];
            let mut cursor = start;
            
            for j in 0..count {
                A_x_vec[cursor] *= acc;
                A_y_vec[cursor] *= acc;
                
                cursor += stepsize;
            }
            acc *= rho;
        }
        let A_x = eval_xi.lagrange_from_vec(A_x_vec);
        let A_y = eval_xi.lagrange_from_vec(A_y_vec);
        A_x + &(A_y * rho_big)
    };
    
    // Do iFFT
    let coeff_F_V = eval_xi.lagrange_to_coeff(F_V.clone());
    
    // Compute s, t1, t2, r1, r2 by evaluating z
    let s  = eval_as_poly(&coeff_F_V[..], z);
    let t1 = eval_as_poly(&coeff_Q_x[..], z);
    let t2 = eval_as_poly(&coeff_Q_y[..], z);    
    let r1 = eval_as_poly(&coeff_W_x[..], z);
    let r2 = eval_as_poly(&coeff_W_y[..], z);

    
    transcript.write_scalar(s);
    transcript.write_scalar(t1);
    transcript.write_scalar(t2);
    transcript.write_scalar(r1);
    transcript.write_scalar(r2);
    
    // Step 4. Get tau and compute F_P
    let tau: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    let tau_sq = tau.square();
    let F_P = F_V + &(A_x * tau) + &(A_y * tau_sq) 
                   + &(Q_x * (tau * tau_sq)) + &(Q_y * tau_sq.square());
    // Do iFFT
    let coeff_F_P = eval_xi.lagrange_to_coeff(F_P);
    
    // Heres for Plonkish
    // Step 0. Pre-computation for ID, Permutation, and Columns
    let (idpolys, permpolys) = get_id_and_permpolys(&eval_xi, D, &perm);
    
    let mut coeff_idpolys: Vec<Polynomial<BScalar, Coeff>> = Vec::with_capacity(5);
    let mut coeff_permpolys: Vec<Polynomial<BScalar, Coeff>> = Vec::with_capacity(5);
    let mut coeff_colpolys: Vec<Polynomial<BScalar, Coeff>> = Vec::with_capacity(7);
    
    let mut ext_idpolys: Vec<Polynomial<BScalar, ExtendedLagrangeCoeff>> = Vec::with_capacity(5);
    let mut ext_permpolys: Vec<Polynomial<BScalar, ExtendedLagrangeCoeff>> = Vec::with_capacity(5);
    let mut ext_colpolys: Vec<Polynomial<BScalar, ExtendedLagrangeCoeff>> = Vec::with_capacity(7);
    
    for poly in idpolys.clone().into_iter() {
        let poly_lag = eval_xi.lagrange_from_vec(poly);
        let poly_coeff = eval_xi.lagrange_to_coeff(poly_lag);
        let poly_ext = eval_xi.coeff_to_extended(poly_coeff.clone());
        let commit = Leopard_commit(&ck_pc, &poly_coeff[..]);
        coeff_idpolys.push(poly_coeff);
        ext_idpolys.push(poly_ext);
        transcript.write_gt(commit);
    }
    
    for poly in permpolys.clone().into_iter() {
        let poly_lag = eval_xi.lagrange_from_vec(poly);
        let poly_coeff = eval_xi.lagrange_to_coeff(poly_lag);
        let poly_ext = eval_xi.coeff_to_extended(poly_coeff.clone());
        let commit = Leopard_commit(&ck_pc, &poly_coeff[..]);
        coeff_permpolys.push(poly_coeff);
        ext_permpolys.push(poly_ext);
        transcript.write_gt(commit);
    }
    
    for poly in ext.clone().into_iter() {
        let poly_lag = eval_xi.lagrange_from_vec(poly);
        let poly_coeff = eval_xi.lagrange_to_coeff(poly_lag);
        let poly_ext = eval_xi.coeff_to_extended(poly_coeff.clone());
        let commit = Leopard_commit(&ck_pc, &poly_coeff[..]);
        coeff_colpolys.push(poly_coeff);
        ext_colpolys.push(poly_ext);
        transcript.write_gt(commit);
    }    
    
    // Step 1. Setting up Permutation Polynomial
    // Get u1, u2
    let u1: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    let u2: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    
    // Construct Z polynomial and Commit it!
    let tmp_Z_poly = get_z_poly(&eval_xi, &ext, &idpolys, &permpolys, u1, u2);
    let tmp_Z_poly = eval_xi.lagrange_from_vec(tmp_Z_poly);
    let coeff_Z_poly = eval_xi.lagrange_to_coeff(tmp_Z_poly);
    let Z_poly = eval_xi.coeff_to_extended(coeff_Z_poly.clone());
    let comm_Z = Leopard_commit(&ck_pc, &coeff_Z_poly[..]);
    
//     assert_eq!(coeff_Z_poly[0], BScalar::ONE);
    
    transcript.write_gt(comm_Z);
    
    // Get u3 and eta
    let u3: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    let eta: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    
    // Step 2: Construct T & Q
    // Precompute: onepoly & H1
    let onepoly = eval_xi.constant_extended(BScalar::ONE);
    let H1 = {
        let mut tmp = vec![BScalar::ZERO; D];
        tmp[0] = BScalar::ONE;
        let tmp = eval_xi.lagrange_from_vec(tmp);
        let tmp = eval_xi.lagrange_to_coeff(tmp);
        eval_xi.coeff_to_extended(tmp)
    };        
    
    let T_left = constructor_T(&eval_xi, &ext_colpolys, &onepoly, eta);
    let (T_mid, T_right) = constructor_Z(
        &eval_xi, &Z_poly, &ext_colpolys[2..7], &ext_idpolys, &ext_permpolys, 
        &onepoly, &H1, u1, u2,
    );
    let T = T_left + &(T_mid * u3) + &(T_right * u3.square());
    
    let T_coeff = eval_xi.extended_to_coeff_wo_trunc(T.clone());
    
    let Q_coeff = {
        let Q = eval_xi.divide_by_vanishing_poly(T);
        eval_xi.extended_to_coeff(Q)
    };
        
    // Chunk T & Q; Send to Verifier
    let mut coeff_Ts: Vec<Polynomial<BScalar, Coeff>> = Vec::with_capacity(6);
    let mut coeff_Qs: Vec<Polynomial<BScalar, Coeff>> = Vec::with_capacity(5);
    
    for Ti in T_coeff[..].chunks(D) {
        let com_T = Leopard_commit(&ck_pc, Ti);
        coeff_Ts.push(eval_xi.coeff_from_vec(Ti.to_vec()));
        transcript.write_gt(com_T);
    }
    
    for Qi in Q_coeff[..].chunks(D) {
        let com_Q = Leopard_commit(&ck_pc, Qi);
        coeff_Qs.push(eval_xi.coeff_from_vec(Qi.to_vec()));
        transcript.write_gt(com_Q);
    }
       
    // Step 3: Get u4
    let u4: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    let wu4 = eval_xi.get_omega() * u4;
    
    // Step 4: Compute Evaluations
    // ts
    for poly in coeff_Ts.iter() {
        let ts = eval_as_poly(poly, u4);
        transcript.write_scalar(ts);
    }
    
    // qs
    for poly in coeff_Qs.iter() {
        let qs = eval_as_poly(poly, u4);
        transcript.write_scalar(qs);
    }
    
    // alphas
    for poly in coeff_colpolys.iter() {
        let alphas = eval_as_poly(poly, u4);
        transcript.write_scalar(alphas);
    }
    
    // beta
    let beta = eval_as_poly(&coeff_Z_poly[..], u4);
    transcript.write_scalar(beta);
    
    // beta_aa
    let beta_aa = eval_as_poly(&coeff_colpolys[2][..], wu4);
    transcript.write_scalar(beta_aa);
    let beta_aa = eval_as_poly(&coeff_colpolys[3][..], wu4);
    transcript.write_scalar(beta_aa);
    let beta_aa = eval_as_poly(&coeff_colpolys[4][..], wu4);
    transcript.write_scalar(beta_aa);
    
    // gamma
    let gamma = eval_as_poly(&coeff_Z_poly[..], wu4);
    transcript.write_scalar(gamma);
    
    // phis
    for poly in coeff_idpolys.iter() {
        let phis = eval_as_poly(&poly[..], u4);
        transcript.write_scalar(phis);
    }    
    
    // psis
    for poly in coeff_permpolys.iter() {
        let psis = eval_as_poly(&poly[..], u4);
        transcript.write_scalar(psis);
    }    
    
    
    // Step 5. Run Batched Eval
    let vecs_z: Vec<Polynomial<BScalar, Coeff>> = vec![coeff_F_P, coeff_W_x, coeff_W_y];

    let vecs_u4: Vec<Polynomial<BScalar, Coeff>> = vec![
        // Ts
        coeff_Ts[0].clone(), coeff_Ts[1].clone(), coeff_Ts[2].clone(), coeff_Ts[3].clone(),
        coeff_Ts[4].clone(), coeff_Ts[5].clone(),
        // Qs
        coeff_Qs[0].clone(), coeff_Qs[1].clone(), coeff_Qs[2].clone(), coeff_Qs[3].clone(), 
        coeff_Qs[4].clone(), 
        // alphas
        coeff_colpolys[0].clone(), coeff_colpolys[1].clone(), 
        coeff_colpolys[2].clone(), coeff_colpolys[3].clone(), coeff_colpolys[4].clone(), 
        coeff_colpolys[5].clone(), coeff_colpolys[6].clone(), 
        // betas
        coeff_Z_poly.clone(),
        // phis
        coeff_idpolys[0].clone(), coeff_idpolys[1].clone(), coeff_idpolys[2].clone(), 
        coeff_idpolys[3].clone(), coeff_idpolys[4].clone(),
        // psis
        coeff_permpolys[0].clone(), coeff_permpolys[1].clone(), coeff_permpolys[2].clone(), 
        coeff_permpolys[3].clone(), coeff_permpolys[4].clone(),
    ];
    
    // Prove on wu4
    let vecs_wu4: Vec<Polynomial<BScalar, Coeff>> = vec![
        // Referring Next Column
        coeff_colpolys[2].clone(), coeff_colpolys[3].clone(), coeff_colpolys[4].clone(),
        // Z
        coeff_Z_poly.clone()
    ];
    
    leopard_merged_prover(
        &eval_xi,
        &ck_pc,
        transcript,
        vec![vecs_z, vecs_u4, vecs_wu4],
        vec![z, u4, wu4],
    );

    // Done!
}
