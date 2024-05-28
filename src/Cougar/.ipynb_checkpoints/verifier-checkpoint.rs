use halo2curves::bn256::{
    G1Affine as BG1Affine, Gt, Fr as BScalar, Fq as GScalar,
};

use halo2curves::grumpkin::{
    G1Affine as GG1Affine
};

use halo2curves::msm::{
    best_multiexp,
};

use halo2curves::msm_gt::{best_multiexp as best_multiexp_gt};

use crate::Cougar::LeopardPCS::leopard_verifier::{
    leopard_merged_verifier,
};

use crate::Cougar::transcript::{
    Transcript,
    TranscriptRead,
    Blake2bRead,
    Challenge255,
};

use crate::Cougar::setup::{
    CougarParams,
};

use crate::Cougar::utils::{
    inner_prod,
};

use group::ff::{
    Field, PrimeField, BatchInvert
};

use crate::Cougar::custom_gates::{
    poly_eval_T,
    poly_eval_Z,
};

use halo2_proofs::poly::EvaluationDomain;


pub fn cougar_verifier(
    CP: &CougarParams,
    transcript: &mut Blake2bRead<&[u8], BG1Affine, Challenge255<BG1Affine>>,
    c: GScalar,
) 
{
    let m = CP.m; let n = CP.n; let d = CP.d; let D = CP.D;
    let lgm = m.next_power_of_two().trailing_zeros() as usize;
    let lgn = n.next_power_of_two().trailing_zeros() as usize;
    let eval_xi: &EvaluationDomain<BScalar> = &CP.eval_xi;

    let ck_pc = &CP.ck_pc;
    
    let mut LRP_x_vec : Vec<Gt> = Vec::new();
    let mut LRP_y_vec : Vec<Gt> = Vec::new();
    
    let mut cl_vec: Vec<GScalar> = Vec::with_capacity(lgm + lgn);
    let mut cr_vec: Vec<GScalar> = Vec::with_capacity(lgm + lgn);
    let mut x_vec: Vec<GScalar> = Vec::with_capacity(lgm + lgn);
        
    // Transcripts from PTROW
    for i in 0..lgm {
        // Receive cL, cR, L, R
        cl_vec.push(transcript.read_base().unwrap());
        cr_vec.push(transcript.read_base().unwrap());
        // L
        LRP_x_vec.push(transcript.read_gt().unwrap());
        LRP_y_vec.push(transcript.read_gt().unwrap());
        // R
        LRP_x_vec.push(transcript.read_gt().unwrap());
        LRP_y_vec.push(transcript.read_gt().unwrap());
        // P
        LRP_x_vec.push(transcript.read_gt().unwrap());
        LRP_y_vec.push(transcript.read_gt().unwrap());        
        
        // Get Challenge
        x_vec.push(*transcript.squeeze_challenge_base::<()>());
    }
    
    // Receive Last P
    LRP_x_vec.push(transcript.read_gt().unwrap());
    LRP_y_vec.push(transcript.read_gt().unwrap());
    
    let mut s_vec = Vec::with_capacity(m);
    s_vec.push(GScalar::ONE);
    
    for i in 1..m {
        let lgi: usize = (32 - 1 - (i as u32).leading_zeros()) as usize;
        let k = 1 << lgi;
        let u = x_vec[lgm - 1 - lgi];
        s_vec.push(s_vec[i-k] * u);
    }    
    
    // Transcripts for PTCOL
    for i in 0..lgn {
        // Receive cL, cR
        cl_vec.push(transcript.read_base().unwrap());
        cr_vec.push(transcript.read_base().unwrap());
        
        // Get Challenge
        x_vec.push(*transcript.squeeze_challenge_base::<()>());
        
        // Receive P
        LRP_x_vec.push(transcript.read_gt().unwrap());
        LRP_y_vec.push(transcript.read_gt().unwrap());
    }
        
    // Final Phase
    let a0 = transcript.read_base().unwrap();
    let b0 = transcript.read_base().unwrap();
    let Lx0 = transcript.read_scalar().unwrap();
    let Ly0 = transcript.read_scalar().unwrap();
    let Rx0 = transcript.read_scalar().unwrap();
    let Ry0 = transcript.read_scalar().unwrap();
    
    // Check c = a * b;
    let mut x_inv_vec = x_vec.clone();    
    let mut s_inv_vec = s_vec.clone();
    
    x_inv_vec.iter_mut()
            .chain(s_inv_vec.iter_mut())
            .batch_invert();
    
    let c_left = inner_prod(&x_inv_vec, &cl_vec) + c + inner_prod(&x_vec, &cr_vec);
    let c_right = a0 * b0;    
    let G_base: GG1Affine = best_multiexp(&s_inv_vec, &CP.G).into();
    let H_base: GG1Affine = best_multiexp(&s_vec, &CP.H).into();
    let P_pub_left = G_base * a0;
    let P_pub_right = H_base * b0;

    assert_eq!(c_left, c_right);
    assert_eq!(P_pub_left.x, Lx0);
    assert_eq!(P_pub_left.y, Ly0);
    assert_eq!(P_pub_right.x, Rx0);
    assert_eq!(P_pub_right.y, Ry0);

    // Step 1. Receive W, Q
    let Wx = transcript.read_gt().unwrap();
    let Wy = transcript.read_gt().unwrap();
    let Qx = transcript.read_gt().unwrap();
    let Qy = transcript.read_gt().unwrap();
    
    // Step 2. Get z, rho and compute V
    let z: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    let rho: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    
    // Step 3. Receive s, t1, t2, r1, r2
    let s: BScalar = transcript.read_scalar().unwrap();
    let t1: BScalar = transcript.read_scalar().unwrap();
    let t2: BScalar = transcript.read_scalar().unwrap();    
    let r1: BScalar = transcript.read_scalar().unwrap();
    let r2: BScalar = transcript.read_scalar().unwrap();
    
    // Step 4. Get tau and compute P & y
    let tau: BScalar = *transcript.squeeze_challenge_scalar::<()>();

    // Transcripts for Plonkish
    let mut comm_IDs: Vec<Gt> = Vec::with_capacity(5);
    let mut comm_Rs: Vec<Gt> = Vec::with_capacity(5);
    
    // Step 0. Receive ID, R, COL polys
    for i in 0..5 {
        comm_IDs.push(transcript.read_gt().unwrap());
    }
    
    for i in 0..5 {
        comm_Rs.push(transcript.read_gt().unwrap());
    }
    
    let mut comm_COLs: Vec<Gt> = Vec::with_capacity(7);
    for i in 0..7 {
        comm_COLs.push(transcript.read_gt().unwrap());
    }
    
    // Step 1. Get u1, u2, u3 and eta
    let u1: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    let u2: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    
    let Z: Gt = transcript.read_gt().unwrap();
    
    let u3: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    let eta: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    
    // Step 2. Receive T, Q polys
    let mut comm_Ts: Vec<Gt> = Vec::with_capacity(6);
    let mut comm_Qs: Vec<Gt> = Vec::with_capacity(5);
    
    for i in 0..6 {
        comm_Ts.push(transcript.read_gt().unwrap());
    }
    
    for i in 0..5 {
        comm_Qs.push(transcript.read_gt().unwrap());
    }
    
    // Step 3. Get u4
    let u4: BScalar = *transcript.squeeze_challenge_scalar::<()>();
    let wu4 = CP.eval_xi.get_omega() * u4;
    
    // Step 4. Get Eval Points
    let mut ts: Vec<BScalar> = Vec::with_capacity(6);
    for i in 0..6 {
        ts.push(transcript.read_scalar().unwrap());
    }
    
    let mut qs: Vec<BScalar> = Vec::with_capacity(5);
    for i in 0..5 {
        qs.push(transcript.read_scalar().unwrap());
    }
    
    let mut alphas: Vec<BScalar> = Vec::with_capacity(7);
    for i in 0..7 {
        alphas.push(transcript.read_scalar().unwrap());
    }
    
    let beta: BScalar = transcript.read_scalar().unwrap();
    
    let mut  beta_aa: Vec<BScalar> = Vec::with_capacity(3);
    for i in 0..3 {
        beta_aa.push(transcript.read_scalar().unwrap());
    }
    
    let gamma: BScalar = transcript.read_scalar().unwrap();
    
    let mut phis: Vec<BScalar> = Vec::with_capacity(5);
    for i in 0..5 {
        phis.push(transcript.read_scalar().unwrap());
    }    
    
    let mut psis: Vec<BScalar> = Vec::with_capacity(5);
    for i in 0..5 {
        psis.push(transcript.read_scalar().unwrap());
    }
    
    // Step 5. Construct and Check Eval Values
    let u4D = u4.pow_vartime([D as u64, 0, 0, 0]);
    let H1_eval = {
        let mut denom = BScalar::from(D as u64) * (u4 - BScalar::ONE);
        denom = denom.invert().unwrap();
        let ret = u4D - BScalar::ONE;
        ret * denom
    };
    
    let rho1 = poly_eval_T(&alphas, &beta_aa, eta) 
             + poly_eval_Z(beta, gamma, &alphas[2..7], &phis, &psis, u1, u2, u3, H1_eval);
    let rho2 = rho1 * ((u4D - BScalar::ONE).invert().unwrap());
    
    // Check Consistency with Truncated T, Qs
    let t_right = ts.iter().rev().fold(BScalar::ZERO, |acc, x| acc * u4D + x);
    let q_right = qs.iter().rev().fold(BScalar::ZERO, |acc, x| acc * u4D + x);
    
    assert_eq!(rho1, t_right);
    assert_eq!(rho2, q_right);

    let zd = z.pow_vartime([(D/d) as u64, 0, 0, 0]);
    let tau_sq = tau.square();

    // Optimization on computing P
    let expon: Vec<BScalar> = {
        let n_LRP = LRP_x_vec.len();
        let mut ret: Vec<BScalar> = Vec::with_capacity(2 * n_LRP + 2);
        let mut acc = BScalar::ONE;
        for i in 0..n_LRP {
            ret.push(acc + tau);
            acc *= rho;
        }
        for i in 0..n_LRP {
            ret.push(acc + tau_sq);
            acc *= rho;
        }
        ret.push(tau_sq * tau);
        ret.push(tau_sq.square());
        ret
    };

    let y = {
        let t1t2tau= t1 + t2 * tau;
        let r1r2tau = r1 + r2 * tau;
        let t2mzhat = tau_sq - zd + BScalar::ONE;
        s + tau * (t1t2tau * t2mzhat + r1r2tau)
    };
    
    // Step 6. Batched Eval Verification
    let mut commits: Vec<Gt> = vec![
        // AggMEC
        Wx, Wy,
        // Check on Ts
        comm_Ts[0], comm_Ts[1], comm_Ts[2], comm_Ts[3], comm_Ts[4], comm_Ts[5],
        // Check on Qs
        comm_Qs[0], comm_Qs[1], comm_Qs[2], comm_Qs[3], comm_Qs[4],
        // Check on alphas
        comm_COLs[0], comm_COLs[1], comm_COLs[2], comm_COLs[3], comm_COLs[4], comm_COLs[5],           
        comm_COLs[6],
        // Check on beta
        Z,
        // Check on phis
        comm_IDs[0], comm_IDs[1], comm_IDs[2], comm_IDs[3], comm_IDs[4],
        // Check on psis
        comm_Rs[0], comm_Rs[1], comm_Rs[2], comm_Rs[3], comm_Rs[4],
        // Open @ wu4
        comm_COLs[2], comm_COLs[3], comm_COLs[4], Z
    ];
    let evals: Vec<BScalar> = vec![
        // AggMEC
        y, r1, r2,
        // Check on Ts
        ts[0], ts[1], ts[2], ts[3], ts[4], ts[5],
        // Check on Qs
        qs[0], qs[1], qs[2], qs[3], qs[4],
        // Check on alphas
        alphas[0], alphas[1], alphas[2], alphas[3], alphas[4], alphas[5], alphas[6],
        // Check on beta
        beta,
        // Check on phis
        phis[0], phis[1], phis[2], phis[3], phis[4],
        // Check on psis
        psis[0], psis[1], psis[2], psis[3], psis[4],        
        // Open @ wu4
        beta_aa[0], beta_aa[1], beta_aa[2], gamma
    ];
    let points: Vec<BScalar> = vec![
        z, z, z,
        u4, u4, u4, u4, u4, u4,
        u4, u4, u4, u4, u4,
        u4, u4, u4, u4, u4 ,u4 ,u4,
        u4,
        u4, u4, u4, u4, u4,
        u4, u4, u4, u4, u4,
        wu4, wu4, wu4, wu4
    ];
    
    leopard_merged_verifier(
        &ck_pc,
        transcript,
        LRP_x_vec, LRP_y_vec, Qx, Qy, expon,
        commits,
        points,
        evals,
    );
}
