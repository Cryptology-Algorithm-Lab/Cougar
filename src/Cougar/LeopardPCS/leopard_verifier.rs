use halo2curves::bn256::{
    G1Affine, G2Affine, G2Prepared, Gt, Fr as Scalar
};

use halo2curves::bn256::{pairing, multi_miller_loop};
use halo2curves::pairing::MillerLoopResult;
use halo2curves::msm::best_multiexp;
use halo2curves::msm_gt::{best_multiexp as best_multiexp_gt};

use halo2_proofs::poly::{
    Polynomial,
    Coeff,
};



use crate::Cougar::transcript::{
    Transcript,
    TranscriptRead,
    Blake2bRead,
    Challenge255,
};

use crate::Cougar::utils::{
    inner_prod,
    ped_commit
};

use crate::Cougar::setup::{
    LeopardParams,
};

use halo2curves::pairing::{PairingCurveAffine, MultiMillerLoop};

use ff::{Field, PrimeField, BatchInvert};

pub fn leopard_verifier(
    ck: &LeopardParams,
    transcript: &mut Blake2bRead<&[u8], G1Affine, Challenge255<G1Affine>>,
    merged_commit_base: Vec<Gt>,
    merged_commit_expon: Vec<Scalar>,
    point: Scalar,
    eval: Scalar,
)
{
    let m = ck.G.len(); let n = ck.H.len();
    let lgm = m.next_power_of_two().trailing_zeros() as usize;
    let lgn = n.next_power_of_two().trailing_zeros() as usize;
    
    let mut L_vec: Vec<Gt> = Vec::with_capacity(lgm + lgn);
    let mut R_vec: Vec<Gt> = Vec::with_capacity(lgm + lgn);
    let mut x_vec = Vec::with_capacity(lgm + lgn);
    
    let r = *transcript.squeeze_challenge_scalar::<()>();
    let U = ck.U * r;
    
    
    for i in 0..(lgm + lgn) {
        L_vec.push(transcript.read_gt().unwrap());
        R_vec.push(transcript.read_gt().unwrap());
        x_vec.push(*transcript.squeeze_challenge_scalar::<()>());
    }
    
    let mut x_inv_vec = x_vec.clone();
    x_inv_vec.batch_invert();
    
    
    let a0 = transcript.read_scalar().unwrap();
    let z0 = transcript.read_scalar().unwrap();
    
    let mut s_L = Vec::with_capacity(m);
    let mut s_R = Vec::with_capacity(n);
    
    s_L.push(Scalar::ONE); s_R.push(Scalar::ONE);
    
    for i in 1..m {
        let lgi: usize = (32 - 1 - (i as u32).leading_zeros()) as usize;
        let k = 1 << lgi;
        let u = x_vec[lgm - 1 - lgi];
        s_L.push(s_L[i-k] * u);
    }

    for i in 1..n {
        let lgi: usize = (32 - 1 - (i as u32).leading_zeros()) as usize;
        let k = 1 << lgi;
        let u = x_vec[(lgm + lgn - 1) - lgi];
        s_R.push(s_R[i-k] * u);
    }
    
    let left = {
        let z_giant = point.pow([n as u64]);
        let mut z1: Vec<Scalar> = Vec::with_capacity(m);
        let mut z2: Vec<Scalar> = Vec::with_capacity(n);
        
        let mut z1_acc = Scalar::ONE;
        let mut z2_acc = Scalar::ONE;
        
        for i in 0..m {
            z1.push(z1_acc);
            z1_acc *= z_giant;
        }
        
        for i in 0..n {
            z2.push(z2_acc);
            z2_acc *= point;
        }
        
        let z_final = inner_prod(&s_L, &z1) * inner_prod(&s_R, &z2);
        let expon: Vec<Scalar> = s_L.clone().into_iter().map(|x| a0 * x).collect();
        let pair_l: G1Affine = best_multiexp(&expon, &ck.G).into();
        let pair_r: G2Affine = best_multiexp(&s_R, &ck.H).into();
        pairing(&pair_l, &pair_r) + U * (a0 * z_final)
            
    };

    let right = best_multiexp_gt(
        &merged_commit_expon.into_iter()
                           .chain(x_vec.into_iter()).chain(x_inv_vec.into_iter())
                           .chain(Some(eval).into_iter()).collect::<Vec<Scalar>>(),
        &merged_commit_base.into_iter()
                           .chain(L_vec.into_iter()).chain(R_vec.into_iter())
                           .chain(Some(U).into_iter()).collect::<Vec<Gt>>(),
    );
    
    assert_eq!(left, right);
}

pub fn leopard_batch_verifier(
    ck: &LeopardParams,
    transcript: &mut Blake2bRead<&[u8], G1Affine, Challenge255<G1Affine>>,
    LRP_x: Vec<Gt>, LRP_y: Vec<Gt>, Qx: Gt, Qy: Gt,
    expon: Vec<Scalar>,    
    commits: Vec<Gt>,
    point: Scalar,
    evals: Vec<Scalar>,
)
{
    let gamma = *transcript.squeeze_challenge_scalar::<()>();
    let gammas = {
        let mut acc = gamma;
        let mut ret: Vec<Scalar> = Vec::with_capacity(commits.len());
        for i in 0..commits.len() {
            ret.push(acc);
            acc *= gamma;
        }
        ret
    };

    let merged_commit_base = LRP_x.into_iter().chain(LRP_y.into_iter())
        .chain(Some(Qx).into_iter()).chain(Some(Qy).into_iter()).chain(commits.into_iter()).collect::<Vec<Gt>>();
    let merged_commit_expon = expon.into_iter().chain(gammas.into_iter()).collect::<Vec<Scalar>>();

    let expected_eval = evals.iter().rev().fold(Scalar::ZERO, |acc, x| acc * gamma + x);
    
    leopard_verifier(
        ck,
        transcript,
        merged_commit_base,
        merged_commit_expon,
        point,
        expected_eval,
    )
}


pub fn leopard_merged_verifier(
    ck: &LeopardParams,
    transcript: &mut Blake2bRead<&[u8], G1Affine, Challenge255<G1Affine>>,
    LRP_x: Vec<Gt>, LRP_y: Vec<Gt>, Qx: Gt, Qy: Gt,
    expon: Vec<Scalar>,
    mut commits: Vec<Gt>,
    points: Vec<Scalar>,
    evals: Vec<Scalar>,
) 
{  
    // Get challenge z, H, r, and Hrz
    let n_points = points.len();
    let z = *transcript.squeeze_challenge_scalar::<()>();
    let H = transcript.read_gt().unwrap();
    let r = *transcript.squeeze_challenge_scalar::<()>();
    commits.push(H);
    
    let mut eval_at_r: Vec<Scalar> = Vec::new();
    for i in 0..(n_points + 1) {
        eval_at_r.push(transcript.read_scalar().unwrap());
    }
    
    // Open merged_commits at a random point
    leopard_batch_verifier(
        ck,
        transcript,
        LRP_x, LRP_y, Qx, Qy,
        expon,
        commits,
        r,
        eval_at_r.clone(),
    );


    // Construct Verification Equation
    let right = eval_at_r.iter().rev().skip(1)
                .zip(points.iter().rev())
                .zip(evals.iter().rev())
                .map(|((evr, exp), eva)| (evr - eva) * ((r - exp).invert().unwrap()))
                .fold(Scalar::ZERO, |acc, x| acc * z + x);
    assert_eq!(eval_at_r[n_points], right);
}
