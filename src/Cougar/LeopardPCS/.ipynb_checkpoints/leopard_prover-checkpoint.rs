use halo2curves::bn256::{
    G1Affine, G2Affine, G2Prepared, Gt, Fr as Scalar
};

use halo2curves::bn256::{pairing, multi_miller_loop};
use halo2curves::pairing::MillerLoopResult;


use crate::Cougar::transcript::{
    Transcript,
    TranscriptWrite,
    Blake2bWrite,
    Challenge255,
};

use crate::Cougar::utils::{
    inner_prod,
    ped_commit,
    Leopard_commit,
    eval_as_poly,
    kate_division,
};

use crate::Cougar::setup::{
    LeopardParams,
};

use halo2curves::pairing::{PairingCurveAffine, MultiMillerLoop};

use ff::{Field, PrimeField};

use halo2_proofs::poly::{
    EvaluationDomain,    
    Polynomial,
    Coeff,
};

pub fn Leopard_commit_(G: &[G1Affine], H: &[G2Affine], msg: &[Scalar]) -> Gt
{
    let n1 = G.len();
    let n2 = H.len();
    
    // 1st layer
    let mut left: Vec<G1Affine> = Vec::with_capacity(n2);
    
    for i in 0..n2 {
        let tmp_msg: Vec<Scalar> = (0..n1).into_iter().map(|j| msg[j * n2 + i]).collect();
        left.push(ped_commit(G, &tmp_msg).into());
    }
    
    let right: Vec<G2Prepared> = H.into_iter().map(|a| G2Prepared::from(*a)).collect();
    let prepared: Vec<(&G1Affine, &G2Prepared)> = left.iter().zip(right.iter()).collect();
    multi_miller_loop(&prepared).final_exponentiation()
}


pub fn leopard_prover(
    ck: &LeopardParams,
    transcript: &mut Blake2bWrite<Vec<u8>, G1Affine, Challenge255<G1Affine>>,
    a_vec: Vec<Scalar>,
    point: Scalar,
)
{
    let mut m = ck.G.len();
    let mut n = ck.H.len();
    
    // Blind Term for U
    let r = *transcript.squeeze_challenge_scalar::<()>();
    let U = ck.U * r;
    
    let eval = eval_as_poly(&a_vec, point);
    let mut z_vec = {
        let mut tmp: Vec<Scalar> = Vec::new();
        let mut acc = Scalar::ONE;
        for i in 0..(m*n) {
            tmp.push(acc);
            acc *= point;
        }
        tmp
    };
    
    let mut a = &mut a_vec.clone()[..];
    let mut G = &mut ck.G.clone()[..];
    let mut H = &mut ck.H.clone()[..];
    let mut z = &mut z_vec[..];

    // PT3_ROW
    while (m > 1) {
        m >>= 1;
        let (a_L, a_R) = a.split_at_mut(m * n);
        let (z_L, z_R) = z.split_at_mut(m * n);
        let (G_L, G_R) = G.split_at_mut(m);
        
        let c_L = inner_prod(a_L, z_R);
        let c_R = inner_prod(a_R, z_L);
                
        let L = Leopard_commit_(G_R, H, a_L) + U * c_L;
        let R = Leopard_commit_(G_L, H, a_R) + U * c_R;
        
        transcript.write_gt(L);
        transcript.write_gt(R);
        
        let x = *transcript.squeeze_challenge_scalar::<()>();
        let x_inv = x.invert().unwrap();
        
        for i in 0..m*n {
            a_L[i] = a_L[i] + a_R[i] * x_inv;
            z_L[i] = z_L[i] + z_R[i] * x;
        }
        
        for i in 0..m {
            G_L[i] = (G_L[i] + G_R[i] * x).into();
        }
        
        a = a_L;
        z = z_L;
        G = G_L;
    }
    
    // PT3_COL
    while (n > 1) {
        n >>=1;
        let (a_L, a_R) = a.split_at_mut(m * n);
        let (z_L, z_R) = z.split_at_mut(m * n);
        let (H_L, H_R) = H.split_at_mut(n);
        
        let c_L = inner_prod(a_L, z_R);
        let c_R = inner_prod(a_R, z_L);
        let L = Leopard_commit_(G, H_R, a_L) + U * c_L;
        let R = Leopard_commit_(G, H_L, a_R) + U * c_R;
        
        transcript.write_gt(L);
        transcript.write_gt(R);
        
        let x = *transcript.squeeze_challenge_scalar::<()>();
        let x_inv = x.invert().unwrap();
        
        for i in 0..n {
            a_L[i] = a_L[i] + a_R[i] * x_inv;
            z_L[i] = z_L[i] + z_R[i] * x;
            H_L[i] = (H_L[i] + H_R[i] * x).into();
        }
        
        a = a_L;
        z = z_L;
        H = H_L;
    }
    
    // Final
    transcript.write_scalar(a[0]);
    transcript.write_scalar(z[0]);
}

// Single Point, Many Polynomials
pub fn leopard_batch_prover(
    ck: &LeopardParams,
    transcript: &mut Blake2bWrite<Vec<u8>, G1Affine, Challenge255<G1Affine>>,
    a_vecs: Vec<Polynomial<Scalar, Coeff>>,
    point: Scalar,
)
{
    // Aggregation
    let gamma = *transcript.squeeze_challenge_scalar::<()>();
    let n_vecs = a_vecs.len();
    let N = a_vecs[0].len();
    
    let agg_poly = a_vecs.iter()
                          .rev()
                          .skip(1)
                          .fold(a_vecs[n_vecs-1].clone(), |acc, x| acc * gamma + x);
    
    let _ = leopard_prover(ck, transcript, (&agg_poly[..]).to_vec(), point);
}


pub fn leopard_merged_prover(
    D: &EvaluationDomain<Scalar>,
    ck: &LeopardParams,
    transcript: &mut Blake2bWrite<Vec<u8>, G1Affine, Challenge255<G1Affine>>,
    a_vecs: Vec<Vec<Polynomial<Scalar, Coeff>>>,
    points: Vec<Scalar>,
)
{
    let mut agg_polys: Vec<Polynomial<Scalar, Coeff>> = Vec::new();
    
    let ext_points = {
        let mut acc: Vec<Scalar> = Vec::new();
        for i in 0..points.len() {
            for j in 0..a_vecs[i].len() {
                acc.push(points[i]);
            };              
        }
        acc
    };
    for i in 0..a_vecs.len() {
        for j in 0..a_vecs[i].len() {
            agg_polys.push(a_vecs[i][j].clone());
        }
    }
    
    
    let mut q_polys = Vec::with_capacity(agg_polys.len());
    for i in 0..agg_polys.len() {
        let q_poly = kate_division(&agg_polys[i], ext_points[i]);
        q_polys.push(D.coeff_from_vec(q_poly));
    }    
    // compute h(X, z) and commit H
    let z = *transcript.squeeze_challenge_scalar::<()>();
    let h_poly = q_polys.iter().rev().skip(1).fold(q_polys[q_polys.len()-1].clone(), |acc, x| acc * z + x);
    transcript.write_gt(Leopard_commit(ck, &h_poly));

    // Get a random challenge r
    let r = *transcript.squeeze_challenge_scalar::<()>();
    agg_polys.push(h_poly);
    // Compute f(r) & h(r, z) and send it!
    for i in 0..agg_polys.len() {
        transcript.write_scalar(eval_as_poly(&agg_polys[i], r));
    }
    
    
    // Open agg_polys & H at r
    leopard_batch_prover(
        ck,
        transcript,
        agg_polys,
        r,
    );
}