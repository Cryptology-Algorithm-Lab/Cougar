use ff::{Field, PrimeField};
use group::Group;
use halo2curves::pairing::{PairingCurveAffine, MultiMillerLoop};
use halo2curves::msm::best_multiexp;
use halo2curves::CurveAffine;

// Abstraction
pub fn inner_prod<T: Field>(a: &[T], b:&[T]) -> T {
    a.into_iter().zip(b.into_iter()).map(|(x,y)| *x*y).sum()
}

pub fn ped_commit<C: CurveAffine>(ck: &[C], msg: &[C::Scalar]) -> C::Curve {
    best_multiexp(msg, ck)
}

// Inner Commitment
pub fn inner_commit<C: CurveAffine>(ck: &[C], msg:&[C::Scalar]) -> Vec<C> {
    let n1 = ck.len();
    let n2 = msg.len() / n1;
    
    let mut ret: Vec<C> = Vec::with_capacity(n2);
    
    for i in 0..n2 {
        let tmp_msg: Vec<C::Scalar> = (0..n1).into_iter().map(|j| msg[j * n2 + i]).collect();
        ret.push(ped_commit(ck, &tmp_msg).into());
    }
    ret
}

pub fn eval_as_poly<F: PrimeField>(poly: &[F], point: F) -> F
{
    poly.iter().rev().fold(F::ZERO, |acc, x| acc * point + x)    
}

// Instantiated by Grumkin-BN254 Cycle

use halo2curves::bn256::{G1Affine, G2Affine, G2Prepared, Gt, Fr as BScalar};
use halo2curves::grumpkin::{G1Affine as GG1Affine, Fr as GScalar, Fq as GBase};
use halo2curves::bn256::{pairing, multi_miller_loop};
use halo2curves::pairing::MillerLoopResult;


use crate::Cougar::setup::{LeopardParams, CougarParams};

// Leopard PCS.commit
pub fn Leopard_commit(ck_pc: &LeopardParams, msg: &[BScalar]) -> Gt
{
    let n1 = ck_pc.G.len();
    let n2 = ck_pc.H.len();
    
    // 1st layer
    let mut left: Vec<G1Affine> = Vec::with_capacity(n2);
    
    for i in 0..n2 {
        let tmp_msg: Vec<BScalar> = (0..n1).into_iter().map(|j| msg[j * n2 + i]).collect();
        left.push(ped_commit(&ck_pc.G, &tmp_msg).into());
    }
    
    let right: Vec<G2Prepared> = ck_pc.H.clone().into_iter().map(|a| G2Prepared::from(a)).collect();
    let prepared: Vec<(&G1Affine, &G2Prepared)> = left.iter().zip(right.iter()).collect();
    multi_miller_loop(&prepared).final_exponentiation()
}

// Kate Division
pub fn kate_division<F:Field>(a: &[F], b: F) -> Vec<F>
{
    let b = -b;

    let mut q = vec![F::ZERO; a.len() - 1];

    let mut tmp = F::ZERO;
    for (q, r) in q.iter_mut().rev().zip(a.iter().rev()) {
        let mut lead_coeff = *r;
        lead_coeff.sub_assign(&tmp);
        *q = lead_coeff;
        tmp = lead_coeff;
        tmp.mul_assign(&b);
    }
    
    q.resize(a.len(), F::ZERO);
    q
}