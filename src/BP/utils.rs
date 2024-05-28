use halo2curves::secp256k1::{
    Secp256k1Affine as G1Affine, Fq as Scalar,
};

use ff::{Field, PrimeField};

use halo2curves::msm::best_multiexp;

pub fn inner_prod<F: PrimeField>(a: &[F], b: &[F]) -> F {
    a.iter().zip(b.iter()).fold(F::ZERO, |acc, (x,y)| acc + *x*y)
}

pub fn BP_commit(g: &[G1Affine], h:&[G1Affine], a: &[Scalar], b: &[Scalar], U: G1Affine, c: Scalar) -> G1Affine 
{
    (best_multiexp(a, g) + best_multiexp(b, h) + U * c).into()
}