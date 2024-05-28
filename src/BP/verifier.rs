use crate::BP::transcript::{
    Transcript,
    TranscriptRead,
    Blake2bRead,
    Challenge255,
};

use halo2curves::secp256k1::{
    Secp256k1Affine as G1Affine, Fq as Scalar                            
};

use halo2curves::msm::best_multiexp;

use crate::BP::setup::{
    BPParams,
};

use ff::BatchInvert;

pub fn BP_verifier(
    BP: &BPParams,
    transcript: &mut Blake2bRead<&[u8], G1Affine, Challenge255<G1Affine>>,
)
{
    let n = BP.n;
    let lgn = n.next_power_of_two().trailing_zeros() as usize;

    let gamma = *transcript.squeeze_challenge_scalar::<()>();
    let U = BP.U * gamma;
    let P = transcript.read_point().unwrap();

    let mut L_vec : Vec<G1Affine> = Vec::with_capacity(lgn);
    let mut R_vec : Vec<G1Affine> = Vec::with_capacity(lgn);
    let mut x_vec : Vec<Scalar> = Vec::with_capacity(lgn);

    for i in 0..lgn {
        L_vec.push(transcript.read_point().unwrap());
        R_vec.push(transcript.read_point().unwrap());
        x_vec.push(*transcript.squeeze_challenge_scalar::<()>());
    }

    let a = transcript.read_scalar().unwrap();
    let b = transcript.read_scalar().unwrap();
        
    let mut x_inv_vec = x_vec.clone();
    x_inv_vec.batch_invert();
    
    let mut s_vec = Vec::with_capacity(n);
    let init = x_inv_vec.iter().fold(Scalar::one(), |acc, x| acc * x);
    s_vec.push(init);

    for (x, y) in x_vec.iter_mut().zip(x_inv_vec.iter_mut()) {
        *x = x.square();
        *y = y.square();
    }

    for i in 1..n {
        let lgi: usize = (32 - 1 - (i as u32).leading_zeros()) as usize;
        let k = 1<<lgi;
        let u = x_vec[lgn - 1 - lgi];
        s_vec.push(s_vec[i-k] * u);
    }

    let left: G1Affine = {
        let g_exp: Vec<Scalar> = s_vec.iter().map(|x| a * x).collect();
        let h_exp: Vec<Scalar> = s_vec.iter().rev().map(|x| b * x).collect();

        (best_multiexp(&g_exp, &BP.G) + best_multiexp(&h_exp, &BP.H) + (U * (a*b))).into()
    };

    let right: G1Affine = {
        (P + best_multiexp(&x_vec, &L_vec) + best_multiexp(&x_inv_vec, &R_vec)).into()
    };

    assert_eq!(left, right);
}