use crate::BP::transcript::{
    Transcript,
    TranscriptWrite,
    Blake2bWrite,
    Challenge255,
};

use crate::BP::setup::{
    BPParams,
};

use halo2curves::secp256k1::{
    Secp256k1Affine as G1Affine, Fq as Scalar
};

use crate::BP::utils::{
    inner_prod, BP_commit,
};

pub fn BP_prover (
    BP: &BPParams,
    transcript: &mut Blake2bWrite<Vec<u8>, G1Affine, Challenge255<G1Affine>>,
    a_vec: Vec<Scalar>,
    b_vec: Vec<Scalar>,
)
{
    let mut n = BP.n;
    let mut a = &mut a_vec.clone()[..];
    let mut b = &mut b_vec.clone()[..];
    let mut G = &mut BP.G.clone()[..];
    let mut H = &mut BP.H.clone()[..];

    let gamma: Scalar = *transcript.squeeze_challenge_scalar::<()>();
    let U: G1Affine = (BP.U * gamma).into();
    let c = inner_prod(a, b);
    let P = BP_commit(G, H, a, b, U, c);
    transcript.write_point(P);
    let lgn = n.next_power_of_two().trailing_zeros();
    
    for i in 0..lgn {
        n >>= 1;
        let (a_L, a_R) = a.split_at_mut(n);
        let (b_L, b_R) = b.split_at_mut(n);
        let (G_L, G_R) = G.split_at_mut(n);
        let (H_L, H_R) = H.split_at_mut(n);

        let c_L = inner_prod(a_L, b_R);
        let c_R = inner_prod(a_R, b_L);

        let L = BP_commit(G_R, H_L, a_L, b_R, U, c_L);
        let R = BP_commit(G_L, H_R, a_R, b_L, U, c_R);

        transcript.write_point(L);
        transcript.write_point(R);

        let x: Scalar = *transcript.squeeze_challenge_scalar::<()>();
        let x_inv = x.invert().unwrap();

        for i in 0..n {
            a_L[i] = x * a_L[i] + x_inv * a_R[i];
            b_L[i] = x_inv * b_L[i] + x * b_R[i];
            G_L[i] = (G_L[i] * x_inv + G_R[i] * x).into();
            H_L[i] = (H_L[i] * x + H_R[i] * x_inv).into();
        }

        a = a_L;
        b = b_L;
        G = G_L;
        H = H_L;
    }

    transcript.write_scalar(a[0]);
    transcript.write_scalar(b[0]);
}