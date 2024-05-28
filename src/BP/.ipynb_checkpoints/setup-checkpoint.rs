use halo2curves::secp256k1::{Secp256k1Affine as G1Affine, Fq as Scalar};
use rand::thread_rng;
use group::Group;

pub struct BPParams {
    pub n: usize,
    pub G: Vec<G1Affine>,
    pub H: Vec<G1Affine>,
    pub U: G1Affine, 
}


impl BPParams {
    pub fn new(n: usize) -> Self {
        let G = make_vec(n);
        let H = make_vec(n);
    
        let mut rng = thread_rng();
        let U = G1Affine::random(&mut rng);
        Self {
            n, G, H, U
        }
    }
}

pub fn make_vec(n: usize) -> Vec<G1Affine> {
    let mut rng = thread_rng();
    let ret: Vec<G1Affine> = (0..n).into_iter().map(|_| G1Affine::random(&mut rng)).collect();
    ret
}