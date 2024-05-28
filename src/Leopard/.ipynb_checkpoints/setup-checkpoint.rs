use blstrs::*;
use rand::thread_rng;
use group::Group;

pub struct LeopardParams {
    pub n: usize,
    pub m: usize,
    pub g: Vec<G1Affine>,
    pub h: Vec<G1Affine>,
    pub H: Vec<G2Affine>,
    pub U: Gt, 
}

use halo2curves::pairing;
use halo2curves::CurveAffine;


impl LeopardParams {
    pub fn new(N: usize) -> Self {
        let lgN = N.next_power_of_two().trailing_zeros();
        let lgm = (lgN / 2) + (lgN % 2);
        let lgn = lgN - lgm;

        let m = 1<<lgm;
        let n = 1<<lgn;
        let mut rng = thread_rng();
        let g: Vec<G1Affine> = (0..m).into_iter().map(|_| G1Projective::random(&mut rng).into()).collect();
        let h: Vec<G1Affine> = (0..m).into_iter().map(|_| G1Projective::random(&mut rng).into()).collect();
        let H: Vec<G2Affine> = (0..n).into_iter().map(|_| G2Projective::random(&mut rng).into()).collect();
        
        let mut rng = thread_rng();
        let U = Gt::random(&mut rng);
        Self {
            n, m, g, h, H, U,
        }
    }
}