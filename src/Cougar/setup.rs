
use halo2curves::grumpkin::{
    G1Affine as GG1Affine, Fr as GScalar, G1 as GG1
};

use halo2curves::bn256::{
    G1Affine as BG1Affine, G2Affine as BG2Affine, G1 as BG1, G2 as BG2,
    Fr as BScalar, Gt, G2Prepared as BG2Prepared
};

use halo2curves::bn256::{pairing, multi_miller_loop};
use halo2curves::pairing::MillerLoopResult;
use halo2curves::{CurveAffine, CurveExt};
use halo2curves::group::{Group, Curve};


use rand::thread_rng;

use halo2_proofs::poly::{
    EvaluationDomain,
    ExtendedLagrangeCoeff,
    Polynomial,
};

// Cougar Parameters
pub struct CougarParams {
    // Parameters
    pub m: usize, pub n: usize, pub d: usize, pub D: usize,
    // Evaluation domains
    pub eval_xi: EvaluationDomain<BScalar>,
    pub eval_zeta: EvaluationDomain<BScalar>,
    // G, H for 1st layer 
    pub G: Vec<GG1Affine>,
    pub H: Vec<GG1Affine>,
    // Commit Key for PC
    pub ck_pc: LeopardParams,
}

// Commit Keys for LeopardPC
pub struct LeopardParams {
    pub G: Vec<BG1Affine>,
    pub H: Vec<BG2Affine>,
    pub U: Gt,
}

// Parameter Selector
pub fn param_selector(N: usize) -> (usize, usize, usize, usize) {
    let lgN = N.next_power_of_two().trailing_zeros() as usize;
    let c1c2 = 4;
    let m3 = (c1c2 * c1c2) * 64 * (3 * N * lgN + 2 * N);
    let lgm3 = m3.next_power_of_two().trailing_zeros() as usize;
    let lgm = lgm3/3;
    let mut m = 1<<lgm;
    let mut n = N/m as usize;
    if (m >= N) {
        m = N;
        n = 1;
    } else {
        n = N / m;
    }
    let lgm = m.next_power_of_two().trailing_zeros() as usize;
    let d = 512;
    let D = 6 * n * (lgm) * d + 2 * n * d  + d * (2 * n - 2);
    let D = D.next_power_of_two() as usize;
    assert_eq!(N, m * n);
    (m as usize, n as usize, d as usize, D as usize)
}


// Setup Algorithm
impl CougarParams {
    pub fn new(N: usize) -> Self {
        // Choose m, n, D through Automatic Parameter Selector
        // m: Length of commit key on 1st Layer
        // n: Length of commit key on 2nd Layer
        // d: Length of sparcified L, R.
        // D: Length of commit key for Leopard
        let (m, n, d, D) = param_selector(N);
        let logD = D.next_power_of_two().trailing_zeros();
        let logd = d.next_power_of_two().trailing_zeros();
        let eval_xi: EvaluationDomain<BScalar> = EvaluationDomain::new(6, logD);    
        let eval_zeta: EvaluationDomain<BScalar> = EvaluationDomain::new(d as u32, logD-logd);
        let G = make_param::<GG1Affine>(m);
        let H = make_param::<GG1Affine>(m);
        let ck_pc = LeopardParams::new(D);
        
        CougarParams{
            m,n,d,D,
            eval_xi, eval_zeta,
            G, H,
            ck_pc
        }
    }
}

// Setup Algorithm
impl LeopardParams {
    pub fn new(N: usize) -> Self {
        // Automatic Parameter Selector
        let logN = N.next_power_of_two().trailing_zeros();
        let (m, n) = (N >> (logN / 2), N >> ((logN / 2) + logN%2));
        
        let G = make_param::<BG1Affine>(m);
        let H = make_param::<BG2Affine>(n);
        
        let mut rng = thread_rng();
        let U = Gt::random(&mut rng);
        LeopardParams {
            G, H, U
        }
    }
}

pub fn make_param<C: CurveAffine>(n:usize) -> Vec<C>
{
    
    let mut rng = thread_rng();
    let mut ret: Vec<C> = Vec::with_capacity(n);
    
    for i in 0..n {
        ret.push(C::CurveExt::random(&mut rng).into());
    }
    ret
}