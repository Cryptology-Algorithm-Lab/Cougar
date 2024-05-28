use ff::Field;
use blstrs::*;
use pairing::{MillerLoopResult, MultiMillerLoop};


pub fn inner_prod<F: Field>(a: &[F], b: &[F]) -> F {
    a.iter().zip(b.iter()).fold(F::ZERO, |acc, (x,y)| acc + *x*y)
}

// Calculate L, R on Column Reduction
pub fn calc_LR1(
    g: &[G1Affine], 
    h: &[G1Affine], 
    H: &[G2Affine],
    a: &[Scalar], 
    b: &[Scalar]) -> Gt {
        
        // Prepare for pairing
        // Left Term (G1)
        let mut left: Vec<G1Affine> = Vec::with_capacity(H.len());
        let n = H.len();
        let m = g.len();
        
        
        // Preprocessing 
        for j in 0..n {
            
            // Setup for base & exponents
            let base: Vec<G1Projective> = g.into_iter().
                        map(|x| x.into()).chain(h.into_iter().map(|y| y.into())).collect();
            
            let expon: Vec<Scalar> = (0..m).into_iter().map(|i| a[i*n+j])
                                            .chain((0..m)
                                            .into_iter().map(|i| b[i*n+j])).collect();
            
            let pip_res: G1Affine = G1Projective::multi_exp(&base, &expon).into();
            left.push(pip_res);
        }
        
        // Right Term (G2)
        let right: Vec<G2Prepared> = H.into_iter().map(|x| G2Prepared::from(*x)).collect();
        
        let pre_pairing: Vec<(&G1Affine, &G2Prepared)> = left.iter().zip(right.iter()).collect();
        Bls12::multi_miller_loop(&pre_pairing).final_exponentiation()
}

// Calculate L, R on Row Reduction
pub fn calc_LR2(    
    g: &[G1Affine], 
    h: &[G1Affine], 
    gH: &[G2Affine],
    hH: &[G2Affine],
    a: &[Scalar], 
    b: &[Scalar]) -> Gt {
        
        // Left Term (G1)
        let left: Vec<G1Affine> = a.into_iter().map(|x| (g[0]*x).into())
                                        .chain(b.into_iter().map(|x| (h[0]*x).into())).collect();
        // Right Term (G2)
        let right: Vec<G2Prepared> = gH.into_iter().map(|x| G2Prepared::from(*x))
                                   .chain(hH.into_iter().map(|x| G2Prepared::from(*x))).collect();
        
        let pre_pairing: Vec<(&G1Affine, &G2Prepared)> = left.iter().zip(right.iter()).collect();
        Bls12::multi_miller_loop(&pre_pairing).final_exponentiation()
}

// Batch Inversion with Montgomery's Trick
pub fn batch_invert(lg_m:usize, lg_n:usize, scalars: &mut [Scalar]) -> (Scalar, Scalar) {
    // In-place batch invert
    let lg_mn = scalars.len();
    let mut prods = vec![Scalar::ONE; lg_mn];
    let mut acc: Scalar = Scalar::ONE;
    
    
    // Forward Pass
    for (input, prod) in scalars.iter().zip(prods.iter_mut()) {
        *prod = acc;
        acc = acc * *input;
    }
    
    acc = acc.invert().unwrap();
    
    let inv_left = prods[lg_m];
    let inv_right = acc * inv_left.invert().unwrap();
    
    // backward pass
    for (input, prod) in scalars.iter_mut().rev().zip(prods.into_iter().rev()) {
        let tmp = acc * *input;
        *input = acc * prod;
        acc = tmp;
    }
    
    (inv_left, inv_right)
    
}