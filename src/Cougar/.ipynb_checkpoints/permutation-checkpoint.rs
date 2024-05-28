use halo2_proofs::poly::{EvaluationDomain, LagrangeCoeff as LC, ExtendedLagrangeCoeff as ELC};
use halo2_proofs::poly::{Evaluator, new_evaluator, AstLeaf, AstMul, Ast, Basis, Rotation, Polynomial};
use halo2curves::grumpkin::{Fq as Base};
use ff::{Field, PrimeField, BatchInvert};

// Construct IDPoly and Permuted Poly
pub fn get_id_and_permpolys(
    D: &EvaluationDomain<Base>,
    // Column Length
    N: usize,
    perms: &[Vec<(usize, usize)>],
) -> (Vec<Vec<Base>>, Vec<Vec<Base>>)
{
    // Discard Selectors
    let M = perms.len() - 2;
    let mut idpolys: Vec<Vec<Base>> = vec![vec![Base::ONE; N]; M];
    
    
    
    let omega = D.get_omega();
    let mut deltaomega = Base::ONE;
    
    assert_eq!(omega.pow_vartime([N as u64, 0, 0, 0]), Base::ONE);
    
    for i in 0..M {
        for j in 0..N {
            idpolys[i][j] *= deltaomega;
            deltaomega *= omega;
        }
        deltaomega *= Base::DELTA;
    };
    
    let mut permpolys = idpolys.clone();
    
    // Swap by permutation
    for i in 0..M {
        for (left, right) in perms[i+2].clone().into_iter() {
            let tmp = permpolys[i][left].clone();
            permpolys[i][left] = permpolys[i][right].clone();
            permpolys[i][right] = tmp;
        }
    }
    
    (idpolys, permpolys)
}

// Permutation Polynomial
pub fn get_z_poly(
    D: &EvaluationDomain<Base>,
    cols: &[Vec<Base>],
    idpolys: &[Vec<Base>],
    permpolys: &[Vec<Base>],
    u1: Base, u2: Base
) -> Vec<Base> {
    let M = cols.len();
    let N = cols[0].len();
    
    let mut ret: Vec<Base> = vec![Base::ONE; N];
    
    // Denominator First!
    for i in 2..M {
        let mut acc = Base::ONE;
        for j in 0..N {
            ret[j] *= acc;
            acc *= (cols[i][j] + u1 * permpolys[i-2][j] + u2);
        }
    }
    
    ret.batch_invert();
    
    // Final Product
    for i in 2..M {
        let mut acc = Base::ONE;
        for j in 0..N {
            ret[j] *= acc;
            acc *= (cols[i][j] + u1 * idpolys[i-2][j] + u2);
        }
    }
    
    ret
}