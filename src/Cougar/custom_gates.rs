// Custom Gate Constructors
// All constraints are combined by a random challenge gamma

use halo2_proofs::poly::{EvaluationDomain, LagrangeCoeff as LC, ExtendedLagrangeCoeff as ELC, Coeff};
use halo2_proofs::poly::{Evaluator, new_evaluator, AstLeaf, AstMul, Ast, Basis, Rotation, Polynomial};
use halo2curves::grumpkin::{Fq as Base};

use ff::Field;


// Polynomial Constructor for Prover;
// sels: Selector
// cols: Columns
// onepoly: Helper Polynomial; A NNT-ed 1 polynomial.
pub fn constructor_T(D: &EvaluationDomain<Base>, cols: &[Polynomial<Base, ELC>], 
    onepoly: &Polynomial<Base, ELC>, eta: Base) -> Polynomial<Base, ELC>
{
    let mut eval: Evaluator<_, Base, ELC> = new_evaluator(|| {});
    let s1 = eval.register_poly(cols[0].clone());
    let s2 = eval.register_poly(cols[1].clone());
    let px = eval.register_poly(cols[2].clone());
    let py = eval.register_poly(cols[3].clone());
    let pz = eval.register_poly(cols[4].clone());
    let gx = eval.register_poly(cols[5].clone());
    let gy = eval.register_poly(cols[6].clone());
    
    let px_nxt = px.with_rotation(Rotation::next());
    let py_nxt = py.with_rotation(Rotation::next());
    let pz_nxt = pz.with_rotation(Rotation::next());
    
    // Encode into Ast.
    let s1 = Ast::<_, Base, _>::from(s1);
    let s2 = Ast::<_, Base, _>::from(s2);
    let px = Ast::<_, Base, _>::from(px);
    let py = Ast::<_, Base, _>::from(py);
    let pz = Ast::<_, Base, _>::from(pz);
    let gx = Ast::<_, Base, _>::from(gx);
    let gy = Ast::<_, Base, _>::from(gy);
    let px_nxt = Ast::<_, Base, _>::from(px_nxt);
    let py_nxt = Ast::<_, Base, _>::from(py_nxt);
    let pz_nxt = Ast::<_, Base, _>::from(pz_nxt);
    let onepoly = Ast::<_, Base, _>::from(eval.register_poly(onepoly.clone()));
    
    // Grumpkin Curve Constant
    let b3 = Base::from(3) * (-Base::from(17));
    
    // Selector 11: Mixed Addition
    let t0 = &px * &gx; let t1 = &py * &gy; let t3 = &gx + &gy;
    let t4 = &px + &py; let t3 = &t3 * &t4; let t4 = &t0 + &t1;
    let t3 = &t3 - &t4; let t4 = &gy * &pz; let t4 = &t4 + &py;
    let y3 = &gx * &pz; let y3 = &y3 + &px; let x3 = &t0 + &t0;
    let t0 = &x3 + &t0; let t2 = &pz * b3;  let z3 = &t1 + &t2;
    let t1 = &t1 - &t2; let y3 = &y3 * b3;  let x3 = &t4 * &y3;
    let t2 = &t3 * &t1; let x3 = &t2 - &x3; let y3 = &y3 * &t0;
    let t1 = &t1 * &z3; let y3 = &t1 + &y3; let t0 = &t0 * &t3;
    let z3 = &z3 * &t4; let z3 = &z3 + &t0;
    
    let w1 = &px_nxt - &x3;
    let w2 = &py_nxt - &y3;
    let w3 = &pz_nxt - &z3;
    let G1 = w1 + w2 * eta + w3 * eta.square();
    
    // Selector 01: Identity Constraints
    let w1 = &px_nxt - &px;
    let w2 = &py_nxt - &py;
    let w3 = &pz_nxt - &pz;
    let G2 = w1 + w2 * eta + w3 * eta.square();
    
    // Selector 10: Doubling
    let t0 = &py * &py; let z3 = &t0 + &t0; let z3 = &z3 + &z3;
    let z3 = &z3 + &z3; let t1 = &py * &pz; let t2 = &pz * &pz;
    let t2 = &t2 * b3;  let x3 = &t2 * &z3; let y3 = &t0 + &t2;
    let z3 = &t1 * &z3; let t1 = &t2 + &t2; let t2 = &t1 + &t2;
    let t0 = &t0 - &t2; let y3 = &t0 * &y3; let y3 = &x3 + &y3;
    let t1 = &px * &py; let x3 = &t0 * &t1; let x3 = &x3 + &x3;
    
    let w1 = &px_nxt - &x3;
    let w2 = &py_nxt - &y3;
    let w3 = &pz_nxt - &z3;
    let G3 = &w1 + &(&w2 * eta) + (&w3 * eta.square());
    
    // Selector 00: Projection to Affine
    let w1 = &(&gx * &pz) - &px;
    let w2 = &(&gy * &pz) - &py;
    let G4 = &w1 + &(&w2 * eta);
    
    // Selector Combination
    let s1m = &onepoly - &s1; let s2m = &onepoly - &s2;
    let s11 = &s1 * &s2; let s10 = &s1 * &s2m; let s01 = &s1m * &s2; let s00 = &s1m * &s2m;
    let G1 = &s11 * &G1; let G2 = &s01 * &G2; let G3 = &s10 * &G3; let G4 = &s00 * &G4;
    
    // Total: Degree 6 Circuit.
    eval.evaluate(&(G1 + G2 + G3 + G4), &D)
}


pub fn constructor_Z(
    D: &EvaluationDomain<Base>, 
    z_poly: &Polynomial<Base, ELC>,
    cols: &[Polynomial<Base, ELC>],
    ids: &[Polynomial<Base, ELC>],
    perms: &[Polynomial<Base, ELC>],
    onepoly: &Polynomial<Base, ELC>,
    H1: &Polynomial<Base, ELC>,
    u1: Base, u2:Base,
) -> (Polynomial<Base, ELC>, Polynomial<Base, ELC>)
{
    let mut eval: Evaluator<_, Base, ELC> = new_evaluator(|| {});
    let z = eval.register_poly(z_poly.clone());
    let z_nxt = z.clone().with_rotation(Rotation::next());
    
    let z = Ast::<_, Base, _>::from(z.clone());
    let z_nxt = Ast::<_, Base, _>::from(z_nxt.clone());
    
    let v0 = eval.register_poly(cols[0].clone());
    let v1 = eval.register_poly(cols[1].clone());
    let v2 = eval.register_poly(cols[2].clone());
    let v3 = eval.register_poly(cols[3].clone());
    let v4 = eval.register_poly(cols[4].clone());
    
    let id0 = eval.register_poly(ids[0].clone());
    let id1 = eval.register_poly(ids[1].clone());
    let id2 = eval.register_poly(ids[2].clone());
    let id3 = eval.register_poly(ids[3].clone());
    let id4 = eval.register_poly(ids[4].clone());

    let r0 = eval.register_poly(perms[0].clone());
    let r1 = eval.register_poly(perms[1].clone());
    let r2 = eval.register_poly(perms[2].clone());
    let r3 = eval.register_poly(perms[3].clone());
    let r4 = eval.register_poly(perms[4].clone());
    
    let v0 = Ast::<_, Base, _>::from(v0);
    let v1 = Ast::<_, Base, _>::from(v1);
    let v2 = Ast::<_, Base, _>::from(v2);
    let v3 = Ast::<_, Base, _>::from(v3);
    let v4 = Ast::<_, Base, _>::from(v4);
    
    let id0 = Ast::<_, Base, _>::from(id0);
    let id1 = Ast::<_, Base, _>::from(id1);
    let id2 = Ast::<_, Base, _>::from(id2);
    let id3 = Ast::<_, Base, _>::from(id3);
    let id4 = Ast::<_, Base, _>::from(id4);
    
    let r0 = Ast::<_, Base, _>::from(r0);
    let r1 = Ast::<_, Base, _>::from(r1);
    let r2 = Ast::<_, Base, _>::from(r2);
    let r3 = Ast::<_, Base, _>::from(r3);
    let r4 = Ast::<_, Base, _>::from(r4);
   
    let one = Ast::<_, Base, _>::from(eval.register_poly(onepoly.clone()));
    let u2poly = Ast::<_, Base, _>::from(eval.register_poly(onepoly.clone()));
    let u2poly = &u2poly * u2;
    let H1 = Ast::<_, Base, _>::from(eval.register_poly(H1.clone()));
    
    // T_mid
    let acc_left = {
        let add = &v0 + &u2poly; let mul = &id0 * u1; let acc = &add + &mul;
        let add = &v1 + &u2poly; let mul = &id1 * u1; 
        let tmp = &add + &mul; let acc = &acc * &tmp;
        let add = &v2 + &u2poly; let mul = &id2 * u1; 
        let tmp = &add + &mul; let acc = &acc * &tmp;
        let add = &v3 + &u2poly; let mul = &id3 * u1; 
        let tmp = &add + &mul; let acc = &acc * &tmp;
        let add = &v4 + &u2poly; let mul = &id4 * u1; 
        let tmp = &add + &mul; let acc = &acc * &tmp;
        let acc = &acc * &z;
        acc
    };
    
    let acc_right = {
        let add = &v0 + &u2poly; let mul = &r0 * u1; let acc = &add + &mul;
        let add = &v1 + &u2poly; let mul = &r1 * u1; 
        let tmp = &add + &mul; let acc = &acc * &tmp;
        let add = &v2 + &u2poly; let mul = &r2 * u1; 
        let tmp = &add + &mul; let acc = &acc * &tmp;
        let add = &v3 + &u2poly; let mul = &r3 * u1; 
        let tmp = &add + &mul; let acc = &acc * &tmp;
        let add = &v4 + &u2poly; let mul = &r4 * u1; 
        let tmp = &add + &mul; let acc = &acc * &tmp;
        let acc = &acc * &z_nxt;
        acc
    };
    
    let T_mid = &acc_left - &acc_right;
    let T_mid = eval.evaluate(&T_mid, D);
    
    // T_right
    let z_1 = &z - &one;
    let T_right = &z_1 * &H1;
    let T_right = eval.evaluate(&T_right, D);
    
    (T_mid, T_right)
}


// Polynomial Evaluation Results for Verifier
pub fn poly_eval_T(cols: &[Base], cols_nxt: &[Base], eta: Base) -> Base
{
    // Encoding
    let s1 = cols[0]; let s2 = cols[1];
    let s11 = s1 * s2; let s10 = s1 * (Base::ONE - s2); 
    let s01 = (Base::ONE - s1) * s2; let s00 = (Base::ONE - s1) * (Base::ONE - s2);
    
    let px = cols[2]; let py = cols[3]; let pz = cols[4];
    let gx = cols[5]; let gy = cols[6];
    let px_nxt = cols_nxt[0]; let py_nxt = cols_nxt[1]; let pz_nxt = cols_nxt[2];
    
    // Grumpkin Curve Constant
    let b3 = Base::from(3) * (-Base::from(17));
    
    // Mixed Addition
    let t0 = px * gx; let t1 = py * gy; let t3 = gx + gy;
    let t4 = px + py; let t3 = t3 * t4; let t4 = t0 + t1;
    let t3 = t3 - t4; let t4 = gy * pz; let t4 = t4 + py;
    let y3 = gx * pz; let y3 = y3 + px; let x3 = t0 + t0;
    let t0 = x3 + t0; let t2 = pz * b3; let z3 = t1 + t2;
    let t1 = t1 - t2; let y3 = y3 * b3; let x3 = t4 * y3;
    let t2 = t3 * t1; let x3 = t2 - x3; let y3 = y3 * t0;
    let t1 = t1 * z3; let y3 = t1 + y3; let t0 = t0 * t3;
    let z3 = z3 * t4; let z3 = z3 + t0;
    
    let w1 = px_nxt - x3;
    let w2 = py_nxt - y3;
    let w3 = pz_nxt - z3;
    let G1 = w1 + w2 * eta + w3 * eta.square();
    
    // Identity
    let w1 = px_nxt - px;
    let w2 = py_nxt - py;
    let w3 = pz_nxt - pz;
    let G2 = w1 + w2 * eta + w3 * eta.square();
    
    // Doubling
    let t0 = py * py; let z3 = t0 + t0; let z3 = z3 + z3;
    let z3 = z3 + z3; let t1 = py * pz; let t2 = pz * pz;
    let t2 = t2 * b3; let x3 = t2 * z3; let y3 = t0 + t2;
    let z3 = t1 * z3; let t1 = t2 + t2; let t2 = t1 + t2;
    let t0 = t0 - t2; let y3 = t0 * y3; let y3 = x3 + y3;
    let t1 = px * py; let x3 = t0 * t1; let x3 = x3 + x3;
    
    let w1 = px_nxt - x3;
    let w2 = py_nxt - y3;
    let w3 = pz_nxt - z3;
    let G3 = w1 + w2 * eta + w3 * eta.square();
    
    // Selector 00: Projection to Affine
    let w1 = (gx * pz) - px;
    let w2 = (gy * pz) - py;
    let G4 = w1 + w2 * eta;
    
    // Selector Combination
    s11 * G1 + s01 * G2 + s10 * G3 + s00 * G4
        
}

pub fn poly_eval_Z(
    z: Base,
    z_nxt: Base,
    cols: &[Base],
    ids: &[Base],
    perms: &[Base],
    u1: Base, u2:Base, u3:Base, H1_eval: Base,
) -> Base
{
    let acc_left = cols.iter().zip(ids.iter())
                               .fold(Base::ONE, |acc, (x,y)| acc * (x + u1 * y + u2));
    let acc_right = cols.iter().zip(perms.iter())
                               .fold(Base::ONE, |acc, (x,y)| acc * (x + u1 * y + u2));    
    
    let T_mid_eval = z * acc_left - z_nxt * acc_right;
    let T_right_eval = (z - Base::ONE) * H1_eval;
    
    T_mid_eval * u3 + T_right_eval * u3.square()
}
