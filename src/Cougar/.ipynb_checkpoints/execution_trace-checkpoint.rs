// Execution Trace Constructor
use halo2curves::grumpkin::{
    G1, G1Affine, Fr as Scalar, Fq as Base,
};

use group::Group;
use ff::{Field, PrimeField};
use subtle::{Choice, ConditionallySelectable};

// Helper Function
pub fn bit2base(bit: Choice) -> Base {
    if bit.into() { Base::ONE} 
    else { Base::ZERO }
}

// P_{i+1} = x^{-1} * L_{i} + P + x * R_{i}
pub fn op_PTROW(
    // Group elements
    L_vec: &[G1Affine], R_vec: &[G1Affine], P_vec: &[G1Affine],
    // Scalars
    x: Scalar, x_inv: Scalar,
    // Execution Trace
    ext: &mut [Vec<Base>], perm: &mut[Vec<(usize, usize)>], offset: usize,
) -> (Vec<G1Affine>, usize)
{
    let n = L_vec.len();
    let mut ret: Vec<G1Affine> = Vec::with_capacity(n);
    let mut cursor = offset;
    
    // Helper
    let one = Base::ONE; let zero = Base::ZERO;
    
    for i in 0..n {
        // Compute x_inv * L_i
        let mut acc = G1::identity();
        let mut cnt = 0;
        
        for bit in x_inv
                .to_repr()
                .iter()
                .rev().flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
        {
            acc = acc.double();
            
            ext[0][cursor] = bit2base(bit); ext[1][cursor] = one;
            ext[2][cursor] = acc.x; ext[3][cursor] = acc.y; ext[4][cursor] = acc.z;
            ext[5][cursor] = L_vec[i].x; ext[6][cursor] = L_vec[i].y;
            cursor += 1;
            
            acc = G1::conditional_select(&acc, &(acc + L_vec[i]), bit);
            
            ext[0][cursor] = one; ext[1][cursor] = zero;
            ext[2][cursor] = acc.x; ext[3][cursor] = acc.y; ext[4][cursor] = acc.z;
            ext[5][cursor] = L_vec[i].x; ext[6][cursor] = L_vec[i].y;
            cursor += 1;
        }
        
        // Overlapping
        ext[0][cursor-1] = zero; ext[1][cursor-1] = zero;
        ext[2][cursor-1] = acc.x; ext[3][cursor-1] = acc.y; ext[4][cursor-1] = acc.z;
        let L_hat = acc;
        let L_hat_aff: G1Affine = L_hat.into();
        ext[5][cursor-1] = L_hat_aff.x; ext[6][cursor-1] = L_hat_aff.y;        
        
        // We record current position for L.
        let perm_L_start = cursor - 1;
        
        // Compute x * R_i
        let mut acc = G1::identity();
        for bit in x
                .to_repr()
                .iter()
                .rev().flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
        {
            acc = acc.double();
            
            ext[0][cursor] = bit2base(bit); ext[1][cursor] = one;
            ext[2][cursor] = acc.x; ext[3][cursor] = acc.y; ext[4][cursor] = acc.z;
            ext[5][cursor] = R_vec[i].x; ext[6][cursor] = R_vec[i].y;
            cursor += 1;
            
            acc = G1::conditional_select(&acc, &(acc + R_vec[i]), bit);
            
            ext[0][cursor] = one; ext[1][cursor] = zero;
            ext[2][cursor] = acc.x; ext[3][cursor] = acc.y; ext[4][cursor] = acc.z;
            ext[5][cursor] = R_vec[i].x; ext[6][cursor] = R_vec[i].y;
            cursor += 1;
        }        
        
        // Overlapping
        ext[0][cursor-1] = zero; ext[1][cursor-1] = zero;
        ext[2][cursor-1] = acc.x; ext[3][cursor-1] = acc.y; ext[4][cursor-1] = acc.z;
        let R_hat: G1Affine = acc.into();
        ext[5][cursor-1] = R_hat.x; ext[6][cursor-1] = R_hat.y;        
        
        // Record r's position
        let perm_R_start = cursor - 1;
        
        // Compute (x * L_i) + P_i
        let LP = L_hat + P_vec[i];
        ext[0][cursor] = one; ext[1][cursor] = one;
        ext[2][cursor] = L_hat.x; ext[3][cursor] = L_hat.y; ext[4][cursor] = L_hat.z;
        ext[5][cursor] = P_vec[i].x; ext[6][cursor] = P_vec[i].y;        
        
        // Record Permutations
        perm[2].push((perm_L_start, cursor));
        perm[3].push((perm_L_start, cursor));
        perm[4].push((perm_L_start, cursor));
        
        
        cursor += 1;
        
        // Compute (x * L_i) + P_i + (x^{-1} * R_i)
        let LPR = LP + R_hat;
        ext[0][cursor] = one; ext[1][cursor] = one;
        ext[2][cursor] = LP.x; ext[3][cursor] = LP.y; ext[4][cursor] = LP.z;
        ext[5][cursor] = R_hat.x; ext[6][cursor] = R_hat.y;        
        
        // Record Permutations
        perm[5].push((perm_R_start, cursor));
        perm[6].push((perm_R_start, cursor));
        
        cursor += 1; 
        
        ext[0][cursor] = zero; ext[1][cursor] = zero;
        ext[2][cursor] = LPR.x; ext[3][cursor] = LPR.y; ext[4][cursor] = LPR.z;
        let LPR_aff: G1Affine = LPR.into();
        ext[5][cursor] = LPR_aff.x; ext[6][cursor] = LPR_aff.y;        
        
        cursor += 1;
        
        // Add Dummy Operations 
        for j in 0..(512-3) {
            ext[0][cursor] = zero; ext[1][cursor] = zero;
            ext[2][cursor] = LPR.x; ext[3][cursor] = LPR.y; ext[4][cursor] = LPR.z;
            ext[5][cursor] = LPR_aff.x; ext[6][cursor] = LPR_aff.y;        
            
            cursor += 1;
        }
        
        
        
        ret.push(LPR_aff);
    }
    let diff = cursor - offset;
    (ret, diff )
}

// P_{i+1} = (P1 + xP2 || P3 + x^{-1}P4)
pub fn op_PTCOL(
    // Group elements
    P1_vec: &[G1Affine], P2_vec: &[G1Affine], P3_vec: &[G1Affine], P4_vec: &[G1Affine],
    // Scalars
    x: Scalar, x_inv: Scalar,
    // Execution Trace
    ext: &mut [Vec<Base>], offset: usize,
) -> (Vec<G1Affine>, usize)
{
    let n = P1_vec.len();
    let mut ret: Vec<G1Affine> = Vec::with_capacity(2 * n);
    
    // Helper
    let one = Base::ONE; let zero = Base::ZERO;
    let mut cursor = offset;
    
    // First Parts
    // We will use 511 rows & 1 dummies.
    for i in 0..n {
        // Compute x * P_2
        let mut acc = G1::identity();
        let mut cnt = 0;
        
        for bit in x
                .to_repr()
                .iter()
                .rev().flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
        {
            if (cnt == 0) {
                cnt += 1;
                continue;
            } 
            else {
                acc = acc.double();

                ext[0][cursor] = bit2base(bit); ext[1][cursor] = one;
                ext[2][cursor] = acc.x; ext[3][cursor] = acc.y; ext[4][cursor] = acc.z;
                ext[5][cursor] = P2_vec[i].x; ext[6][cursor] = P2_vec[i].y;
                cursor += 1;

                acc = G1::conditional_select(&acc, &(acc + P2_vec[i]), bit);

                ext[0][cursor] = one; ext[1][cursor] = zero;
                ext[2][cursor] = acc.x; ext[3][cursor] = acc.y; ext[4][cursor] = acc.z;
                ext[5][cursor] = P2_vec[i].x; ext[6][cursor] = P2_vec[i].y;
                cursor += 1;
            }
        }

        // Overlapping
        ext[0][cursor-1] = one; ext[1][cursor-1] = one;
        ext[2][cursor-1] = acc.x; ext[3][cursor-1] = acc.y; ext[4][cursor-1] = acc.z;
        ext[5][cursor-1] = P1_vec[i].x; ext[6][cursor-1] = P1_vec[i].y;      
        
        
        // Do operation
        let P_hat = acc + P1_vec[i];
        ext[0][cursor] = zero; ext[1][cursor] = zero;
        ext[2][cursor] = P_hat.x; ext[3][cursor] = P_hat.y; ext[4][cursor] = P_hat.z;
        let P_hat_aff: G1Affine = P_hat.into();
        ext[5][cursor] = P_hat_aff.x; ext[6][cursor] = P_hat_aff.y;              
        cursor += 1;
        
        // Dummy Data
        ext[0][cursor] = zero; ext[1][cursor] = zero;
        ext[2][cursor] = P_hat.x; ext[3][cursor] = P_hat.y; ext[4][cursor] = P_hat.z;
        ext[5][cursor] = P_hat_aff.x; ext[6][cursor] = P_hat_aff.y;              
        cursor += 1;        
        
        ret.push(P_hat_aff);
    }
    
    
    // Second Parts
    for i in 0..n {
        // Compute x * P_2
        let mut acc = G1::identity();
        let mut cnt = 0;
        
        for bit in x_inv
                .to_repr()
                .iter()
                .rev().flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
        {
            if (cnt == 0) {
                cnt += 1;
                continue;
            } 
            else {
                acc = acc.double();

                ext[0][cursor] = bit2base(bit); ext[1][cursor] = one;
                ext[2][cursor] = acc.x; ext[3][cursor] = acc.y; ext[4][cursor] = acc.z;
                ext[5][cursor] = P4_vec[i].x; ext[6][cursor] = P4_vec[i].y;
                cursor += 1;

                acc = G1::conditional_select(&acc, &(acc + P4_vec[i]), bit);

                ext[0][cursor] = one; ext[1][cursor] = zero;
                ext[2][cursor] = acc.x; ext[3][cursor] = acc.y; ext[4][cursor] = acc.z;
                ext[5][cursor] = P4_vec[i].x; ext[6][cursor] = P4_vec[i].y;
                cursor += 1;
            }
        }

        // Overlapping
        ext[0][cursor-1] = one; ext[1][cursor-1] = one;
        ext[2][cursor-1] = acc.x; ext[3][cursor-1] = acc.y; ext[4][cursor-1] = acc.z;
        ext[5][cursor-1] = P3_vec[i].x; ext[6][cursor-1] = P3_vec[i].y;      
        
        
        // Do operation
        let P_hat = acc + P3_vec[i];
        ext[0][cursor] = zero; ext[1][cursor] = zero;
        ext[2][cursor] = P_hat.x; ext[3][cursor] = P_hat.y; ext[4][cursor] = P_hat.z;
        let P_hat_aff: G1Affine = P_hat.into();
        ext[5][cursor] = P_hat_aff.x; ext[6][cursor] = P_hat_aff.y;              
        cursor += 1;
        
        // Dummy Data
        ext[0][cursor] = zero; ext[1][cursor] = zero;
        ext[2][cursor] = P_hat.x; ext[3][cursor] = P_hat.y; ext[4][cursor] = P_hat.z;
        ext[5][cursor] = P_hat_aff.x; ext[6][cursor] = P_hat_aff.y;              
        cursor += 1;        
        
        ret.push(P_hat_aff);
    }
    let diff = cursor - offset;
    
    (ret, diff)
}