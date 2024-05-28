use crate::Cougar::prover::{
    cougar_prover
};


use crate::Cougar::verifier::{
    cougar_verifier
};

use crate::Cougar::setup::{
    CougarParams,
};

use ff::Field;
use group::Group;
use rand::thread_rng;


use crate::Cougar::utils::{
    inner_prod, Leopard_commit
};

use crate::Cougar::transcript::{
    Transcript,
    TranscriptWrite,
    TranscriptRead,
    Blake2bWrite,
    Blake2bRead,
    Challenge255,
};



pub fn test_single(lgN: u32) -> (f32, f32, usize) {
    use halo2curves::bn256::{
        Fr as BScalar, Fq as GScalar,
        G1Affine,
    };
    use std::time::Instant;

    // Setup
    let N = 1<<lgN;
    let ck = CougarParams::new(N);
    
    println!("N: 2^{}", lgN);
    println!("Setup is DONE!");
    println!("(m, n, D, d) = ({}, {}, {}, {})", ck.m, ck.n, ck.D, ck.d);
    
    let mut rng = thread_rng();
    let a_vec: Vec<GScalar> = (0..N).into_iter().map(|_| GScalar::random(&mut rng)).collect();
    let b_vec: Vec<GScalar> = (0..N).into_iter().map(|_| GScalar::random(&mut rng)).collect();
    
    let c = inner_prod(&a_vec, &b_vec);
    
    // Prover
    let mut transcript = Blake2bWrite::<Vec<u8>, G1Affine, Challenge255<G1Affine>>::init(vec![]);
        
    let time_proof = Instant::now();
    let proof = {
        cougar_prover(&ck, &mut transcript, a_vec, b_vec);
        let proof = transcript.finalize();
        proof
    };
    let proof_elapsed = time_proof.elapsed();

    
    println!("Proof Generation DONE!");
    println!("Time for Proof: {} sec", proof_elapsed.as_micros() as f32/ 1_000_000 as f32);
    println!("Proof Size: {} bytes", proof.len());
    
    // Verifier
    let mut transcript = Blake2bRead::<&[u8], G1Affine, Challenge255<G1Affine>>::init(&proof[..]);
    let time_verify = Instant::now();
    cougar_verifier(&ck, &mut transcript, c);
    let verify_elapsed = time_verify.elapsed();
    println!("Verification DONE!");
    println!("Time for Verify: {} sec", verify_elapsed.as_micros() as f32 / 1_000_000 as f32);

    println!("Yay!");
    (
        proof_elapsed.as_micros() as f32/ 1_000_000 as f32,
        verify_elapsed.as_micros() as f32 / 1_000_000 as f32,
        proof.len(),        
    )
}