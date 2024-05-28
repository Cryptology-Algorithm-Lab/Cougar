use crate::BP::prover::{
    BP_prover
};


use crate::BP::verifier::{
    BP_verifier
};

use crate::BP::setup::{
    BPParams,
};

use ff::Field;
use group::Group;
use rand::thread_rng;


use crate::BP::utils::{
    inner_prod, BP_commit,
};

use crate::BP::transcript::{
    Transcript,
    TranscriptWrite,
    TranscriptRead,
    Blake2bWrite,
    Blake2bRead,
    Challenge255,
};


pub fn test_single(lgN: u32) -> (f32, f32, usize) {
    use halo2curves::secp256k1::{
        Secp256k1Affine as G1Affine,
        Fq as Scalar,
    };
    use std::time::Instant;
    
    // Setup
    let N = 1<<lgN;
    let ck = BPParams::new(N);
    
    println!("N: 2^{}", lgN);
    println!("Setup is DONE!");
    
    let mut rng = thread_rng();
    let a_vec: Vec<Scalar> = (0..N).into_iter().map(|_| Scalar::random(&mut rng)).collect();
    let b_vec: Vec<Scalar> = (0..N).into_iter().map(|_| Scalar::random(&mut rng)).collect();
    
    // Prover
    let mut transcript = Blake2bWrite::<Vec<u8>, G1Affine, Challenge255<G1Affine>>::init(vec![]);
        
    let time_proof = Instant::now();
    let proof = {
        BP_prover(&ck, &mut transcript, a_vec, b_vec);
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
    BP_verifier(&ck, &mut transcript);
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