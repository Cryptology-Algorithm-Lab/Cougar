use crate::Leopard::IPA::InnerProductProof;
use merlin::Transcript;
use crate::Leopard::transcript::*;

use crate::Leopard::setup::{LeopardParams};

use rand::thread_rng;

use std::time::Instant;

use std::mem::size_of;

use ff::Field;

use blstrs::*;

// A helper algorithm for calculating proof size.
pub fn get_proof_size(IP: &InnerProductProof) -> usize {
    let size_gt = size_of::<Gt>();
    let size_Fr = size_of::<Scalar>();
        
    let n_gts = IP.L_vec.len() * 2;
    let n_frs = 2;

    size_gt * n_gts + size_Fr * n_frs
}

pub fn test_single(lgN: u32) -> (f32, f32, usize) {
    // Setup
    let N = 1<<lgN;
    let ck = LeopardParams::new(N);

    println!("N: 2^{}", lgN);
    println!("Setup is DONE!");
    
    let mut rng = thread_rng();
    let a_vec: Vec<Scalar> = (0..N).into_iter().map(|_| Scalar::random(&mut rng)).collect();
    let b_vec: Vec<Scalar> = (0..N).into_iter().map(|_| Scalar::random(&mut rng)).collect();

    // Prover
    let mut trans_prover = Transcript::new(b"Asiacrypt24_Cougar");
    
    let time_proof = Instant::now();
    let IP: InnerProductProof = InnerProductProof::proof(
                                &mut trans_prover, &ck, a_vec ,b_vec);
    let proof_elapsed = time_proof.elapsed();

    
    println!("Proof Generation DONE!");
    println!("Time for Proof: {} sec", proof_elapsed.as_micros() as f32/ 1_000_000 as f32);
    let proof_size = get_proof_size(&IP);
    println!("Proof Size: {} bytes", proof_size);

    // Verifier
    let mut trans_verifier = Transcript::new(b"Asiacrypt24_Cougar");
    let time_verify = Instant::now();
    IP.verify(&mut trans_verifier, &ck);
    let verify_elapsed = time_verify.elapsed();
    println!("Time for Verify: {} sec", verify_elapsed.as_micros() as f32 / 1_000_000 as f32);
    
    println!("\n");
    (
        proof_elapsed.as_micros() as f32/ 1_000_000 as f32,
        verify_elapsed.as_micros() as f32 / 1_000_000 as f32,
        proof_size,        
    )
    
}