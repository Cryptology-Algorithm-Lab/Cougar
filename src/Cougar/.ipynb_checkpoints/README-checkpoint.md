### Implemenatation of Cougar

### Introduction

This is an implementation of the proposed Cougar using half-pairing cycle of curves: Grumpkin and BN254.

We used modified `halo2curves` and `halo2_proofs` libraries in the main root.

### Description

The overall structure of this codebase is as follows:

'''
├ LeopardPCS 
    ├ leopard_prover.rs      # Leopard PCS (Batched) Prover 
    └ leopard_verifier.rs    # Leopard PCS (Batched) Verifier
├ utils.rs                   # Utilities
├ custom_gate.rs             # Tools for constructing custom gates.
├ execution_trace.rs         # Tools for constructing execution trace.
├ permutation.rs             # Tools for permutation argument.
├ setup.rs                   # Setup
├ prover.rs                  # Proof Generation
├ verifier.rs                # Verification
├ transcript.rs              # Utility for handling transcript
├ test_IPA.rs                # Test Codes
├ LeopardPCS.rs
└ mod.rs
'''

Enjoy!