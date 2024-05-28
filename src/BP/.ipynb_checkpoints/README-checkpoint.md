### Implemenatation of BulletProofs

### Introduction

This is an implementation of BulletProofs [BBB+18] using Secp256k1 Curve.

We used a forked `halo2curve` library in the main root.

### Description

The overall structure of this codebase is as follows:

'''
├ utils.rs         # Utilities
├ setup.rs         # Setup
├ prover.rs        # Proof Generation
├ verifier.rs      # Verification
├ transcript.rs    # Utility for handling transcript
├ test_IPA.rs      # Test Codes
└ mod.rs
'''

Enjoy!