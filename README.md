### <span style="font-family:Sans Serif">Cougar</span>: Cubic Root Verifier Inner Product Argument under Discrete Logarithm Assumption

- Hyeonbum Lee, Seunghun Paik, Hyunjung Son, and Jae Hong Seo
- https://eprint.iacr.org/2024/616

### Introduction.

This repository is our Rust implementation of `Cougar` and other prior IPAs based on the DL assumption: `BulletProofs` and `Leopard`.

### Experimental Environment

Every experiment was done in the following experimental environment.

- OS: Ubuntu 20.04 LTS
- CPU: AMD EPYC 7543P (Single Core & Single threaded; Core speed: 2.8GHz)
- RAM: 512GB

We also attached raw data for reproducing Figure 7 and Table 3. Each data is recorded in order (Prover time, Verifier time, proof size) for logN = 10 to 20.

### Description

The overall structure of this codebase is as follows:

├ halo2curves          # Forked from "https://github.com/privacy-scaling-explorations/halo2curves", version 0.6.1,
│                        with slight modifications on the target group of BN254 (More functionalities & MSM)
├ halo2_proofs         # Forked from "https://github.com/zcash/halo2/tree/main/halo2_proofs", version 0.3.1,
│                        with slight modifications on dealing with polynomials
├ RawData              # Raw data obtained from our experimental environment, including BulletProofs, Leopard and Cougar.
├ src
    ├ BulletProofs     # Implementation of BulletProofs
    ├ Leopard          # Implementation of Leopard
    ├ Cougar           # Implementation of Cougar
    └ main.rs
└ plot_data.iypnb      # Visualization for Raw data

For more details on our implementations, please check `README.md` files in each implementation folder.

Enjoy!
