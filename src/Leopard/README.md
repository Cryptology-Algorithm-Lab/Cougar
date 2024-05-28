### Implemenatation of Leopard

### Introduction

This is an implementation of Leopard [KLS22, KLLS23] using BLS12-381 Curve.

We used `blstrs` library of version 0.7.1.

### Description

The overall structure of this codebase is as follows:

'''
├ utils.rs            # Utilities.
├ setup.rs            # Setup.
├ IPA.rs              # Code for doing IPA.
├ msm_pippenger.rs    # Implemenation of Pippenger algorithm for the target group.
├ transcript.rs       # Utility for handling transcript.
├ test_IPA.rs         # Test Code.
└ mod.rs
'''

Enjoy!