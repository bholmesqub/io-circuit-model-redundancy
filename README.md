# io-circuit-model-redundancy
A symbolic analysis method for detecting parameter redundancy in circuit models used in identification.

## Overview

The central contribution of this repository is the function `detectRedundancy.m`. This function performs symbolic analysis upon two vectors: `theta`, a model's parameters, and `kappa`, how the parameters appear in the model. The symbolic analysis determines how many of the parameters in `theta` can be estimated given the contents of `kappa` by checking the rank of the jacobian matrix.

This code is largely based on the paper *"Determining the parametric structure of models" - D. J. Cole et. al.,
2010, Mathematical Biosciences*. The Maple code at [D. Cole's site](https://www.kent.ac.uk/smsas/personal/djc24/index.htm) was also invaluable in development.

Much of the extended functionality of the Maple code has been excluded as the main objective of this function is only to *detect* the redundancy. The Maple code offers a method of determining what combinations can be estimated, but this requires a PDE solver which MATLAB does not have.
