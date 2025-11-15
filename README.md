# Total Variation Regularized Polynomial Chaos Expansion (TV-PCE)

This repository contains the MATLAB code used in **Paper D of my PhD thesis**, where I develop a Total Variation–regularized Polynomial Chaos Expansion (TV-PCE) method to improve uncertainty quantification for systems with discontinuities. Standard Polynomial Chaos Expansion (PCE) suffers from Gibbs oscillations when the model response is discontinuous. The TV-PCE method suppresses these oscillations by applying total variation regularization to the PCE coefficients.

The code in this repository reproduces the numerical experiments and examples presented in the paper, including the stochastic McKendrick–von Foerster equation with discontinuous mortality.

---

## How to Run the MATLAB Code

1. Open MATLAB and navigate to the project directory:
    ```matlab
    cd /path/to/paper-D-tv-pce
    ```

2. Run the main script for generating PCE approximations:
    ```matlab
    run('main.m')
    ```

3. To reproduce the numerical examples from the paper, run the corresponding example scripts:
    ```matlab
    run('example1.m')
    run('example2.m')
    run('example3.m')
    run('example4.m')
    ```

---

## Purpose of the Code

The scripts implement:
- Standard Polynomial Chaos Expansion (PCE)
- Total Variation regularization applied to PCE coefficients
- Iterative solution of the variational TV-PCE problem
- Numerical solution of the stochastic McKendrick–von Foerster PDE
- Reproduction of all examples and figures from the paper

---

## Citation

If you use this code, please cite:

Jokar, M. (2025).  
"Total Variation Method to Improve the Convergence of Polynomial Chaos Expansion Approximations in Discontinuities."  
PhD Thesis, Western Norway University of Applied Sciences (HVL).

---

## Contact

Mahmood Jokar  
PhD Researcher, HVL / Institute of Marine Research  
Email: jmah@hvl.no
