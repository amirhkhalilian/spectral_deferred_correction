# Spectral Deferred Correction for Non-stiff ODEs

This project implements the **Spectral Deferred Correction (SDC)** method for solving non-stiff ordinary differential equations (ODEs). The method enhances the accuracy of low-order solvers using deferred corrections over Chebyshev nodes [DGR00], achieving high-order accuracy in a computationally efficient manner.

## 📘 Overview

Traditional solvers like Runge-Kutta face stability or cost limitations at high orders. Spectral Deferred Correction (SDC), introduced in [DGR00], is an iterative scheme that accelerates low-order methods through a correction sweep. The details and background are explained in sdc_method.pdf.

This repo includes:
- An implementation of the SDC method with **fixed and adaptive step-size**.
- Utilities for **spectral integration** using Chebyshev polynomials and FFT.
- Numerical experiments demonstrating convergence and accuracy.

## 📈 Results

- **Empirical convergence** up to order ~J+1, where J is the number of correction sweeps.
- **Adaptive step-size** improves efficiency by dynamically adjusting interval size based on error estimates.
- Spectral methods outperform basic Euler in accuracy and stability for the same number of function evaluations.

## 🛠 Requirements

- MATLAB (tested with R2021a and later)
- Symbolic Math Toolbox (optional for analytic comparisons)

## 🚀 Getting Started

1. Run a demo:
    ```matlab
    % Fixed step-size demo
    Demo_SDC_IVP;

    % Adaptive step-size demo
    Demo_SDCA_ellipj;
    ```

2. Modify `FE_SDC_solver.m` for custom ODEs or options.

### ⚙️ Solver Options

You can configure the solver via the `opts` struct:

| Option       | Description                            | Default |
|--------------|----------------------------------------|---------|
| `k`          | Initial step size                      | 0.1     |
| `m`          | Number of Chebyshev nodes per interval | 5       |
| `J`          | Number of correction sweeps            | 4       |
| `adaptive`   | Enable adaptive step-size (true/false) | false   |
| `tol`        | Tolerance for adaptive control         | NaN     |
| `verbose`    | Print debugging info                   | false   |

## 📖 Reference

[DGR00] Dutt, A., Greengard, L., & Rokhlin, V. (2000). *Spectral deferred correction methods for ordinary differential equations*. BIT Numerical Mathematics, 40(2), 241–266.





