# Distributed Machine Learning for Localization in Cell-Free Massive MIMO Systems

This repository contains the code associated with the paper:  
**"Distributed Machine Learning Approach for Low-Latency Localization in Cell-Free Massive MIMO Systems."**

---

## Overview

The simulations were executed in **MATLAB**, with results saved using `diary` logs for reproducibility. Most figures were generated using **Python**, which offers greater flexibility and control for publication-quality plots.  
- **Exception**:  
  - **Fig. 9** (CDF of localization error) was generated directly in MATLAB.  
  - **Fig. 10** was derived from CRLB outputs printed by `Fig_3_to_6_localization_across_N.m`.

To ensure reproducibility, please **set the random number generator (`rng`) seed** to the value used in this repository.

> **Note**: Full simulation runs can take several days. For convenience, **final results are included in the repository**.

---

## Generating Figure 9 (CDF of Localization Errors)

1. **Run** `Fig9_CDF_Localization_error.m` to compute localization error data.
2. **Then run** `Fig9_CDF_Localization_error_plot.m` to generate the CDF plot.  
   Be sure to **update the variable name and plot title** to match the localization method being visualized.
3. **Alternatively**, use the precomputed file `Fig9_CDF_Localization_error_plot.mat` to generate the plot directly without recomputation. This `.mat` file includes CDF data for **all methods**.

---

## Additional Methods (Not Discussed in Paper)

The repository includes code for a **Distributed Mean + Z-Score filtering method** (not covered in the paper).  
This method filters out outliers using z-score thresholding **before** averaging the position estimates. It achieves:
- **Lower localization error** than the simple mean,
- **Slightly higher error ellipse area** than the mean,
- But **still significantly lower ellipse area** than the median.

---

## FCNN Baseline Notes

The referenced FCNN baseline paper specifies only the **number of layers**, not the full architecture.  
Therefore, the following hyperparameters were selected empirically based on performance at **K = 225** training points:
- Number of nodes per layer
- Activation functions
- Optimizer
- Number of training epochs

---

## Shadowing Model with Spatial Correlation

In our setup, shadowing effects are modeled using a **spatially correlated Gaussian process** to reflect realistic wireless environments.

- The correlation coefficient between two points is defined as:  
  \[
  \rho(d) = 2^{-\frac{d}{\text{decorr}}}
  \]
  where:
  - *d* is the distance between points
  - *decorr* is a decorrelation parameter (scenario-specific, per 3GPP)

- A **correlation matrix** is built incrementally:
  - Each new reference point (RP) computes distances to previous RPs
  - Conditional **mean** and **standard deviation** are calculated using Theorem 10.2 from *Steven Kay's Estimation Theory* (p. 325 — snippet included in the repo)
  - Shadowing terms are generated accordingly using **conditional Gaussian statistics**

- During the **online phase**, the same model is used to compute shadowing at test points (TPs), maintaining spatial consistency.

For a deeper understanding, see:
- [`shadowing_realisation_example.m`](shadowing_realisation_example.m) – illustrates the shadowing process and is suitable for experimentation
- [arXiv:2108.02541](https://arxiv.org/pdf/2108.02541) – especially page 66
- MATLAB implementation of the full procedure is included in the repo

---

## License

This repository is intended for academic and research use only. Please cite the original paper if you use the code.

---

For questions or contributions, feel free to open an issue or submit a pull request.
