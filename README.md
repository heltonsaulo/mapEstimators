# Exponential-Family Transformation Analysis

This project implements:

1. **Density and estimation functions** for a transformed exponential-family model  
2. **Sampling** from various special cases (Gamma, Inverse Gamma, Weibull, Inverse Weibull)  
3. **Analytic** and **numeric MAP** estimators for parameters \(\mu\) and \(\sigma\)  
4. **Maximum likelihood** estimators for the Gamma case  
5. A **Monte Carlo study** comparing relative bias and MSE of the proposed analytic estimator, MAP and ML

---

## Structure

- `R/functions.R`  
  - Core density function  
  - General analytic estimator  
  - Model-specific wrappers  
  - RNG for each special case  
  - Numeric MAP estimator  
  - MLE for the Gamma model

- `R/simulation.R`  
  - Runs the Monte Carlo comparison  
  - Produces two ggplot2 figures: relative bias and MSE

---

## Usage

1. **Install dependencies**  
   ```r
   install.packages(c("ggplot2","tidyverse","dplyr","MASS","nleqslv"))
