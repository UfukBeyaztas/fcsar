# fcsar <img src="https://img.shields.io/badge/R-%3E=3.5.0-1f425f.svg" alt="R (>= 3.5.0)" align="right" height="20"/>

**fcsar** provides tools for fitting **classical**, **robust**, and **Fisher-consistent redescending** estimators in **spatial scalar-on-function regression** models. The package combines functional principal component analysis with spatial autoregressive modeling, allowing users to analyze scalar outcomes that depend on spatially structured observations and functional predictors.

The package includes simulation-based data generation, FPCA- and RFPCA-based estimation, Fisher-consistent redescending robust estimation, and out-of-sample prediction under spatial dependence.

## 🚀 Key Features

- **Spatial scalar-on-function regression:** fits regression models where the response is scalar, the predictor is functional, and spatial dependence is present through a SAR structure.
- **Multiple estimation strategies:** includes
  - classical FPCA + likelihood-based estimation,
  - robust FPCA + robust spatial estimation,
  - Fisher-consistent redescending estimation based on **Andrews** and **Danish** loss functions.
- **Built-in simulation tools:** generates clean and contaminated datasets for Monte Carlo experiments under spatial dependence.
- **Robust functional dimension reduction:** supports both classical and robust FPCA-based representations of functional predictors.
- **Prediction for new spatial units:** provides out-of-sample prediction using the estimated spatial multiplier and the functional regression effect.
- **User manual included:** the package is accompanied by a user manual describing the methodology, arguments, returned objects, and example workflows.

---

## 📦 Installation

You can install the development version of **fcsar** from GitHub:

```r
install.packages("remotes")
remotes::install_github("UfukBeyaztas/fcsar")


## 💻 Quick Start Example
1. Generate a contaminated sample
library(fcsar)
set.seed(123)
dat <- data_generation(
  n = 100,
  j = 101,
  rho = 0.5,
  sig.e = 1,
  m.o = 10,
  out.p = 0.05
)

2. Fit the classical FPCA-based model
fit_classical <- fsac_pca(
  y = dat$y,
  x = dat$x,
  wei_mat = dat$w,
  method.type = "classical"
)

3. Fit the robust FPCA-based model
fit_robust <- fsac_pca(
  y = dat$y,
  x = dat$x,
  wei_mat = dat$w,
  method.type = "robust"
)

4. Fit the Fisher-consistent redescending model
Andrews-loss version
fit_andrews <- fsac_pcaM(
  y = dat$y,
  x = dat$x,
  wei_mat = dat$w,
  Type = "Andrews"
)

Danish-loss version
fit_danish <- fsac_pcaM(
  y = dat$y,
  x = dat$x,
  wei_mat = dat$w,
  Type = "Danish"
)

5. Inspect estimated quantities
fit_andrews$rho
fit_danish$rho
head(fit_andrews$fitted.values)

6. Predict for new spatial units
dat_test <- data_generation(
  n = 50,
  j = 101,
  rho = 0.5,
  sig.e = 1,
  m.o = 10,
  out.p = 0
)

preds <- predict_fsac(
  object = fit_andrews,
  xnew = dat_test$x,
  wnew = dat_test$w
)
