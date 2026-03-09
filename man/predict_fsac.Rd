\name{predict_fsac}
\alias{predict_fsac}
\title{
Prediction for fitted spatial scalar-on-function regression models
}
\description{
Obtains predicted response values for new observations from a fitted
spatial scalar-on-function regression model. The function is designed for
objects returned by \code{fsac_pca()} and \code{fsac_pcaM()} and computes
out-of-sample predictions using the estimated spatial dependence parameter,
intercept, and coefficient vector associated with the retained functional
principal component scores.
}
\usage{
predict_fsac(object, xnew, wnew)
}
\arguments{
  \item{object}{
  A fitted model object returned by \code{\link{fsac_pca}} or
  \code{\link{fsac_pcaM}}.
  }

  \item{xnew}{
  A numeric matrix containing discretized new functional predictor observations.
  Each row corresponds to one new curve and each column corresponds to one grid
  point of the functional support.
  }

  \item{wnew}{
  A numeric spatial weight matrix for the new observations. This matrix is
  typically row-standardized and should have zero diagonal. Its dimension must
  match the number of rows of \code{xnew}.
  }
}
\details{
Let the fitted spatial scalar-on-function regression model be
\deqn{
Y_i = \beta_0 + \rho \sum_{j=1}^n w_{ij}Y_j + \int_{\mathcal{I}} X_i(t)\beta(t)\,dt + \varepsilon_i,
\qquad i=1,\dots,n,
}
where \eqn{Y_i} is the scalar response, \eqn{X_i(t)} is the functional predictor,
\eqn{\beta_0} is the intercept, \eqn{\rho} is the spatial autoregressive parameter,
and \eqn{\beta(t)} is the coefficient function.

After functional principal component decomposition, the model is approximated by
a finite-dimensional spatial autoregressive regression of the form
\deqn{
\mathbf{Y} = \rho \mathbf{W}\mathbf{Y} + \mathbf{Z}\boldsymbol{\gamma} + \boldsymbol{\varepsilon},
}
where \eqn{\mathbf{Z} = [\mathbf{1},\mathbf{\Xi}]} and \eqn{\mathbf{\Xi}} contains
the retained principal component scores. If \eqn{\widehat{\rho}},
\eqn{\widehat{\beta}_0}, and \eqn{\widehat{\mathbf{b}}} denote the fitted parameter
estimates, then for a new sample with functional predictor matrix \eqn{\mathbf{X}^{new}}
and spatial weight matrix \eqn{\mathbf{W}^{new}}, prediction is based on
\deqn{
\widehat{\mathbf{Y}}^{new}
=
(\mathbf{I}-\widehat{\rho}\mathbf{W}^{new})^{-1}
\left(
\widehat{\beta}_0\mathbf{1} + \mathbf{\Xi}^{new}\widehat{\mathbf{b}}
\right),
}
where \eqn{\mathbf{\Xi}^{new}} is the matrix of principal component scores for the
new functional observations.

The function first extracts the fitted model components from \code{object},
namely:
\itemize{
  \item the estimated spatial parameter \code{rho},
  \item the estimated intercept \code{b0},
  \item the estimated coefficient vector \code{b},
  \item the stored FPCA or RFPCA decomposition \code{fpca}.
}

The new functional data \code{xnew} are then projected onto the principal
component space determined from the training sample through the auxiliary function
\code{getPCA_test()}. This ensures that the new functional observations are
represented in the same finite-dimensional basis that was used when fitting the
original model.

Finally, the predicted values are computed through the reduced-form spatial
representation
\deqn{
\widehat{\mathbf{Y}}^{new}
=
(\mathbf{I}-\widehat{\rho}\mathbf{W}^{new})^{-1}
\widehat{\boldsymbol{\mu}}^{new},
}
with
\deqn{
\widehat{\boldsymbol{\mu}}^{new}
=
\widehat{\beta}_0\mathbf{1} + \mathbf{\Xi}^{new}\widehat{\mathbf{b}}.
}

This function is intended for genuine out-of-sample prediction. Therefore,
\code{xnew} and \code{wnew} should correspond to a new set of observations and
their associated spatial configuration. In applications, \code{wnew} may be
constructed by inverse-distance weights, neighborhood relations, or any other
user-specified spatial weighting rule that is appropriate for the new sample.

The prediction formula used here is fully compatible with both classical and
robust fits produced by \code{fsac_pca()} and with the Fisher-consistent
redescending fits produced by \code{fsac_pcaM()}.
}
\value{
A numeric vector or matrix of predicted response values for the new observations.
The length (or number of rows) equals the number of rows in \code{xnew}.
}
\author{
Muge Mutis, Ufuk Beyaztas, Han Lin Shang
}
\note{
The spatial weight matrix \code{wnew} must be conformable with \code{xnew}.
If \code{xnew} has \eqn{m} rows, then \code{wnew} must be an \eqn{m\times m}
matrix.

No uncertainty quantification is returned. The function only computes point
predictions.
}
\examples{
\donttest{
set.seed(123)

## Generate training data
dat_train <- data_generation(
  n = 100, j = 101, rho = 0.5, sig.e = 1,
  m.o = 10, out.p = 0.05
)

## Fit a robust model
fit <- fsac_pcaM(
  y = dat_train$y,
  x = dat_train$x,
  wei_mat = dat_train$w,
  Type = "Andrews"
)

## Generate test data
dat_test <- data_generation(
  n = 50, j = 101, rho = 0.5, sig.e = 1,
  m.o = 10, out.p = 0
)

## Obtain predictions
yhat <- predict_fsac(
  object = fit,
  xnew = dat_test$x,
  wnew = dat_test$w
)
}
}
