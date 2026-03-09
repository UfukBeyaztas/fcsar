\name{fsac_pcaM}
\alias{fsac_pcaM}
\title{
Robust FPCA-based Fisher-consistent redescending estimation for spatial scalar-on-function regression
}
\description{
Fits a spatial scalar-on-function regression model by combining robust functional
principal component analysis (RFPCA) with a Fisher-consistent redescending
M-estimation step for the finite-dimensional spatial autoregressive model.
The function provides Andrews-loss and Danish-loss versions of the estimator.
}
\usage{
fsac_pcaM(y, x, nbasis = NULL, gpx = NULL, wei_mat,
          Type = c("Andrews", "Danish"),
          maxIter = 100, tol = 1e-6)
}
\arguments{
  \item{y}{
  A numeric vector of length \code{n} containing the scalar response values.
  }

  \item{x}{
  A numeric matrix of dimension \code{n x p_x} containing discretized observations
  of the functional predictor. Each row corresponds to one functional observation,
  and each column corresponds to one grid point on the functional support.
  }

  \item{nbasis}{
  An optional integer specifying the number of basis functions used in the RFPCA
  decomposition. If \code{NULL}, the default is \code{round(min(n/4, 40))}.
  }

  \item{gpx}{
  An optional numeric vector giving the grid points at which the functional predictor
  is observed. If \code{NULL}, an equally spaced grid over \eqn{(0,1)} is used.
  }

  \item{wei_mat}{
  A numeric \code{n x n} spatial weight matrix. It is typically row-standardized and
  should have zero diagonal.
  }

  \item{Type}{
  A character string specifying the redescending loss function used in the robust
  spatial estimation step. Possible values are \code{"Andrews"} and \code{"Danish"}.
  }

  \item{maxIter}{
  Maximum number of iterations allowed in the robust spatial M-estimation algorithm.
  }

  \item{tol}{
  A positive numeric convergence tolerance for the iterative estimation algorithm.
  }
}
\details{
The function fits the spatial scalar-on-function regression model
\deqn{
Y_i = \beta_0 + \rho \sum_{j=1}^n w_{ij}Y_j +
\int_{\mathcal{I}} X_i(t)\beta(t)\,dt + \varepsilon_i,
\qquad i=1,\dots,n,
}
where \eqn{Y_i} is a scalar response, \eqn{X_i(t)} is a functional predictor,
\eqn{\beta_0} is an intercept, \eqn{\rho} is the spatial autoregressive parameter,
\eqn{\mathbf{W} = (w_{ij})} is the spatial weight matrix, and \eqn{\beta(t)} is an
unknown coefficient function.

Since the predictor is functional, the model is first reduced to a
finite-dimensional form by robust functional principal component analysis.
If the functional predictor is approximated by
\deqn{
X_i(t) \approx \mu_R(t) + \sum_{k=1}^{K}\xi_{ik}\phi_k(t),
}
and the coefficient function is approximated by
\deqn{
\beta(t) \approx \sum_{k=1}^{K} b_k \phi_k(t),
}
then
\deqn{
\int_{\mathcal{I}} X_i(t)\beta(t)\,dt \approx \sum_{k=1}^{K}\xi_{ik}b_k.
}
Therefore, after RFPCA, the model becomes the finite-dimensional spatial
autoregressive regression
\deqn{
\mathbf{Y} = \rho \mathbf{W}\mathbf{Y} + \mathbf{Z}\boldsymbol{\gamma} +
\boldsymbol{\varepsilon},
}
where \eqn{\mathbf{Z} = [\mathbf{1}_n,\mathbf{\Xi}]} contains the intercept and the
retained robust principal component scores, and
\eqn{\boldsymbol{\gamma} = (\beta_0,b_1,\dots,b_K)^\top}.

The function first calls \code{getPCA()} with \code{method.type = "robust"} to obtain:
\itemize{
  \item robust principal component scores,
  \item robust principal directions,
  \item basis evaluation matrices needed to reconstruct the coefficient function.
}

In the second stage, the finite-dimensional spatial model is fitted by
\code{robestM()}, which implements a Fisher-consistent redescending M-estimation
procedure for spatial autoregressive regression. In particular, the estimator is
constructed from a bias-corrected estimating system so that the resulting estimating
equations are centered at the true parameter under the model. This improves both
robustness and theoretical coherence relative to ordinary likelihood-based or
monotone-loss M-estimation procedures.

Two redescending loss functions are available:
\itemize{
  \item \code{Type = "Andrews"}: Andrews sine loss,
  \item \code{Type = "Danish"}: Danish redescending loss.
}
These correspond to the Andrews-based and Danish-based versions of the proposed
FCSAR estimator in the associated methodology. In applications, the Andrews version
is often more aggressive in rejecting extreme outliers, while the Danish version may
provide a somewhat smoother redescending behavior.

After estimation, the coefficient function is reconstructed on the grid \code{gpx} as
\deqn{
\widehat{\beta}(t) = \sum_{k=1}^{K}\widehat{b}_k \widehat{\phi}_k(t),
}
and fitted values are computed from the reduced-form spatial representation
\deqn{
\widehat{\mathbf{Y}} =
(\mathbf{I}_n-\widehat{\rho}\mathbf{W})^{-1}
(\widehat{\beta}_0\mathbf{1}_n + \widehat{\mathbf{\Xi}}\widehat{\mathbf{b}}).
}

Thus, this function provides a robust FPCA-based implementation of spatial
scalar-on-function regression in which both functional outliers and scalar-response
outliers are handled through a two-layer robustification strategy:
first at the functional dimension-reduction stage, and second at the spatial
regression stage.
}
\value{
A list with the following components:
\item{b}{
A numeric vector containing the estimated coefficients associated with the retained
robust principal component scores.
}
\item{b0}{
The estimated intercept.
}
\item{bhat}{
A numeric vector containing the estimated coefficient function evaluated on the grid
specified by \code{gpx}.
}
\item{rho}{
The estimated spatial autoregressive parameter.
}
\item{sig}{
The estimated scale parameter.
}
\item{fitted.values}{
A numeric vector of fitted response values.
}
\item{residuals}{
A numeric vector of model residuals.
}
\item{fpca}{
The object returned by \code{getPCA()}, containing the robust FPCA decomposition
used in the model fitting step.
}
\item{ncomp}{
The number of retained robust principal components.
}
}
\author{
Muge Mutis, Ufuk Beyaztas, Han Lin Shang
}
\examples{
\donttest{
set.seed(123)

## Generate contaminated data
dat <- data_generation(
  n = 100,
  j = 101,
  rho = 0.5,
  sig.e = 1
)

## Andrews-loss version
fitA <- fsac_pcaM(
  y = dat$y,
  x = dat$x,
  wei_mat = dat$w,
  Type = "Andrews"
)

## Danish-loss version
fitD <- fsac_pcaM(
  y = dat$y,
  x = dat$x,
  wei_mat = dat$w,
  Type = "Danish"
)

plot(fitA$bhat, type = "l",
     main = "Estimated coefficient function",
     xlab = "Grid index", ylab = expression(hat(beta)(t)))
}
}
