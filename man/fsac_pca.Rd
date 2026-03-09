\name{fsac_pca}
\alias{fsac_pca}
\title{
FPCA-based estimation for spatial scalar-on-function regression
}
\description{
Fits a spatial scalar-on-function regression model using a functional principal
component analysis (FPCA) representation of the functional predictor. The function
implements both a classical estimation approach and a robust alternative based on
robust FPCA and robust spatial autocorrelation estimation.
}
\usage{
fsac_pca(y, x, nbasis = NULL, gpx = NULL, wei_mat,
         method.type = c("classical", "robust"),
         c1 = 1.4, c2 = 2.4, c3 = 1.65, tol = 1e-04)
}
\arguments{
  \item{y}{
  A numeric vector of length \code{n} containing the scalar response values.
  }

  \item{x}{
  A numeric matrix of dimension \code{n x p_x} containing the discretized functional
  predictor values. Each row corresponds to one observation curve, and each column
  corresponds to one grid point on the functional support.
  }

  \item{nbasis}{
  An optional integer specifying the number of basis functions used in the FPCA
  approximation. If \code{NULL}, the default is \code{round(min(n/4, 40))}.
  }

  \item{gpx}{
  An optional numeric vector giving the grid points at which the functional predictor
  is observed. If \code{NULL}, an equally spaced grid over \eqn{(0,1)} is used.
  }

  \item{wei_mat}{
  A numeric \code{n x n} spatial weight matrix. Typically, this matrix should be
  row-standardized and have zero diagonal.
  }

  \item{method.type}{
  A character string specifying the estimation method. Possible values are
  \code{"classical"} for the classical FPCA + likelihood-based spatial estimator,
  and \code{"robust"} for the robust FPCA + robust spatial estimator.
  }

  \item{c1}{
  A tuning constant used in the weight function \code{alpha_fun()} within the robust
  estimation step. Ignored when \code{method.type = "classical"}.
  }

  \item{c2}{
  A tuning constant used in the scale update through \code{psi_fun()} and
  \code{hfun()} in the robust estimation step. Ignored when
  \code{method.type = "classical"}.
  }

  \item{c3}{
  A tuning constant used in the robust objective function for updating the spatial
  dependence parameter. Ignored when \code{method.type = "classical"}.
  }

  \item{tol}{
  A positive numeric convergence tolerance for the iterative robust estimation
  algorithm.
  }
}
\details{
The function fits a spatial scalar-on-function regression model of the form
\deqn{
Y_i = \beta_0 + \rho \sum_{j=1}^n w_{ij}Y_j + \int_{\mathcal{I}} X_i(t)\beta(t)\,dt + \varepsilon_i,
\qquad i=1,\dots,n,
}
or, in matrix notation,
\deqn{
\mathbf{Y} = \rho \mathbf{W}\mathbf{Y} + \beta_0 \mathbf{1}_n +
\int_{\mathcal{I}} \mathbf{X}(t)\beta(t)\,dt + \boldsymbol{\varepsilon},
}
where \eqn{\mathbf{W} = (w_{ij})} is the spatial weight matrix, \eqn{\rho} is the
spatial autoregressive parameter, and \eqn{\beta(t)} is the unknown coefficient
function.

Because the predictor is functional, the model is estimated after dimension reduction
through functional principal component analysis. If
\deqn{
X_i(t) \approx \mu(t) + \sum_{k=1}^{K} \xi_{ik}\phi_k(t),
}
then the slope function is approximated by
\deqn{
\beta(t) \approx \sum_{k=1}^{K} b_k \phi_k(t),
}
so that
\deqn{
\int_{\mathcal{I}} X_i(t)\beta(t)\,dt \approx \sum_{k=1}^{K}\xi_{ik}b_k.
}
Hence, after FPCA, the problem reduces to a finite-dimensional spatial autoregressive
regression model:
\deqn{
\mathbf{Y} = \rho \mathbf{W}\mathbf{Y} + \mathbf{Z}\boldsymbol{\gamma} + \boldsymbol{\varepsilon},
}
where \eqn{\mathbf{Z}} contains an intercept and the retained FPCA scores, and
\eqn{\boldsymbol{\gamma}} contains the corresponding regression coefficients.

The function first computes a functional principal component decomposition by calling
\code{getPCA()}, which returns the retained principal component basis, scores, and
basis evaluation matrices.

\strong{Classical estimation.}
If \code{method.type = "classical"}, the finite-dimensional spatial regression model
is fitted by calling \code{sac_ml()}, which returns estimates of the intercept,
component coefficients, spatial dependence parameter, and error variance. This
corresponds to the standard nonrobust FPCA-based approach.

\strong{Robust estimation.}
If \code{method.type = "robust"}, the function starts from the classical estimator
and then iteratively updates the regression coefficients, scale parameter, and spatial
autoregressive parameter. The robust procedure uses:
\itemize{
  \item weighted least squares updates for the regression coefficients,
  \item a robust scale update based on \code{psi_fun()} and \code{hfun()},
  \item numerical optimization over \eqn{\rho} using a robust objective function.
}
This produces a robust FPCA-based spatial estimator that is less sensitive to
outlying curves and anomalous response values than the classical approach.

After estimation, the coefficient function is reconstructed on the observation grid as
\deqn{
\widehat{\beta}(t) = \sum_{k=1}^{K}\widehat{b}_k\widehat{\phi}_k(t),
}
and fitted values are obtained from the reduced-form representation
\deqn{
\widehat{\mathbf{Y}} =
(\mathbf{I}_n-\widehat{\rho}\mathbf{W})^{-1}(\widehat{\beta}_0\mathbf{1}_n +
\widehat{\mathbf{\Xi}}\widehat{\mathbf{b}}).
}
}
\value{
A list with the following components:
\item{b}{
A numeric vector containing the estimated coefficients associated with the retained
principal component scores.
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
The estimated scale (standard deviation) parameter.
}
\item{fitted.values}{
A numeric vector of fitted response values.
}
\item{residuals}{
A numeric vector of model residuals.
}
\item{fpca}{
The object returned by \code{getPCA()}, containing the FPCA decomposition used in the
model fit.
}
\item{ncomp}{
The number of retained principal components.
}
}
\author{
Muge Mutis, Ufuk Beyaztas, Han Lin Shang
}
\note{
The argument \code{method.type = "robust"} performs robust estimation after FPCA-based
dimension reduction. The precise behavior of the robust estimator depends on the
definitions of the internal functions \code{alpha_fun()}, \code{psi_fun()}, and
\code{hfun()}.

The returned component \code{bhat} contains the estimated slope function evaluated on
the observation grid, not an \code{fdata} object.
}
\examples{
\donttest{
set.seed(123)

## Generate a clean sample
dat <- data_generation(
  n = 100,
  j = 101,
  rho = 0.5,
  sig.e = 1
)

## Classical FPCA-based fit
fit_classical <- fsac_pca(
  y = dat$y,
  x = dat$x,
  wei_mat = dat$w,
  method.type = "classical"
)

## Robust FPCA-based fit
fit_robust <- fsac_pca(
  y = dat$y,
  x = dat$x,
  wei_mat = dat$w,
  method.type = "robust"
)

fit_classical$rho
fit_robust$rho

## Plot estimated coefficient function
plot(fit_robust$bhat, type = "l",
     main = "Estimated coefficient function",
     xlab = "Grid index", ylab = expression(hat(beta)(t)))
}
}
