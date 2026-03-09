\name{data_generation}
\alias{data_generation}
\title{Generate synthetic data for spatial scalar-on-function regression}
\description{
Generates synthetic datasets from a spatial scalar-on-function regression model
under clean or contaminated settings. The function returns the scalar response,
functional predictor observations on a grid, the row-normalized spatial weight
matrix, the true coefficient function, and the indices of contaminated
observations when contamination is introduced.
}
\usage{
data_generation(n, j = 101, rho, sig.e, m.o, out.p = 0)
}
\arguments{
  \item{n}{
  Integer. Training sample size.
  }
  \item{j}{
  Integer. Number of equally spaced grid points on the function support
  \eqn{[0,1]}. Default is \code{101}.
  }
  \item{rho}{
  Numeric. Spatial autoregressive parameter.
  }
  \item{sig.e}{
  Numeric. Standard deviation of the clean error term.
  }
  \item{m.o}{
  Numeric. Mean of the contamination error distribution used for generating
  vertical outliers. This argument is only used when \code{out.p > 0}.
  }
  \item{out.p}{
  Numeric in \eqn{[0,1]}. Proportion of contaminated observations. Default is
  \code{0}, corresponding to the clean-data setting.
  }
}
\details{
The function generates functional predictors from the expansion
\deqn{
\mathcal{X}_i(t) = \sum_{k=1}^{5} \kappa_{ik}\phi_k(t),
}
where \eqn{\kappa_{ik} \sim N(0, 4k^{-3/2})} and
\eqn{\phi_k(t) = \sin(k\pi t) - \cos(k\pi t)}.
The true slope function is
\deqn{
\beta(t) = \sin(2\pi t).
}

Let \eqn{m_i = \int_0^1 \mathcal{X}_i(t)\beta(t)\,dt}. The scalar response is
generated from the spatial scalar-on-function regression model
\deqn{
\mathbf{Y} = (\mathbf{I}_n - \rho \mathbf{W})^{-1}(\mathbf{m} + \boldsymbol{\varepsilon}),
}
where \eqn{\mathbf{W}} is a row-normalized inverse-distance spatial weight
matrix constructed on a one-dimensional regular grid.

When \code{out.p > 0}, a proportion of observations is randomly selected and
contaminated in two ways: the functional predictors are replaced by shifted
curves generated with coefficient means equal to one, and the corresponding
errors are replaced by draws from a contaminated normal distribution with mean
\code{m.o} and variance one. This produces both leverage-type functional
outliers and vertical outliers in the scalar response.
}
\value{
A list with the following components:
\item{y}{
A numeric vector containing the generated scalar response values.
}
\item{x}{
A matrix of dimension \code{n x j} containing the discretized functional
predictors.
}
\item{w}{
The \code{n x n} row-normalized spatial weight matrix.
}
\item{tcoefs}{
A vector containing the true coefficient function evaluated on the grid.
}
\item{out.index}{
An integer vector containing the indices of contaminated observations. Returns
\code{NULL} when \code{out.p = 0}.
}
}
\author{
Muge Mutis, Ufuk Beyaztas, and Han Lin Shang
}
\examples{
\donttest{
set.seed(123)

dat <- data_generation(
  n = 100,
  j = 101,
  rho = 0.5,
  sig.e = 1,
  m.o = 10,
  out.p = 0.05
)
}
}
