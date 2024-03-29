% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_exv.R
\name{Exv.VXX}
\alias{Exv.VXX}
\alias{Exv.VES}
\alias{Exv.VRS}
\alias{Exv.VER}
\alias{Exv.VRR}
\alias{Var.VES}
\alias{Var.VRS}
\alias{Var.VER}
\alias{Var.VRR}
\title{Moments of eigenvalue dispersion indices}
\usage{
Exv.VES(
  Sigma,
  n = 100,
  Lambda,
  divisor = c("UB", "ML"),
  m = switch(divisor, UB = n, ML = n + 1),
  drop_0 = FALSE,
  tol = .Machine$double.eps * 100,
  ...
)

Exv.VRS(
  Sigma,
  n = 100,
  Lambda,
  drop_0 = FALSE,
  tol = .Machine$double.eps * 100,
  ...
)

Exv.VER(
  Rho,
  n = 100,
  Lambda,
  tol = .Machine$double.eps * 100,
  tol.hg = 0,
  maxiter.hg = 2000,
  ...
)

Exv.VRR(
  Rho,
  n = 100,
  Lambda,
  tol = .Machine$double.eps * 100,
  tol.hg = 0,
  maxiter.hg = 2000,
  ...
)

Var.VES(
  Sigma,
  n = 100,
  Lambda,
  divisor = c("UB", "ML"),
  m = switch(divisor, UB = n, ML = n + 1),
  drop_0 = FALSE,
  tol = .Machine$double.eps * 100,
  ...
)

Var.VRS(
  Sigma,
  n = 100,
  Lambda,
  drop_0 = FALSE,
  tol = .Machine$double.eps * 100,
  ...
)

Var.VER(Rho, n = 100, Lambda, ...)

Var.VRR(
  Rho,
  n = 100,
  method = c("Pan-Frank", "Konishi"),
  Lambda,
  fun = c("pfd", "pfv", "pfc", "pf", "klv", "kl", "krv", "kr"),
  ...
)
}
\arguments{
\item{Sigma}{Population covariance matrix; assumed to be validly constructed.}

\item{n}{Degrees of freedom (not sample sizes); numeric of length 1 or more.}

\item{Lambda}{Numeric vector of population eigenvalues.}

\item{divisor}{Either \code{"UB"} (default) or \code{"ML"},
to decide the default value of \code{m}.}

\item{m}{Divisor for the sample covariance matrix (\eqn{n*} in Watanabe (2022)).
By default equals \eqn{n}.}

\item{drop_0}{Logical, when \code{TRUE}, eigenvalues smaller than \code{tol} are dropped.}

\item{tol}{For covariance-related functions, this is the tolerance/threshold
to be used with drop_0.  For correlation-related functions,
this is passed to \code{Exv.r2()} along with other arguments.}

\item{...}{In \code{Var.VRR()}, additional arguments are passed to an internal
function which it in turn calls.  Otherwise ignored.}

\item{Rho}{Population correlation matrix; assumed to be validly constructed
(although simple checks are done).}

\item{tol.hg, maxiter.hg}{Passed to \code{Exv.r2()}; see description of that function.}

\item{method}{For \code{Var.VRR()} (and \code{Var.VER()}), determines the method
to obtain approximate variance in non-null conditions.  Either
\code{"Pan-Frank"} (default) or \code{"Konishi"}.  See \dQuote{Details}.}

\item{fun}{For \code{Var.VRR()} (and \code{Var.VER()}), determines the function
to be used to evaluate approximate variance.  See
\dQuote{Details}.  Options allowed are: \code{"pfd"}, \code{"pfv"},
\code{"pfc"}, \code{"pf"}, \code{"klv"}, \code{"kl"}, \code{"krv"},
and \code{"kr"}.}
}
\value{
A numeric vector of the desired moment, corresponding to \code{n}.
}
\description{
Functions to calculate expectation/variance of (relative) eigenvalue variance
of sample covariance/correlation matrices for a given population
covariance/correlation matrix and degrees of freedom \eqn{n}

\code{Exv.VES()}: expectation of eigenvalue variance of covariance matrix
\eqn{\mathrm{E}[V(\mathbf{S})]}{E[V(S)]}

\code{Exv.VRS()}: expectation of relative eigenvalue variance of
covariance matrix \eqn{\mathrm{E}[V_{\mathrm{rel}}(\mathbf{S})]}{E[Vrel(S)]}

\code{Exv.VER()}: expectation of eigenvalue variance of correlation
matrix \eqn{\mathrm{E}[V(\mathbf{R})]}{E[V(R)]}

\code{Exv.VRR()}: expectation of relative eigenvalue variance of
correlation matrix
\eqn{\mathrm{E}[V_{\mathrm{rel}}(\mathbf{R})]}{E[Vrel(R)]}

\code{Var.VES()}: variance of eigenvalue variance of
covariance matrix \eqn{\mathrm{Var}[V(\mathbf{S})]}{Var[V(S)]}

\code{Var.VRS()}: variance of relative eigenvalue variance of
covariance matrix
\eqn{\mathrm{Var}[V_{\mathrm{rel}}(\mathbf{S})]}{Var[Vrel(S)]}

\code{Var.VER()}: variance of eigenvalue variance of
correlation matrix \eqn{\mathrm{Var}[V(\mathbf{R})]}{Var[V(R)]}

\code{Var.VRR()}: variance of relative eigenvalue variance of
correlation matrix
\eqn{\mathrm{Var}[V_{\mathrm{rel}}(\mathbf{R})]}{Var[Vrel(R)]}
}
\details{
\code{Exv.VES()}, \code{Var.VES()}, and \code{Exv.VRR()} return exact
moments.  \code{Exv.VRS()} and \code{Var.VRS()} return approximations based
on the delta method, except under the null condition
(\eqn{\mathbf{\Sigma}}{\Sigma} proportional to the identity matrix)
where exact moments are returned.

\code{Var.VRR()} returns the exact variance when \eqn{p = 2} or
under the null condition (\eqn{\mathbf{P}}{\Rho} is the identity
matrix).  Otherwise, asymptotic variance is calculated with the \code{method}
of choice: either \code{"Pan-Frank"} (default) or \code{"Konishi"}.  They
correspond to Pan and Frank's heuristic approximation and Konishi's
asymptotic theory, respectively (see Watanabe, 2022).

In this case, calculations are handled by one of the internal functions
\code{AVar.VRR_xx()} (\code{xx} is a suffix to specify \R
implementation).  For completeness, it is possible to directly specify the
function to be used with the argument \code{fun}: for the Pan--Frank method,
\code{"pfd"} (default), \code{"pfv"}, and \code{"pf"}; and for the
Konishi method, \code{"klv"} (default), \code{"kl"}, \code{"krv"}, and
\code{"kr"}.  Within each group, these function yield identical results but
differ in speed (the defaults are the fastest).  See
\code{\link{AVar.VRR_xx}} for details of these functions.

The Pan--Frank method takes a substantial amount of time to be executed
when \eqn{p} is large.  Several \proglang{C++} functions are provided in
the extension package
\href{https://github.com/watanabe-j/eigvaldispRcpp}{\pkg{eigvaldispRcpp}} to
speed-up the calculation.  When this package is available, the default
\code{fun} for the Pan--Frank method is set to \code{"pfc"}.  The
option for \proglang{C++} function is controlled by the argument
\code{cppfun} which in turn is passed to \code{AVar.VRR_pfc()}
(see \code{\link{AVar.VRR_xx}}).  When this argument is provided, the
argument \code{fun} is ignored with a warning, unless one of the Konishi
methods is used (in which case \code{cppfun} is ignored with a warning).

Since the eigenvalue variance of a correlation matrix
\eqn{V(\mathbf{R})}{V(R)} is simply \eqn{(p - 1)} times the
relative eigenvalue variance \eqn{V_{\mathrm{rel}}(\mathbf{R})}{Vrel(R)} of
the same matrix, their distributions are identical up to this
scaling.  Hence, \code{Exv.VER()} and \code{Var.VER()} calls \code{Exv.VRR()}
and \code{Var.VRR()} (respectively), whose outputs are scaled and
returned.  These functions are provided for completeness, although
there will be little practical demand for these functions (and
\eqn{V(\mathbf{R})}{V(R)} itself).

As detailed in Watanabe (2022), the distribution of
\eqn{V_{\mathrm{rel}}(\mathbf{R})}{Vrel(R)} cannot be
uniquely specified by eigenvalues alone.  Hence, a full correlation
matrix \code{Rho} should preferably be provided.  Otherwise, a correlation
matrix is constructed from the eigenvalues \code{Lambda} provided using
the function \code{GenCov()} with randomly picked eigenvectors.

On the other hand, the choice of eigenvectors does not matter for
covariance matrices, thus either the full covariance matrix \code{Sigma} or
vector of eigenvalues \code{Lambda} can be provided (although the latter
is slightly faster if available).

When \code{Rho} is provided, some simple checks are done: the matrix is
scaled to have diagonals of 1; and if any of these are unequal,
an error is returned.

These moments are derived under the assumption of multivariate normality
(Watanabe, 2022), although the distributions will remain the same for
correlation matrices in all elliptically contoured distributions
(see Anderson, 2003).
}
\examples{
# Covariance matrix
N <- 20
Lambda <- c(4, 2, 1, 1)
(Sigma <- GenCov(evalues = Lambda, evectors = "random"))
VE(S = Sigma)$VE
VE(S = Sigma)$VR
# Population values of V(Sigma) and Vrel(Sigma)

# From population covariance matrix
Exv.VES(Sigma, N - 1)
Var.VES(Sigma, N - 1)
Exv.VRS(Sigma, N - 1)
Var.VRS(Sigma, N - 1)
# Note the amount of bias from the population value obtained above

# From population eigenvalues
Exv.VES(Lambda = Lambda, n = N - 1)
Var.VES(Lambda = Lambda, n = N - 1)
Exv.VRS(Lambda = Lambda, n = N - 1)
Var.VRS(Lambda = Lambda, n = N - 1)
# Same, regardless of the random choice of eigenvectors

# Correlation matrix
(Rho <- GenCov(evalues = Lambda / sum(Lambda) * 4, evectors = "Givens"))
VE(S = Rho)$VR
# Population value of Vrel(Rho)

Exv.VRR(Rho, N - 1)
Var.VRR(Rho, N - 1)
# These results vary with the choice of eigenvalues
# If interested, repeat from the definition of Rho

# Different choices for asymptotic variance of Vrel(R)
# Variance from Pan-Frank method
Var.VRR(Rho, N - 1, method = "Pan-Frank") # Internally sets fun = "pfd"
Var.VRR(Rho, N - 1, fun = "pf")  # Slow for large p
Var.VRR(Rho, N - 1, fun = "pfv") # Requires too much RAM for large p
\dontrun{Var.VRR(Rho, n = N - 1, fun = "pfc")} # Requires eigvaldispRcpp
\dontrun{Var.VRR(Rho, n = N - 1, fun = "pfc", cppfun = "Cov_r2P")}
# The last two use \proglang{C++} functions provided by extension package eigvaldispRcpp
# fun = "pfc"' can be omitted in the last call.
# The above results are identical (up to rounding error)

# Variance from Konishi's theory
Var.VRR(Rho, N - 1, method = "Konishi") # Internally sets fun = "klv"
Var.VRR(Rho, N - 1, fun = "kl")
Var.VRR(Rho, N - 1, fun = "krv")
Var.VRR(Rho, N - 1, fun = "kr")
# These are identical, but the first one is fast
# On the other hand, these differ from that obtained with
# the Pan-Frank method above

}
\references{
Watanabe, J. (2022) Statistics of eigenvalue dispersion indices:
quantifying the magnitude of phenotypic integration. \emph{Evolution},
\strong{76}, 4--28. \doi{10.1111/evo.14382}.
}
\seealso{
\code{\link{VE}} for estimation

\code{\link{AVar.VRR_xx}} for internal functions of Var.VRR

\code{\link{Exv.rx}} for internal functions for moments of
correlation coefficients

\code{\link{Exv.VXXa}} for moments of \dQuote{bias-corrected} versions
}
