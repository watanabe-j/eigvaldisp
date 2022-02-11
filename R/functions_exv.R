##### Exv.VXX #####
#' Moments of eigenvalue dispersion indices
#'
#' Functions to calculate expectation/variance of (relative) eigenvalue variance
#' of sample covariance/correlation matrices for a given population
#' covariance/correlation matrix and degrees of freedom \eqn{n}.
#'
#' \code{Exv.VES()}, \code{Var.VES()}, and \code{Exv.VRR()} return exact
#' moments.
#' \code{Exv.VRS()} and \code{Var.VRS()} return approximations based
#' on the delta method, except under the null condition
#' (\eqn{\Sigma} proportional to the identity matrix)
#' where exact moments are returned.
#'
#' \code{Var.VRR()} returns the exact variance when \eqn{p = 2} or
#' under the null condition (\eqn{\Rho} is the identity matrix).
#' Otherwise, asymptotic variance is calculated with the \code{method} of
#' choice: either \code{"Pan-Frank"} (default) or \code{"Konishi"}.
#' They correspond to Pan & Frank's heuristic approximation and Konishi's
#' asymptotic theory, respectively (see Watanabe, 2022).
#'
#' In this case, calculations are handled by one of the internal functions
#' \code{AVar.VRR_xx()} (\code{xx} is a suffix to specify R implementation).
#' For completeness, it is possible to directly specify the function
#' to be used with the argument \code{fun}: for the Pan--Frank method,
#' \code{pfd} (default), \code{pfv}, and \code{pf}; and for the
#' Konishi method, \code{klv} (default), \code{kl}, \code{krv}, and \code{kr}.
#' Within each group, these function yield identical results but differ
#' in speed (the defaults are the fastest).
#' See \code{\link{AVar.VRR_xx}} for details of these functions.
#'
#' The Pan--Frank method takes a substantial amount of time to be executed
#' when p is large. Several C++ functions are provided in the extension package
#' [\code{eigvaldispRcpp}](https://github.com/watanabe-j/eigvaldispRcpp) to
#' speed-up the calculation. When this package is available, the default
#' \code{fun} for the Pan--Frank method is set to \code{"pfc"}.
#' The option for \code{C++} function is controlled by the argument
#' \code{cppfun} which in turn is passed to \code{AVar.VRR_pfc()}
#' (see \code{\link{AVar.VRR_xx}}). When this argument is provided, the argument
#' \code{fun} is ignored with a warning, unless one of the Konishi methods
#' is used (in which case \code{cppfun} is ignored with a warning).
#'
#' Since the eigenvalue variance of a correlation matrix \eqn{V(R)} is simply
#' \eqn{(p - 1)} times the relative eigenvalue variance \eqn{Vrel(R)} of
#' the same matrix, their distributions are identical up to this scaling.
#' Hence, \code{Exv.VER()} and \code{Var.VER()} calls \code{Exv.VRR()} and
#' \code{Var.VRR()} (respectively), whose outputs are scaled and returned.
#' These functions are provided for completeness, although there will be little
#' practical demand for these functions (and \eqn{V(R)} itself).
#'
#' As detailed in Watanabe (2022), the distribution of \eqn{Vrel(R)} cannot be
#' uniquely specified by eigenvalues alone. Hence, a full correlation
#' matrix \code{Rho} should preferably be provided. Otherwise, a correlation
#' matrix is constructed from the eigenvalues \code{Lambda} provided using
#' the function \code{GenCov()} with randomly picked eigenvectors.
#'
#' On the other hand, the choice of eigenvectors does not matter for
#' covariance matrices, thus either the full covariance matrix \code{Sigma} or
#' vector of eigenvalues \code{Lambda} can be provided (although the latter
#' is slightly faster if available).
#'
#' When \code{Rho} is provided, some simple checks are done: the matrix is
#' scaled to have diagonals of 1; and if any of these are unequal,
#' an error is returned.
#'
#' These moments are derived under the assumption of multivariate normality
#' (Watanabe, 2022), although the distributions will remain the same for
#' correlation matrices in all elliptically contoured distributions
#' (see Anderson, 2003).
#'
# #' For covariance matrices, the divisor of \eqn{n}
# #' (which gives the ordinary unbiased estimator) is assumed by default.
# #'
# #' \code{Exv.VRR()} calls \code{Exv.r2()}, which in turn calls \code{hgf()}.
# #' (see \code{\link{Exv.rx}}).
# #'
#' @name Exv.VXX
#'
#' @param Sigma
#'   Population covariance matrix; assumed to be validly constructed.
#' @param Rho
#'   Population correlation matrix; assumed to be validly constructed
#'   (although simple checks are done).
#' @param n
#'   Degrees of freedom (not sample sizes); numeric of length 1 or more.
#' @param Lambda
#'   Numeric vector of population eigenvalues.
#' @param divisor
#'   Either \code{"UB"} (default) or \code{"ML"},
#'   to decide the default value of \code{m}.
#' @param m
#'   Divisor for the sample covariance matrix (\eqn{n*} in Watanabe (2022)).
#'   By default equals \eqn{n}.
#' @param drop_0
#'   Logical, when \code{TRUE}, eigenvalues smaller than \code{tol} are dropped.
#' @param tol
#'   For covariance-related functions, this is the tolerance/threshold
#'   to be used with drop_0. For correlation-related functions,
#'   this is passed to \code{Exv.r2()} along with other arguments.
#' @param tol.hg,maxiter.hg
#'   Passed to \code{Exv.r2()}; see description of that function.
#' @param method
#'   For \code{Var.VRR()} (and \code{Var.VER()}), determines the method
#'   to obtain approximate variance in non-null conditions.
#'   Either \code{"Pan-Frank"} (default) or \code{"Konishi"}. See Details.
#' @param fun
#'   For \code{Var.VRR()} (and \code{Var.VER()}), determines the function
#'   to be used to evaluate approximate variance. See Details.
#'   Options allowed are: \code{"pfd"}, \code{"pfv"}, \code{"pfc"}, \code{"pf"},
#'   \code{"klv"}, \code{"kl"}, \code{"krv"}, and \code{"kr"}.
#' @param ...
#'   In \code{Var.VRR()}, additional arguments are passed to an internal
#'   function which it in turn calls. Otherwise ignored.
#'
#' @return
#' A numeric vector of the desired moment, corresponding to \code{n}.
#'
#' @references
#' Watanabe, J. (2022). Statistics of eigenvalue dispersion indices:
#'  quantifying the magnitude of phenotypic integration. *Evolution*,
#'  **76**, 4--28. doi:[10.1111/evo.14382](https://doi.org/10.1111/evo.14382).
#'
#' @seealso
#' \code{\link{VE}} for estimation
#'
#' \code{\link{AVar.VRR_xx}} for internal functions of Var.VRR
#'
#' \code{\link{Exv.rx}} for internal functions for moments of
#'   correlation coefficients
#'
#' \code{\link{Exv.VXXa}} for moments of ``bias-corrected'' versions
#'
#' @examples
#' # Covariance matrix
#' N <- 20
#' Lambda <- c(4, 2, 1, 1)
#' (Sigma <- GenCov(evalues = Lambda, evectors = "random"))
#' VE(S = Sigma)$VE
#' VE(S = Sigma)$VR
#' # Population values of V(Sigma) and Vrel(Sigma)
#'
#' # From population covariance matrix
#' Exv.VES(Sigma, N - 1)
#' Var.VES(Sigma, N - 1)
#' Exv.VRS(Sigma, N - 1)
#' Var.VRS(Sigma, N - 1)
#' # Note the amount of bias from the population value obtained above
#'
#' # From population eigenvalues
#' Exv.VES(Lambda = Lambda, n = N - 1)
#' Var.VES(Lambda = Lambda, n = N - 1)
#' Exv.VRS(Lambda = Lambda, n = N - 1)
#' Var.VRS(Lambda = Lambda, n = N - 1)
#' # Same, regardless of the random choice of eigenvectors
#'
#' # Correlation matrix
#' (Rho <- GenCov(evalues = Lambda / sum(Lambda) * 4, evectors = "Givens"))
#' VE(S = Rho)$VR
#' # Population value of Vrel(Rho)
#'
#' Exv.VRR(Rho, N - 1)
#' Var.VRR(Rho, N - 1)
#' # These results vary with the choice of eigenvalues
#' # If interested, repeat from the definition of Rho
#'
#' # Different choices for asymptotic variance of Vrel(R)
#' # Variance from Pan-Frank method
#' Var.VRR(Rho, N - 1, method = "Pan-Frank") # Internally sets fun = "pfd"
#' Var.VRR(Rho, N - 1, fun = "pf")  # Slow for large p
#' Var.VRR(Rho, N - 1, fun = "pfv") # Requires too much RAM for large p
#' \dontrun{Var.VRR(Rho, n = N - 1, fun = "pfc")} # Requires eigvaldispRcpp
#' \dontrun{Var.VRR(Rho, n = N - 1, fun = "pfc", cppfun = "Cov_r2P")}
#' # The last two use C++ functions provided by extension package eigvaldispRcpp
#' # fun = "pfc"' can be omitted in the last call.
#' # The above results are identical (up to rounding error)
#'
#' # Variance from Konishi's theory
#' Var.VRR(Rho, N - 1, method = "Konishi") # Internally sets fun = "klv"
#' Var.VRR(Rho, N - 1, fun = "kl")
#' Var.VRR(Rho, N - 1, fun = "krv")
#' Var.VRR(Rho, N - 1, fun = "kr")
#' # These are identical, but the first one is fast
#' # On the other hand, these differ from that obtained with
#' # the Pan-Frank method above
#'
NULL

##### Exv.VES #####
#' Expectation of eigenvalue variance of covariance matrix
#'
#' \code{Exv.VES()}: expectation of eigenvalue variance of covariance matrix
#' \eqn{E[V(S)]}.
#'
#' @rdname Exv.VXX
#'
#' @export
#'
Exv.VES <- function(Sigma, n = 100, Lambda, divisor = c("UB", "ML"),
                    m = switch(divisor, UB = n, ML = n + 1), drop_0 = FALSE,
                    tol = .Machine$double.eps * 100, ...) {
    divisor <- match.arg(divisor)
    if(missing(Lambda)) {
        Lambda <- svd(Sigma, nu = 0, nv = 0)$d
    }
    if(drop_0) {
        Lambda <- Lambda[Lambda > tol]
    } else {
        Lambda[Lambda < tol] <- 0
    }
    p <- length(Lambda)
    t1 <- sum(Lambda)
    t2 <- sum(Lambda ^ 2)
    n * ((p - n) * t1 ^ 2 + (p * n + p - 2) * t2) / (p ^ 2 * m ^ 2)
}

##### Exv.VRS #####
#' Expectation of relative eigenvalue variance of covariance matrix
#'
#' \code{Exv.VRS()}: expectation of relative eigenvalue variance of
#' covariance matrix \eqn{E[Vrel(S)]}.
#'
#' @rdname Exv.VXX
#'
#' @export
#'
Exv.VRS <- function(Sigma, n = 100, Lambda, drop_0 = FALSE,
                    tol = .Machine$double.eps * 100, ...) {
    if(missing(Lambda)) {
        Lambda <- svd(Sigma, nu = 0, nv = 0)$d
    }
    if(drop_0) {
        Lambda <- Lambda[Lambda > tol]
    } else {
        Lambda[Lambda < tol] <- 0
    }
    p <- length(Lambda)
    t1 <- sum(Lambda)
    t2 <- sum(Lambda ^ 2)
    t3 <- sum(Lambda ^ 3)
    t4 <- sum(Lambda ^ 4)
    EF <- (t1 ^ 2 + (n + 1) * t2) / (n * t1 ^ 2 + 2 * t2) -
          8 * ((n + 2) * (n - 1) * (3 * t4 * t1 ^ 2 - 2 * t3 * t2 * t1 - t2 ^ 3
                                    + n * t3 * t1 ^ 3 - n * t2 ^ 2 * t1 ^ 2)) /
          (n * ((n * t1 ^ 2 + 2 * t2)) ^ 3)
    (p * EF - 1) / (p - 1)
}

##### Exv.VER #####
#' Expectation of eigenvalue variance of correlation matrix
#'
#' \code{Exv.VER()}: expectation of eigenvalue variance of correlation
#' matrix \eqn{E[V(R)]}.
#'
#' @rdname Exv.VXX
#'
#' @export
#'
Exv.VER <- function(Rho, n = 100, Lambda, tol = .Machine$double.eps * 100,
                    tol.hg = 0, maxiter.hg = 2000, ...) {
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    exv.VRR <- Exv.VRR(Rho = Rho, n = n, tol = tol,
                       tol.hg = tol.hg, maxiter.hg = maxiter.hg, ...)
    exv.VRR * (p - 1)
}

##### Exv.VRR #####
#' Expectation of relative eigenvalue variance of correlation matrix
#'
#' \code{Exv.VRR()}: expectation of relative eigenvalue variance of
#' correlation matrix \eqn{E[Vrel(R)]}.
#'
#' @rdname Exv.VXX
#'
#' @export
#'
Exv.VRR <- function(Rho, n = 100, Lambda, tol = .Machine$double.eps * 100,
                    tol.hg = 0, maxiter.hg = 2000, ...) {
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    R2 <- Rho[lower.tri(Rho)] ^ 2
    exv_r2 <- matrix(sapply(n, Exv.r2, x = R2, do.square = FALSE, tol = tol,
                            tol.hg = tol.hg, maxiter.hg = maxiter.hg),
                     ncol = length(n))
    2 * colSums(exv_r2) / (p * (p - 1))
}

##### Var.VES #####
#' Variance of eigenvalue variance of covariance matrix
#'
#' \code{Var.VES()}: variance of eigenvalue variance of
#' covariance matrix \eqn{Var[V(S)]}.
#'
#' @rdname Exv.VXX
#'
#' @export
#'
Var.VES <- function(Sigma, n = 100, Lambda, divisor = c("UB", "ML"),
                    m = switch(divisor, UB = n, ML = n + 1), drop_0 = FALSE,
                    tol = .Machine$double.eps * 100, ...) {
    divisor <- match.arg(divisor)
    if(missing(Lambda)) {
        Lambda <- svd(Sigma, nu = 0, nv = 0)$d
    }
    if(drop_0) {
        Lambda <- Lambda[Lambda > tol]
    } else {
        Lambda[Lambda < tol] <- 0
    }
    p <- length(Lambda)
    t1 <- sum(Lambda)
    t2 <- sum(Lambda ^ 2)
    t3 <- sum(Lambda ^ 3)
    t4 <- sum(Lambda ^ 4)
    4 * n * ((2 * p^2 * n^2 + 5 * p^2 * n + 5 * p^2 - 12 * p * n - 12 * p + 12)
             * t4 + 4 * (p - n) * (p * n + p - 2) * t3 * t1 +
             (p^2 * n + p^2 - 4 * p + 2 * n) * t2^2 +
             2 * (p - n)^2 * t2 * t1^2) / (p ^ 4 * m ^ 4)
}

##### Var.VRS #####
#' Variance of relative eigenvalue variance of covariance matrix
#'
#' \code{Var.VRS()}: variance of relative eigenvalue variance of
#' covariance matrix \eqn{Var[Vrel(S)]}.
#'
#' @rdname Exv.VXX
#'
#' @export
#'
Var.VRS <- function(Sigma, n = 100, Lambda, drop_0 = FALSE,
                    tol = .Machine$double.eps * 100, ...) {
    if(missing(Lambda)) {
        Lambda <- svd(Sigma, nu = 0, nv = 0)$d
    }
    if(drop_0) {
        Lambda <- Lambda[Lambda > tol]
    } else {
        Lambda[Lambda < tol] <- 0
    }
    p <- length(Lambda)
    t1 <- sum(Lambda)
    t2 <- sum(Lambda ^ 2)
    t3 <- sum(Lambda ^ 3)
    t4 <- sum(Lambda ^ 4)
    if(isTRUE(all.equal(Lambda / Lambda[1], rep.int(1, p)))) {
        4 * p^2 * (p + 2) * (n - 1) * (n + 2) /
             (p - 1) / (p * n + 2)^2 / (p * n + 4) / (p * n + 6)
    } else {
        4 * p^2 / (p - 1)^2 * (n - 1) * (n + 2) *
            (- 4 * t4 * t2^2 - 4 * n * t4 * t2 * t1^2
             + (2 * n^2 + 3 * n - 6) * t4 * t1^4
             - 4 * (n - 1) * (n + 2) * t3 * t2 * t1^3
             + 2 * (n + 1) * t2^4 + 2 * n * (n + 1) * t2^3 * t1^2
             + n * t2^2 * t1^4) / n / (2 * t2 + n * t1^2)^4
    }
}

##### Var.VER #####
#' Variance of eigenvalue variance of correlation matrix
#'
#' \code{Var.VER()}: variance of eigenvalue variance of
#' correlation matrix \eqn{Var[V(R)]}.
#'
#' @rdname Exv.VXX
#'
#' @export
#'
Var.VER <- function(Rho, n = 100, Lambda, ...) {
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    var.VRR <- Var.VRR(Rho = Rho, n = n, ...)
    var.VRR * (p - 1) ^ 2
}

##### Var.VRR #####
#' Variance of relative eigenvalue variance of correlation matrix
#'
#' \code{Var.VRR()}: variance of relative eigenvalue variance of
#' correlation matrix \eqn{Var[Vrel(R)]}.
#'
#' @rdname Exv.VXX
#'
#' @export
#'
Var.VRR <- function(Rho, n = 100, method = c("Pan-Frank", "Konishi"), Lambda,
                    fun = c("pfd", "pfv", "pfc", "pf",
                            "klv", "kl", "krv", "kr"), ...) {
    fun_missing <- missing(fun)
    if(fun_missing) {
        method <- match.arg(method)
        if(method == "Pan-Frank") {
            if(requireNamespace("eigvaldispRcpp", quietly = TRUE)) {
                fun <- "pfc"
            } else {
                fun <- "pfd"
            }
        } else {
            fun <- "klv"
        }
    } else {
        fun <- match.arg(fun)
    }
    if(any(grepl("^cp", names(list(...))))) {
        if(grepl("pf", fun)) {
            fun <- "pfc"
            if(!fun_missing) {
                warning("The argument 'fun' was ignored as another argument ",
                        "that seems to match 'cppfun' was provided")
            }
        } else {
            warning("The argument 'cppfun' was ignored; it is used ",
                    "only when method = 'Pan-Frank'")
        }
    }
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    if(isTRUE(all.equal(Rho, diag(p)))) {
        4 * (n - 1) / (p * (p - 1) * n^2 * (n + 2))
    } else if(p == 2) {
        R2 <- Rho[lower.tri(Rho)] ^ 2
        var_r2 <- matrix(sapply(n, Var.r2, x = R2, do.square = FALSE),
                         ncol = length(n))
        4 * colSums(var_r2) / (p * (p - 1)) ^ 2
    } else {
        Fun <- switch(fun,
                      klv = AVar.VRR_klv, kl = AVar.VRR_kl, kr = AVar.VRR_kr,
                      krv = AVar.VRR_krv, pf = AVar.VRR_pf, pfv = AVar.VRR_pfv,
                      pfd = AVar.VRR_pfd, pfc = AVar.VRR_pfc)
        Fun(Rho = Rho, n = n, ...)
    }
}

##### AVar.VRR_xx #####
#' Approximate variance of relative eigenvalue variance of correlation matrix
#'
#' Functions to obtain approximate variance of relative eigenvalue variance
#' of correlation matrix \eqn{Var[Vrel(R)]}. There are several versions
#' for each of two different expressions: \code{pf*} and \code{k*} families.
#'
#' Watanabe (2022) presented two approaches to evaluate approximate variance
#' of the relative eigenvalue variance of a correlation matrix \eqn{Vrel(R)}.
#' One is Pan & Frank's (2004) heuristic approximation (eqs. 28 and 36--38 in
#' Watanabe 2022). The other is based on Konishi's (1979) asymptotic
#' theory (eq. 39 in Watanabe 2022). Simulations showed that the former tends
#' to be more accurate, but the latter is much faster. This is mainly because
#' the Pan--Frank approach involves evaluation of covariances in \eqn{~p^4 / 4}
#' pairs of (squared) correlation coefficients.
#' (That said, the speed will not be a practical concern unless \eqn{p}
#' exceeds a few hundreds.)
#'
#' The Pan--Frank approach is at present implemented in several functions
#' which yield (almost) identical results:
#' \describe{
#'   \item{\code{AVar.VRR_pf()}}{Prototype version. Simplest implementation.}
#'   \item{\code{AVar.VRR_pfv()}}{Vectorized version. Much faster, but
#'     requires a large RAM space as \code{p} grows.}
#'   \item{\code{AVar.VRR_pfd()}}{Improvement over \code{AVar.VRR_pfv()}.
#'     Faster and more RAM-efficient. This is the default to be called in
#'     \code{Var.VRR(..., method = "Pan-Frank")}, unless the extension package
#'     \code{eigvaldispRcpp} is installed.}
#'   \item{\code{AVar.VRR_pfc()}}{Fast version using \code{Rcpp}.
#'     Requires the extension package \code{eigvaldispRcpp}; this is
#'     the default when this package is installed (and detected).}
#' }
#' \code{AVar.VRR_pfc()} implements the same algorithm as the others,
#' but makes use of \code{C++} API via the package \code{Rcpp} for evaluation of
#' the sum of covariance across pairs of squared correlation coefficients.
#' This version works much faster than vectorized \code{R} implementations.
#' Note that the output can slightly differ from those of pure \code{R}
#' implementations (by the order of ~1e-9).
#'
#' The Konishi approach is implemented in several functions:
#' \describe{
#'   \item{\code{AVar.VRR_kl()}}{From Konishi (1979: corollary 2.2):
#'     \eqn{Vrel(R)} as function of eigenvalues. Prototype version.}
#'   \item{\code{AVar.VRR_klv()}}{Vectorized version of \code{AVar.VRR_kl()}.
#'     This is the default when \code{Var.VRR(..., method = "Konishi")}.}
#'   \item{\code{AVar.VRR_kr()}}{From Konishi (1979: theorem 6.2):
#'     \eqn{Vrel(R)} as function of correlation coefficients.}
#'   \item{\code{AVar.VRR_krv()}}{Vectorized version of \code{AVar.VRR_kr()};
#'     slightly faster for moderate \eqn{p}, but not particularly fast
#'     for large \eqn{p} as the number of elements to be summed becomes large.}
#' }
#' Empirically, these all yield the same result, but
#' \code{AVar.VRR_klv()} is by far the fastest.
#'
#' The \code{AVar.VRR_pfx()} family functions by default return exact variance
#' when \eqn{p = 2},
#' If asymptotic result is desired, use \code{mode.var2 = "asymptotic"}.
#'
#' Options for \code{mode} in \code{AVar.VRR_pf()} and \code{AVar.VRR_pfd()}:
#' \describe{
#'   \item{\code{"nested.for"}}{Only for \code{AVar.VRR_pf()}. Uses nested
#'     for loops, which is straifhgforward and RAM efficient but slow.}
#'   \item{\code{"for.ind"} (default)/\code{"lapply"}}{Run the iteration along
#'     an index vector to shorten computational time,
#'     with \code{for} loop and \code{lapply()}, respectively.}
#'   \item{\code{"mclapply"}/\code{"parLapply"}}{Only for \code{AVar.VRR_pfd()}.
#'     Parallelize the same iteration by
#'     forking and socketing, respectively, with the named functions in the
#'     package \code{parallel}. Note that the former doesn't work in the
#'     Windows environment. See \code{vignette("parallel")} for details.}
#' }
#'
#' \code{AVar.VRR_pfv()} internally generates vectors and matrices
#' whose lengths are about \eqn{p^4 / 8} and \eqn{p^4 / 4}. These take about
#' \eqn{2*p^4} bytes of RAM; this can be prohibitively large for large \eqn{p}.
#'
#' \code{AVar.VRR_pfd()} divides the index vector \code{b} (used in
#' \code{AVar.VRR_pfv()}) into a list \code{bd} using the internal function
#' \code{eigvaldisp:::divInd()}. The calculations are then done on each element
#' of this list to save RAM space.
#' This process takes some time when \eqn{p} is large (~10 sec
#' for \eqn{p = 1024}).
#' Alternatively, this list can be provided as the argument \code{bd}
#' (which should exactly match the one to be generated; use
#' \code{eigvaldisp:::divInd()}).
#' The argument \code{max.size} controls the maximum size of resulting vectors;
#' at least \code{max.size * (2 * length(n) + 6) * 8} bytes of RAM is required
#' for storing temporary results (and more during computation);
#' e.g., ~2e7 seems good for 16 GB RAM, ~4e8 for 256 GB.
#' However, performance does not seem to improve past 1e6--1e7 presumably
#' because memory allocation takes substantial time for large objects.
#' The iteration can be parallelized with \code{mode = "mclapply"} or
#' \code{"parLapply"}, but be careful about RAM limitations.
#'
#' \code{AVar.VRR_pfc()} provides a faster implementation with one of the
#' \code{C++} functions defined in the extension package \code{eigvaldispRcpp}
#' (which is required to run this function).
#' The \code{C++} function is specified by the argument \code{cppfun}:
#' \describe{
#'   \item{\code{"Cov_r2C"} (default)}{Serial evaluation with base
#'     \code{Rcpp} functionalities.}
#'   \item{\code{"Cov_r2A"} or \code{"Armadillo"}}{Using \code{RcppArmadillo}.
#'     Parallelized with OpenMP when the environment allows.}
#'   \item{\code{"Cov_r2E"} or \code{"Eigen"}}{Using \code{RcppEigen}.
#'     Parallelized with OpenMP when the environment allows.}
#'   \item{\code{"Cov_r2P"} or \code{"Parallel"}}{Using \code{RcppParallel}.
#'     Parallelized with IntelTBB when the environment allows.}
#' }
#' The default option would be sufficiently fast for up to p = 100 or so.
#' The latter three options aim at speeding-up the calculation via
#' parallelization with other \code{Rcpp}-related packages. Although these
#' would have similar performance in most environments, \code{"Cov_r2E"} seems
#' the fastest in the development environment, closely followed by
#' \code{"Cov_r2A"}.
#'
#' @name AVar.VRR_xx
#'
#' @inheritParams Exv.VXX
#'
#' @param exv1.mode
#'   Whether \code{"exact"} or \code{"asymptotic"} expression is used for
#'   \eqn{E(r)}.
#' @param var2.mode
#'   Whether \code{"exact"} or \code{"asymptotic"} expression is used for
#'   \eqn{Var(r^2)}.
#' @param var1.mode
#'   Whether \code{"exact"} or \code{"asymptotic"} expression is used for
#'   \eqn{Cov(rij, rkl)}. (At present, only \code{"asymptotic"} is allowed.)
#' @param order.exv1,order.var2
#'   Used to specify the order of evaluation for asymptotic expressions of
#'   \eqn{E(r)} and \eqn{Var(r^2)} when \code{exv1.mode} and \code{var2.mode} is
#'   \code{"asymptotic"}; see \code{\link{Exv.rx}}.
#' @param mode
#'   In \code{AVar.VRR_pf()} and \code{AVar.VRR_pfd()},
#'   specifies the mode of iterations (see Details).
#' @param mc.cores
#'   Number of cores to be used (numeric/integer). When \code{"auto"} (default),
#'   set to \code{min(c(ceiling(p / 2), max.cores))}, which usually works well.
#'   (Used only when \code{mode = "mclapply"}, or \code{"parLapply"})
#' @param max.cores
#'   Maximum number of cores to be used.
#'   (Used only when \code{mode = "mclapply"}, or \code{"parLapply"})
#' @param do.mcaffinity
#'   Whether to run \code{parallel::mcaffinity()}, which seems required in
#'   some Linux environments to assign threads to multiple cores.
#'   (Used only when \code{mode = "mclapply"}, or \code{"parLapply"})
#' @param affinity_mc
#'   Argument of \code{parallel::mcaffinity()} to specify assignment of threads.
#'   (Used only when \code{mode = "mclapply"}, or \code{"parLapply"})
#' @param cl
#'   A cluster object (made by \code{parallel::makeCluster()}); when already
#'   created, one can be specified with this argument. Otherwise, one is created
#'   within function call, which is turned off on exit.
#'   (Used only when \code{mode = "parLapply"})
#' @param max.size
#'   Maximum size of vectors created internally (see Details).
#' @param bd
#'   List of indices used for iteration (see Details).
#' @param verbose
#'   When \code{"yes"} or \code{"inline"}, pogress of iteration is printed
#'   on console. \code{"no"} (default) turns off the printing.
#'   To be used in \code{AVar.VRR_pfd()} for large \eqn{p} (hundreds or more).
#' @param cppfun
#'   Option to specify the C++ function to be used (see Details).
#' @param nthreads
#'   Integer to specify the number of threads used in OpenMP parallelization
#'   in \code{"Cov_r2A"} and \code{"Cov_r2E"}. By default (0), the number of
#'   threads is automatically set to one-half of that of (logical) processors
#'   detected (by the \code{C++} function \code{omp_get_num_procs()}).
#'   Setting this beyond the number of physical processors can result in
#'   poorer performance, depending on the environment.
#' @param ...
#'   In the \code{pf} family functions, passed to \code{Exv.r1()} and
#'   \code{Var.r2()} (when the corresponding modes are \code{"exact"}).
#'   Otherwise ignored.
#'
#' @return
#' A numeric vector of \eqn{Var[Vrel(R)]}, corresponding to \code{n}.
#'
#' @references
#' Konishi, S. (1979). Asymptotic expansions for the distributions of statistics
#'  based on the sample correlation matrix in principal componenet analysis.
#'  *Hiroshima Mathematical Journal* **9**, 647--700.
#'  doi:[10.32917/hmj/1206134750](https://doi.org/10.32917/hmj/1206134750).
#'
#' Pan, W. & Frank, K. A. (2004). An approximation to the distribution of the
#'  product of two dependent correlation coefficients. *Journal of Statistical
#'  Computation and Simulation* **74**, 419--443.
#'  doi:[10.1080/00949650310001596822](https://doi.org/10.1080/00949650310001596822).
#'
#' Watanabe, J. (2022). Statistics of eigenvalue dispersion indices:
#'  quantifying the magnitude of phenotypic integration. *Evolution*,
#'  **76**, 4--28. doi:[10.1111/evo.14382](https://doi.org/10.1111/evo.14382).
#'
#' @seealso \code{\link{Exv.VXX}} for main moment functions
#'
#' @examples
#' # See also examples of Exv.VXX
#' # Correlation matrix
#' N <- 20
#' Lambda <- c(4, 2, 1, 1)
#' (Rho <- GenCov(evalues = Lambda / sum(Lambda) * 4, evectors = "Givens"))
#'
#' # Different choices for asymptotic variance of Vrel(R)
#' # Variance from Pan-Frank method
#' Var.VRR(Rho, n = N - 1) # By default, method = "Pan-Frank" and AVar.VRR_pfd() is called
#' eigvaldisp:::AVar.VRR_pfd(Rho, n = N - 1) # Same as above
#' Var.VRR(Rho, n = N - 1, fun = "pf")  # Calls AVar.VRR_pf(), which is slow for large p
#' Var.VRR(Rho, n = N - 1, fun = "pfv") # Calls AVar.VRR_pfv(), which requires much RAM for large p
#' # Various implementations with Rcpp (require eigvaldispRcpp):
#' \dontrun{Var.VRR(Rho, n = N - 1, fun = "pfc")} # By default, cppfun = "Cov_r2C" is used
#' \dontrun{Var.VRR(Rho, n = N - 1, cppfun = "Cov_r2A")}
#' \dontrun{Var.VRR(Rho, n = N - 1, cppfun = "Cov_r2E")}
#' \dontrun{Var.VRR(Rho, n = N - 1, cppfun = "Cov_r2P")}
#' # When the argument cppfun is provided, fun need not be specified
#' # The above results are identical
#'
#' # Variance from Konishi's theory
#' Var.VRR(Rho, n = N - 1, method = "Konishi")  # By default for this method, AVar.VRR_klv() is called
#' eigvaldisp:::AVar.VRR_klv(Rho, n = N - 1)    # Same as above
#' Var.VRR(Rho, n = N - 1, fun = "kl")
#' Var.VRR(Rho, n = N - 1, fun = "kr")
#' Var.VRR(Rho, n = N - 1, fun = "krv")
#' # The results are identical, but the last three are slower
#' # On the other hand, these differ from that obtained with the Pan-Frank method
#'
#' # Example with p = 2
#' Rho2 <- GenCov(evalues = c(1.5, 0.5), evectors = "Givens")
#' Var.VRR(Rho2, n = N - 1)  # When p = 2, this does not call AVar.VRR_pfd()
#' eigvaldisp:::AVar.VRR_pfd(Rho2, n = N - 1)
#' # By default, the above returns the same, exact result
#'
#' eigvaldisp:::AVar.VRR_pfd(Rho2, n = N - 1, var2.mode = "asymptotic")
#' eigvaldisp:::AVar.VRR_klv(Rho2, n = N - 1)
#' # These return different asymptotic expressions
#'
NULL

##### AVar.VRR_pf #####
#' Approximate variance of relative eigenvalue variance of correlation matrix
#'
#' \code{AVar.VRR_pf()}: asymptotic and approximate variance of \eqn{Vrel(R)}
#' based on Pan & Frank's (2004) approach.  Prototype version.
#'
#' @rdname AVar.VRR_xx
#'
AVar.VRR_pf <- function(Rho, n = 100, Lambda, exv1.mode = c("exact", "asymptotic"),
                        var1.mode = "asymptotic",
                        var2.mode = c("exact", "asymptotic"),
                        order.exv1 = 2, order.var2 = 2,
                        mode = c("for.ind", "nested.for", "lapply"), ...) {
                        # mode = c("for.ind", "nested.for", "lapply",
                        #          "mclapply", "parLapply"),
                        # mc.cores = "auto", max.cores = parallel::detectCores(),
                        # do.mcaffinity = TRUE,
                        # affinity_mc = seq_len(max.cores), cl = NULL, ...) {
    exv1.mode <- match.arg(exv1.mode)
    var1.mode <- match.arg(var1.mode)
    var2.mode <- match.arg(var2.mode)
    mode <- match.arg(mode)
    # if(mode == "mclapply" || mode == "parLapply") {
    #     if(!requireNamespace("parallel", quietly = TRUE)) {
    #         stop("Package 'parallel' is required for ",
    #              "mode = 'mclapply' or 'parLapply'")
    #     }
    # }
    getInds <- function(p){
        a <- seq_len(p)
        d <- digit(p)
        A <- outer(as.integer(10 ^ d) * a, a, "+")
        b <- t(A)[lower.tri(A)]
        B <- outer((100 ^ d) * b, b, "+")
        t(B)[lower.tri(B)]
    }
    parseInds <- function(x, d) {
        i <- (x %/% 1000 ^ d) %% 10 ^ d
        j <- (x %/% 100 ^ d) %% 10 ^ d
        k <- (x %/% 10 ^ d) %% 10 ^ d
        l <- x %% 10 ^ d
        c(i, j, k, l)
    }
    Cov_r2s <- function(I = NULL, pI) {
        if(!missing(I)) pI <- parseInds(I, d)
        i <- pI[1]
        j <- pI[2]
        k <- pI[3]
        l <- pI[4]
        Eij <- exv_r1[(j - 1) * p + i - (j - 1) * (2 * p - j + 2) / 2, ]
        Ekl <- exv_r1[(l - 1) * p + k - (l - 1) * (2 * p - l + 2) / 2, ]
        Cijkl <- v1fun(n = n, Rho = Rho, i = i, j = j, k = k, l = l)
        (4 * Eij * Ekl + 2 * Cijkl) * Cijkl
    }
    e1fun <- switch(exv1.mode,
                    exact = function(n, x) Exv.r1(n, x, ...),
                    asymptotic = function(n, x) AExv.r1(n, x,
                                                        order. = order.exv1))
    v1fun <- switch(var1.mode, ACov.r1)
    v2fun <- switch(var2.mode,
                    exact = function(n, x) Var.r2(n, x, do.square = TRUE, ...),
                    asymptotic = function(n, x) AVar.r2(n, x,
                                                        order. = order.var2))
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }

    p <- ncol(Rho)
    l_n <- length(n)
    R_u <- Rho[upper.tri(Rho)]
    var_r2 <- matrix(sapply(n, v2fun, x = R_u), ncol = l_n)
    exv_r1 <- matrix(sapply(n, e1fun, x = R_u), ncol = l_n)
    d <- digit(p)
    if(mode == "nested.for") {
        cov_r2 <- numeric(l_n)
        for(i in 1:(p - 2)) {
            for(j in (i + 1):p) {
                for(k in i:(p - 1)) {
                    for(l in (k + 1):p) {
                        if(i >= k && j >= l) next
                        Eij <- exv_r1[(j - 1) * p +
                                      i - (j - 1) * (2 * p - j + 2) / 2, ]
                        Ekl <- exv_r1[(l - 1) * p +
                                      k - (l - 1) * (2 * p - l + 2) / 2, ]
                        Cijkl <- sapply(n, v1fun, Rho = Rho,
                                        i = i, j = j, k = k, l = l)
                        cov_r2 <- cov_r2 + (4 * Eij * Ekl + 2 * Cijkl) * Cijkl
                    }
                }
            }
        }
    } else {
        Inds <- getInds(p)
        if(mode == "for.ind") {
            cov_r2 <- numeric(l_n)
            for(I in seq_along(Inds)) {
                cov_r2 <- cov_r2 + Cov_r2s(Inds[I])
            }
        } else {
            # if(mode == "lapply") {
            cov_r2 <- lapply(Inds, Cov_r2s)
            # } else {
            #     # require(parallel)
            #     if(mc.cores == "auto") {
            #         mc.cores <- pmin(ceiling(p / 2), max.cores)
            #     }
            #     if(mode == "mclapply") {
            #         if(do.mcaffinity) {
            #             try(invisible(parallel::mcaffinity(affinity_mc)),
            #                 silent = TRUE)
            #         }
            #         cov_r2 <- parallel::mclapply(Inds, Cov_r2s,
            #                                      mc.cores = mc.cores)
            #     } else { # if(mode == "parLapply")
            #         if(is.null(cl)) {
            #             cl <- parallel::makeCluster(mc.cores)
            #             on.exit(parallel::stopCluster(cl))
            #             if(do.mcaffinity) {
            #                 try(parallel::clusterApply(cl, as.list(affinity_mc),
            #                     parallel::mcaffinity), silent = TRUE)
            #             }
            #         }
            #         cov_r2 <- parallel::parLapply(cl = cl, Inds, Cov_r2s)
            #     }
            # }
            cov_r2 <- matrix(unlist(cov_r2), ncol = length(Inds))
            cov_r2 <- rowSums(cov_r2)
        }
    }
    ans <- colSums(var_r2) + 2 * cov_r2
    4 * ans / (p * (p - 1))^2
}

##### AVar.VRR_pfv #####
#' Approximate variance of relative eigenvalue variance of correlation matrix,
#' vectorized
#'
#' \code{AVar.VRR_pfv()}: vectorized version of \code{AVar.VRR_pf()}.
#' Much faster, but requires a large RAM space as \code{p} grows.
#'
#' @rdname AVar.VRR_xx
#'
AVar.VRR_pfv <- function(Rho, n = 100, Lambda, exv1.mode = c("exact", "asymptotic"),
                       var1.mode = "asymptotic",
                       var2.mode = c("exact", "asymptotic"),
                       order.exv1 = 2, order.var2 = 2, ...) {
    exv1.mode <- match.arg(exv1.mode)
    var1.mode <- match.arg(var1.mode)
    var2.mode <- match.arg(var2.mode)
    rep_d <- function(x, from = 1, to = length(x)) {
        x <- x[seq.int(from, length(x))]
        sx <- seq_along(x)[seq_len(to - from + 1)]
        unlist(lapply(sx, function(i) x[-seq_len(i)]))
    }
    e1fun <- switch(exv1.mode,
                    exact = function(n, x) Exv.r1(n, x, ...),
                    asymptotic = function(n, x) AExv.r1(n, x,
                                                        order. = order.exv1))
    c1fun <- switch(var1.mode, function(n, Rij, Rkl, Rik, Rjl, Ril, Rjk) {
        A <- (Rij * Rkl * (Rik^2 + Ril^2 + Rjk^2 + Rjl^2) / 2 + Rik * Rjl +
             Ril * Rjk - (Rij * Rik * Ril + Rij * Rjk * Rjl + Rik * Rjk * Rkl +
                        Ril * Rjl * Rkl))
        outer(A, n, "/")
    })
    v2fun <- switch(var2.mode,
                    exact = function(n, x) Var.r2(n, x, do.square = TRUE, ...),
                    asymptotic = function(n, x) AVar.r2(n, x,
                                                        order. = order.var2))
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    l_n <- length(n)
    R_u <- Rho[upper.tri(Rho)]
    var_r2 <- matrix(sapply(n, v2fun, x = R_u), ncol = l_n)
    d <- digit(p)
    a <- seq_len(p)
    A <- outer(as.integer(10 ^ d) * a, a, "+")
    b <- t(A)[lower.tri(A)]
    IJ <- rep.int(b, seq.int(length(b) - 1, 0))
    KL <- rep_d(b)
    Is <- as.integer((IJ %/% 10 ^ d) %% 10 ^ d)
    Js <- as.integer(IJ %% 10 ^ d)
    Ks <- as.integer((KL %/% 10 ^ d) %% 10 ^ d)
    Ls <- as.integer(KL %% 10 ^ d)
    rm(IJ, KL)
    exv_r1 <- matrix(sapply(n, e1fun, x = R_u), ncol = l_n)
    Eij <- exv_r1[(Js - 1) * p + Is - (Js - 1) * (2 * p - Js + 2) / 2, ]
    Ekl <- exv_r1[(Ls - 1) * p + Ks - (Ls - 1) * (2 * p - Ls + 2) / 2, ]
    Eij_Ekl <- Eij * Ekl
    rm(exv_r1, Eij, Ekl)
    Rij <- Rho[(Js - 1) * p + Is]
    Rkl <- Rho[(Ls - 1) * p + Ks]
    Rik <- Rho[(Ks - 1) * p + Is]
    Rjl <- Rho[(Ls - 1) * p + Js]
    Ril <- Rho[(Ls - 1) * p + Is]
    Rjk <- Rho[(Ks - 1) * p + Js]
    rm(Is, Js, Ks, Ls)
    Cijkl <- c1fun(n, Rij, Rkl, Rik, Rjl, Ril, Rjk)
    rm(Rij, Rkl, Rik, Rjl, Ril, Rjk)
    # cov_r2 <- 8 * Eij_Ekl * Cijkl + 4 * Cijkl ^ 2
    # ans <- colSums(var_r2) + colSums(cov_r2)
    cov_r2 <- 4 * (2 * crossprod(Eij_Ekl, Cijkl) + crossprod(Cijkl))
    ans <- colSums(var_r2) + diag(cov_r2)
    4 * ans / (p * (p - 1))^2
}

##### AVar.VRR_pfd #####
#' Approximate variance of relative eigenvalue variance of correlation matrix,
#' vectorized
#'
#' \code{AVar.VRR_pfd()}: further improvement over \code{AVar.VRR_pfv()}.
#'
#' @rdname AVar.VRR_xx
#'
AVar.VRR_pfd <- function(Rho, n = 100, Lambda, exv1.mode = c("exact", "asymptotic"),
                         var1.mode = "asymptotic",
                         var2.mode = c("exact", "asymptotic"),
                         order.exv1 = 2, order.var2 = 2,
                         mode = c("for.ind", "lapply", "mclapply", "parLapply"),
                         mc.cores = "auto", max.cores = parallel::detectCores(),
                         do.mcaffinity = TRUE,
                         affinity_mc = seq_len(max.cores), cl = NULL,
                         max.size = 2e6, bd = NULL,
                         verbose = c("no", "yes", "inline"), ...) {
    exv1.mode <- match.arg(exv1.mode)
    var1.mode <- match.arg(var1.mode)
    var2.mode <- match.arg(var2.mode)
    mode <- match.arg(mode)
    verbose <- match.arg(verbose)
    if(mode == "mclapply" || mode == "parLapply") {
        if(!requireNamespace("parallel")) {
            stop("Package 'parallel' is required for ",
                 "mode = 'mclapply' or 'parLapply'")
        }
    }
    rep_d <- function(x, from = 1, to = length(x)) {
        x <- x[seq.int(from, length(x))]
        sx <- seq_along(x)[seq_len(to - from + 1)]
        unlist(lapply(sx, function(i) x[-seq_len(i)]))
    }
    bd2I1 <- function(bd) {
        lbd <- length(bd)
        q <- sapply(bd, length)
        qr <- q[seq.int(lbd, 1)]
        p1 <- cumsum(qr)[seq.int(lbd, 1)]
        p2 <- p1 - q + 1
        rbind(p1, p2)
    }
    e1fun <- switch(exv1.mode,
                    exact = function(n, x) Exv.r1(n, x, ...),
                    asymptotic = function(n, x) AExv.r1(n, x,
                                                        order. = order.exv1))
    c1fun <- switch(var1.mode, function(n, Rij, Rkl, Rik, Rjl, Ril, Rjk) {
        A <- (Rij * Rkl * (Rik^2 + Ril^2 + Rjk^2 + Rjl^2) / 2 + Rik * Rjl +
             Ril * Rjk - (Rij * Rik * Ril + Rij * Rjk * Rjl + Rik * Rjk * Rkl +
                        Ril * Rjl * Rkl))
        outer(A, n, "/")
    })
    v2fun <- switch(var2.mode,
                    exact = function(n, x) Var.r2(n, x, do.square = TRUE, ...),
                    asymptotic = function(n, x) AVar.r2(n, x,
                                                        order. = order.var2))
    disp <- switch(verbose,
        yes = function(i) cat(paste0(" ", i, "/", lbd, "; ", Sys.time(), "\n")),
        inline = function(i) {
            if(i %% 10^floor(log(lbd, 10) - 1) == 0) cat(paste0(" ", i, "."))
            else invisible(NULL)},
        function(i) invisible(NULL))
    Cov_r2m <- function(i) {
        disp(i)
        IJ <- rep.int(bd[[i]], seq.int(I1[1, i], I1[2, i]) - 1)
        KL <- rep_d(b, I2[1, i], I2[2, i])
        Is <- as.integer((IJ %/% 10 ^ d) %% 10 ^ d)
        Js <- as.integer(IJ %% 10 ^ d)
        Ks <- as.integer((KL %/% 10 ^ d) %% 10 ^ d)
        Ls <- as.integer(KL %% 10 ^ d)
        # rm(IJ, KL)
        Eij <- exv_r1[(Js - 1) * p + Is - (Js - 1) * (2 * p - Js + 2) / 2, ]
        Ekl <- exv_r1[(Ls - 1) * p + Ks - (Ls - 1) * (2 * p - Ls + 2) / 2, ]
        Eij_Ekl <- Eij * Ekl
        # rm(Eij, Ekl)
        Rij <- Rho[(Js - 1) * p + Is]
        Rkl <- Rho[(Ls - 1) * p + Ks]
        Rik <- Rho[(Ks - 1) * p + Is]
        Rjl <- Rho[(Ls - 1) * p + Js]
        Ril <- Rho[(Ls - 1) * p + Is]
        Rjk <- Rho[(Ks - 1) * p + Js]
        # rm(Is, Js, Ks, Ls)
        Cijkl <- c1fun(n, Rij, Rkl, Rik, Rjl, Ril, Rjk)
        # return(colSums(8 * Eij_Ekl * Cijkl + 4 * Cijkl ^ 2))
        return(4 * diag(2 * crossprod(Eij_Ekl, Cijkl) + crossprod(Cijkl)))
        # rm(Eij_Ekl, Cijkl)
    }
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    l_n <- length(n)
    R_u <- Rho[upper.tri(Rho)]
    var_r2 <- matrix(sapply(n, v2fun, x = R_u), ncol = l_n)
    exv_r1 <- matrix(sapply(n, e1fun, x = R_u), ncol = l_n)
    d <- digit(p)
    if(is.null(bd)) {
        a <- seq_len(p)
        A <- outer(as.integer(10 ^ d) * a, a, "+")
        b <- t(A)[lower.tri(A)]
        bd <- divInd(b, Max = max.size)
    }
    I1 <- bd2I1(bd)
    I2 <- length(b) - I1 + 1
    lbd <- length(bd)
    cov_r2 <- numeric(l_n)
    if(verbose == "inline") cat(paste0(" Out of ", lbd, " iterations:"))
    if(mode == "for.ind") {
        for(i in seq_along(bd)) {
            cov_r2 <- cov_r2 + Cov_r2m(i)
        }
        ans <- colSums(var_r2) + cov_r2
    } else {
        if(mode == "lapply") {
            cov_r2 <- lapply(seq_along(bd), Cov_r2m)
        } else {
            # require(parallel)
            if(mc.cores == "auto") {
                mc.cores <- pmin(lbd, max.cores)
            }
            if(mode == "mclapply") {
                if(do.mcaffinity) {
                    try(invisible(parallel::mcaffinity(affinity_mc)),
                        silent = TRUE)
                }
                cov_r2 <- parallel::mclapply(seq_along(bd), Cov_r2m,
                                             mc.cores = mc.cores)
            } else if(mode == "parLapply") {
                if(is.null(cl)) {
                    cl <- parallel::makeCluster(mc.cores, outfile = "")
                    on.exit(parallel::stopCluster(cl))
                    if(do.mcaffinity) {
                        try(parallel::clusterApply(cl, as.list(affinity_mc),
                            parallel::mcaffinity), silent = TRUE)
                    }
                }
                cov_r2 <- parallel::parLapply(cl = cl, seq_along(bd), Cov_r2m)
            }
        }
        cov_r2 <- matrix(unlist(cov_r2), l_n)
        ans <- colSums(var_r2) + rowSums(cov_r2)
    }
    if(verbose == "inline") cat("\n")
    4 * ans / (p * (p - 1))^2
}

##### AVar.VRR_pfc #####
#' Approximate variance of relative eigenvalue variance of correlation matrix,
#' with Rcpp
#'
#' \code{AVar.VRR_pfc()}: fast version using \code{Rcpp}. Requires
#' the extension package \code{eigvaldispRcpp}.
#'
# #' When the function to be used (specified by cppfun) is not found,
# #' an attempt is made to sourceCpp() the .cpp file in the present directory.
# #' It is recommended to do, e.g., sourceCpp("Cov.r2C.cpp")
# #' to compile the C++ code beforehand as this typically takes several seconds.
#'
#' @rdname AVar.VRR_xx
#'
AVar.VRR_pfc <- function(Rho, n = 100, Lambda, cppfun = "Cov_r2C",
                         nthreads = 0L,
                         exv1.mode = c("exact", "asymptotic"),
                         # var1.mode = "asymptotic",
                         var2.mode = c("exact", "asymptotic"),
                         order.exv1 = 2, order.var2 = 2, ...) {
    if(!requireNamespace("eigvaldispRcpp")) {
        stop("Package 'eigvaldispRcpp' is required for AVar.VRR_pfc() to run.",
             "\n  Install it from github.com/watanabe-j/eigvaldispRcpp")
    }
    cppfun <- match.arg(cppfun, c("Cov_r2C", "Cov_r2A", "Cov_r2E", "Cov_r2P",
                                  "Armadillo", "Eigen", "Parallel"))
    exv1.mode <- match.arg(exv1.mode)
    # var1.mode <- match.arg(var1.mode)
    var2.mode <- match.arg(var2.mode)
    c1fun <- switch(cppfun,
                    Cov_r2C   = eigvaldispRcpp:::Cov_r2C,
                    Cov_r2A   = eigvaldispRcpp:::Cov_r2A,
                    Cov_r2E   = eigvaldispRcpp:::Cov_r2E,
                    Cov_r2P   = eigvaldispRcpp:::Cov_r2P,
                    Armadillo = eigvaldispRcpp:::Cov_r2A,
                    Eigen     = eigvaldispRcpp:::Cov_r2E,
                    Parallel  = eigvaldispRcpp:::Cov_r2P)
    e1fun <- switch(exv1.mode,
                    exact = function(n, x) Exv.r1(n, x, ...),
                    asymptotic = function(n, x) AExv.r1(n, x,
                                                        order. = order.exv1))
    v2fun <- switch(var2.mode,
                    exact = function(n, x) Var.r2(n, x, do.square = TRUE, ...),
                    asymptotic = function(n, x) AVar.r2(n, x,
                                                        order. = order.var2))
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    l_n <- length(n)
    R_u <- Rho[upper.tri(Rho)]
    var_r2 <- matrix(sapply(n, v2fun, x = R_u), ncol = l_n)
    exv_r1 <- matrix(sapply(n, e1fun, x = R_u), ncol = l_n)
    cov_r2 <- c(c1fun(n, Rho, exv_r1, nthreads))
    ans <- colSums(var_r2) + cov_r2
    4 * ans / (p * (p - 1))^2
}

##### AVar.VRR_kl #####
#' Asymptotic variance of relative eigenvalue variance of correlation matrix
#'
#' \code{AVar.VRR_kl()}: asymptotic variance from Konishi's theory:
#' \eqn{Vrel(R)} as function of eigenvalues.
#'
#' @rdname AVar.VRR_xx
#'
AVar.VRR_kl <- function(Rho, n = 100, Lambda, ...) {
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    svd.Rho <- svd(Rho, nu = 0)
    d <- svd.Rho$d
    V2 <- svd.Rho$v ^ 2
    R2 <- Rho ^ 2
    ans <- numeric(1)
    for(i in 1:p) {
        for(j in 1:p) {
            ans <- ans + d[i]^2 * d[j]^2 * drop(as.numeric(i == j)
                   - (d[i] + d[j]) * crossprod(V2[, i], V2[, j])
                   + crossprod(V2[, i], crossprod(R2, V2[, j])))
        }
    }
    8 * ans / (p * (p - 1)) ^ 2 / n
}

##### AVar.VRR_klv #####
#' Asymptotic variance of relative eigenvalue variance of correlation matrix
#'
#' \code{AVar.VRR_klv()}: vectorized version of \code{AVar.VRR_kl()}.
#'
#' @rdname AVar.VRR_xx
#'
AVar.VRR_klv <- function(Rho, n = 100, Lambda, ...) {
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    svd.Rho <- svd(Rho, nu = 0)
    d <- svd.Rho$d
    d2 <- d ^ 2
    V2 <- svd.Rho$v ^ 2
    R2 <- Rho ^ 2
    G <- diag(p) - outer(d, d, "+") * crossprod(V2) +
         crossprod(V2, crossprod(R2, V2))
    F <- (d2 %*% t(d2)) * G
    ans <- sum(F)
    8 * ans / (p * (p - 1)) ^ 2 / n
}

##### AVar.VRR_kr #####
#' Asymptotic variance of relative eigenvalue variance of correlation matrix
#'
#' \code{AVar.VRR_kr()}: asymptotic variance from Konishi's theory:
#' \eqn{Vrel(R)} as function of correlation coefficients
#'
#' @rdname AVar.VRR_xx
#'
AVar.VRR_kr <- function(Rho, n = 100, Lambda, ...) {
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    ans <- numeric(1)
    for(i in 1:p) {
        for(j in (1:p)[-i]) {
            for(k in i:p) {
                for(l in (1:p)[-k]) {
                    ans <- ans + (Rho[j, k] - Rho[i, j] * Rho[i, k]) *
                                 (Rho[i, l] - Rho[i, k] * Rho[k, l]) *
                                 Rho[i, j] * Rho[k, l]
                }
            }
        }
    }
    16 * ans / (p * (p - 1)) ^ 2 / n
}

##### AVar.VRR_krv #####
#' Asymptotic variance of relative eigenvalue variance of correlation matrix
#'
#' \code{AVar.VRR_krv()}: vectorized version of \code{AVar.VRR_kr()}.
#'
#' @rdname AVar.VRR_xx
#'
AVar.VRR_krv <- function(Rho, n = 100, Lambda, ...) {
    # digit <- function(x) {
    #     i <- 1L
    #     while(x %/% 10 >= 1) {
    #         x <- x %/% 10
    #         i <- i + 1
    #     }
    #     return(i)
    # }
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    d <- digit(p)
    a <- seq_len(p)
    A <- outer(as.integer(10 ^ d) * a, a, "+")
    b <- A[lower.tri(A) | upper.tri(A)]
    B <- matrix(rep.int(b, length(b)), length(b))
    IJ <- t(B)[lower.tri(B)]
    KL <- B[lower.tri(B)]
    Is <- as.integer((IJ %/% 10 ^ d) %% 10 ^ d)
    Js <- as.integer(IJ %% 10 ^ d)
    Ks <- as.integer((KL %/% 10 ^ d) %% 10 ^ d)
    Ls <- as.integer(KL %% 10 ^ d)
    Rij <- Rho[(Js - 1) * p + Is]
    Rkl <- Rho[(Ls - 1) * p + Ks]
    Rik <- Rho[(Ks - 1) * p + Is]
    Rjl <- Rho[(Ls - 1) * p + Js]
    Ril <- Rho[(Ls - 1) * p + Is]
    Rjk <- Rho[(Ks - 1) * p + Js]
    ans <- (Rjk - Rij * Rik) * (Ril - Rik * Rkl) * Rij * Rkl
    ans <- sum(ans)
    16 * ans / (p * (p - 1)) ^ 2 / n
}

##### Exv.VXXa #####
#' Moments of ``bias-corrected'' eigenvalue dispersion indices
#'
#' Functions to calculate expectation/variance of eigenvalue dispersion indices
#' of covariance/correlation matrices.
#'
#' Usage is identical to that of the corresponding unadjusted versions
#' (see \code{\link{Exv.VXX}}), which are in most cases called internally.
#'
#' @name Exv.VXXa
#'
#' @inheritParams Exv.VXX
#'
#' @return
#' A numeric vector of the desired moment, corresponding to \code{n}.
#'
#' @references
#' Watanabe, J. (2022). Statistics of eigenvalue dispersion indices:
#'  quantifying the magnitude of phenotypic integration. *Evolution*,
#'  **76**, 4--28. doi:[10.1111/evo.14382](https://doi.org/10.1111/evo.14382).
#'
#' @seealso
#' \code{\link{VXXa}} for ``bias-corrected'' estimators
#'
#' \code{\link{Exv.VXX}} for moments of unajusted versions
#'
#' @examples
#' # See also examples of Exv.VXX
#' # Covariance matrix
#' N <- 20
#' Lambda <- c(4, 2, 1, 1)
#' (Sigma <- GenCov(evalues = Lambda, evectors = "random"))
#' VE(S = Sigma)$VE
#' VE(S = Sigma)$VR
#' # Population values of V(Sigma) and Vrel(Sigma)
#'
#' # Moments of bias-corrected eigenvalue variance of covariance matrix
#' Exv.VESa(Sigma, n = N - 1)
#' Var.VESa(Sigma, n = N - 1)
#' # The expectation is equal to the population value (as it should be)
#'
#' # Moments of adjusted relative eigenvalue variance of covariance matrix
#' Exv.VRSa(Sigma, n = N - 1)
#' Var.VRSa(Sigma, n = N - 1)
#' # Slight underestimation is expected
#' # All these are the same with Lambda = Lambda is specified instead of Sigma.
#'
#' # Correlation matrix
#' (Rho <- GenCov(evalues = Lambda / sum(Lambda) * 4, evectors = "Givens"))
#' VE(S = Rho)$VR
#' # Population value of Vrel(Rho), identical to Vrel(Sigma) as it should be
#'
#' Exv.VRRa(Rho, n = N - 1)
#' Var.VRRa(Rho, n = N - 1)
#' # Slight underestimation is expected
#' # These results vary with the choice of eigenvalues
#' # If interested, repeat from the definition of Rho
#'
#' # All options for Var.VRR() are accommodated
#' Var.VRRa(Rho, n = N - 1, fun = "pfd") # Pan-Frank method; default
#' Var.VRRa(Rho, n = N - 1, fun = "klv") # Konishi's theory
#'
NULL

##### Exv.VESa #####
#' Expectation of bias-corrected eigenvalue variance of covariance matrix
#'
#' \code{Exv.VESa()}: expectation of unbiased eigenvalue variance of
#' covariance matrix. Of little practical use because
#' this is just the population value \eqn{V(\Sigma)}.
#'
#' @rdname Exv.VXXa
#'
#' @export
#'
Exv.VESa <- function(Sigma, n = 100, Lambda, drop_0 = FALSE,
                     tol = .Machine$double.eps * 100, ...) {
    if(missing(Lambda)) {
        Lambda <- svd(Sigma, nu = 0, nv = 0)$d
    }
    if(drop_0) {
        Lambda <- Lambda[Lambda > tol]
    } else {
        Lambda[Lambda < tol] <- 0
    }
    p <- length(Lambda)
    t1 <- sum(Lambda)
    t2 <- sum(Lambda ^ 2)
    t2 / p - t1 ^ 2 / p^2
}

##### Exv.VRSa #####
#' Expectation of adjusted relative eigenvalue variance of covariance matrix
#'
#' \code{Exv.VRSa()}: expectation of adjusted relative eigenvalue variance
#' of covariance matrix.
#'
#' @rdname Exv.VXXa
#'
#' @export
#'
Exv.VRSa <- function(Sigma, n = 100, Lambda, drop_0 = FALSE,
                     tol = .Machine$double.eps * 100, ...) {
    if(missing(Lambda)) {
        Lambda <- svd(Sigma, nu = 0, nv = 0)$d
    }
    if(drop_0) {
        Lambda <- Lambda[Lambda > tol]
    } else {
        Lambda[Lambda < tol] <- 0
    }
    p <- length(Lambda)
    E <- Exv.VRS(Lambda = Lambda, n = n, drop_0 = drop_0, tol = tol, ...)
    En <- (p + 2) / (p * n + 2)
    1 - (1 - E) / (1 - En)
}

##### Exv.VRRa #####
#' Expectation of adjusted relative eigenvalue variance of correlation matrix
#'
#' \code{Exv.VRRa()}: expectation of adjusted relative eigenvalue variance
#' of correlation matrix.
#'
#' @rdname Exv.VXXa
#'
#' @export
#'
Exv.VRRa <- function(Rho, n = 100, Lambda, tol = .Machine$double.eps * 100,
                     tol.hg = 0, maxiter.hg = 2000, ...) {
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    E <- Exv.VRR(Rho = Rho, n = n, tol = tol,
                 tol.hg = tol.hg, maxiter.hg = maxiter.hg, ...)
    En <- 1 / n
    1 - (1 - E) / (1 - En)
}

##### Var.VESa #####
#' Variance of bias-corrected eigenvalue variance of covariance matrix
#'
#' \code{Var.VESa()}: variance of unbiased eigenvalue variance of
#' covariance matrix.
#'
#' @rdname Exv.VXXa
#'
#' @export
#'
Var.VESa <- function(Sigma, n = 100, Lambda, drop_0 = FALSE,
                     tol = .Machine$double.eps * 100, ...) {
    if(missing(Lambda)) {
        Lambda <- svd(Sigma, nu = 0, nv = 0)$d
    }
    if(drop_0) {
        Lambda <- Lambda[Lambda > tol]
    } else {
        Lambda[Lambda < tol] <- 0
    }
    p <- length(Lambda)
    t1 <- sum(Lambda)
    t2 <- sum(Lambda ^ 2)
    t3 <- sum(Lambda ^ 3)
    t4 <- sum(Lambda ^ 4)
    4 * ((2 * p^2 * n^2 + 3 * p^2 * n - 6 * p^2 - 4 * p * n - 4) * t4 +
         - 4 * p * (n - 1) * (n + 2) * t3 * t1 +
         (p^2 * n + 4 * p + 2 * n + 2) * t2^2 +
         2 * (n - 1) * (n + 2) * t2 * t1^2) / (p ^ 4 * n * (n - 1) * (n + 2))
}

##### Var.VRSa #####
#' Variance of adjusted relative eigenvalue variance of covariance matrix
#'
#' \code{Var.VRSa()}: variance of adjusted relative eigenvalue variance of
#' covariance matrix.
#'
#' @rdname Exv.VXXa
#'
#' @export
#'
Var.VRSa <- function(Sigma, n = 100, Lambda, drop_0 = FALSE,
                 tol = .Machine$double.eps * 100, ...) {
    # divisor <- match.arg(divisor)
    if(missing(Lambda)) {
        Lambda <- svd(Sigma, nu = 0, nv = 0)$d
    }
    if(drop_0) {
        Lambda <- Lambda[Lambda > tol]
    } else {
        Lambda[Lambda < tol] <- 0
    }
    p <- length(Lambda)
    V <- Var.VRS(Lambda = Lambda, n = n, drop_0 = drop_0, tol = tol, ...)
    En <- (p + 2) / (p * n + 2)
    V / (1 - En)^2
}

##### Var.VRRa #####
#' Variance of adjusted relative eigenvalue variance of covariance matrix
#'
#' \code{Var.VRRa()}: variance of adjusted relative eigenvalue variance
#' of correlation matrix.
#'
#' @rdname Exv.VXXa
#'
#' @export
#'
Var.VRRa <- function(Rho, n = 100, Lambda, tol = .Machine$double.eps * 100,
                     tol.hg = 0, maxiter.hg = 2000, ...) {
    if(missing(Rho)) {
        Rho <- GenCov(evalues = Lambda, evectors = "Givens")
        if(Lambda[2] != Lambda[length(Lambda)]) {
            warning("Rho was generated from the eigenvalues provided \n  ",
                    "Expectation may vary even if eigenvalues are fixed")
        }
    } else if(any(diag(Rho) != 1)) {
        Rho <- Rho / Rho[1, 1]
        if(any(diag(Rho) != 1)) {
            stop("Provide a valid correlation matrix, or its eigenvalues")
        }
        warning("Rho was scaled to have diagonal elements of unity")
    }
    p <- ncol(Rho)
    V <- Var.VRR(Rho = Rho, n = n, tol = tol,
                 tol.hg = tol.hg, maxiter.hg = maxiter.hg, ...)
    En <- 1 / n
    V / (1 - En)^2
}
