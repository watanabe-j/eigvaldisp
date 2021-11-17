##### VE #####
#' Calculate eigenvalue dispersion indices
#'
#' Function to calculate eigenvalue variance \eqn{V} and
#' relative eigenvalue variance \eqn{V_{rel}}
#' of a covariance/correlation matrix, either from a data matrix \code{X}
#' or a covariance/correlation matrix \code{S}.
#'
#' Provide either a data matrix (\code{X}), covariance/correlation matrix
#' (\code{S}), or vector of eigenvalues (\code{L}).
#' When \code{X} is given, the default divisor is \eqn{N - 1} where \eqn{N}
#' is sample size.
#'
#' Sometimes it might be desirable to evaluate eigenvalue dispersion
#' in a selected subspace, rather than in the full space.
#' For this, provide the argument \code{sub} to restrict calculations to
#' the subspace corresponding to the specified eigenvalues/vectors.
#' Alternatively, set \code{drop_0 = TRUE} to drop
#' zero eigenvalues from calculation.
#' The former way would be more useful when the subspace of interest is known
#' a priori. The latter is ad hoc, automatically dropping zero eigenvalues
#' whose magnitudes are below the specified tolerance.
#'
#' @param X
#'   Data matrix from which covariance/corrrelation matrix is obtained.
#' @param S
#'   Covariance/correlation matrix.
#' @param L
#'   Vector of eigenvalues.
#' @param center
#'   Logical to specify whether sample-mean-centering should be done.
#' @param scale.
#'   Logical to specify whether SD-scaling should be done
#'   (that is, when \code{TRUE}, the analysis is on the correlation matrix).
#'   When \code{S} is provided (but \code{X} is not), this is converted to
#'   a correlation matrix.
#' @param divisor
#'   Either \code{"UB"} (default) or \code{"ML"},
#'   to decide the default value of \code{m}.
#' @param m
#'   Divisor for the sample covariance matrix (\eqn{n*} in Watanabe (2021)).
#' @param nv
#'   Numeric. Specify how many eigenvectors are to be retained; default 0.
#' @param sub
#'   Numeric/integer vector to specify the range of eigenvalue indices
#'   to be involved; used to exclude some subspace.
#' @param drop_0
#'   Logical, when \code{TRUE}, eigenvalues smaller than \code{tol} are dropped.
#' @param tol
#'   Tolerance to be used with \code{drop_0}.
#'
#' @return
#'   A list containing the following:
#'   \describe{
#'     \item{VE}{Eigenvalue variance (\eqn{V})}
#'     \item{VR}{Relative eigenvalue variance (\eqn{V_{rel}})}
#'     \item{meanL}{Mean (average) of the eigenvalues}
#'     \item{L}{Vector of eigenvalues}
#'     \item{U}{Matrix of eigenvectors, only when \code{nv > 0}}
#'   }
#'
#' @references
#' Cheverud, J. M., Rutledge, J. J., & Atchley, W. R. (1983). Quantitative
#'  genetics of development: genetic correlations among age-specific trait
#'  values and the evolution of ontogeny. *Evolution* **37**, 5--42.
#'  doi:[10.1111/j.1558-5646.1983.tb05619.x](https://doi.org/10.1111/j.1558-5646.1983.tb05619.x).
#'
#' Haber, A. (2011). A comparative analysis of integration indices.
#'  *Evolutionary Biology* **38**, 476--488.
#'  doi:[10.1007/s11692-011-9137-4](https://doi.org/10.1007/s11692-011-9137-4).
#'
#' Pavlicev, M., Cheverud, J. M., & Wagner, G. P. (2009). Measuring
#'  morphological integration using eigenvalue variance. *Evolutionary Biology*
#'  **36**, 157--170.
#'  doi:[10.1007/s11692-008-9042-7](https://doi.org/10.1007/s11692-008-9042-7).
#'
#' Van Valen, L. (1974). Multivariate structural statistics in natural history.
#'  *Journal of Theoretical Biology* **45**, 235--247.
#'  doi:[10.1016/0022-5193(74)90053-8](https://doi.org/10.1016/0022-5193(74)90053-8).
#'
#' Wagner, G. P. (1984). On the eigenvalue distribution of genetic and
#'  phenotypic dispersion matrices: evidence for a nonrandom organization
#'  of quantitative character variation. *Journal of Mathematical Biology*
#'  **7**, 77--95. doi:[10.1007/BF00275224](https://doi.org/10.1007/BF00275224).
#'
#' Watanabe, J. (2021). Statistics of eigenvalue dispersion indices:
#'  quantifying the magnitude of phenotypic integration. *Evolution*,
#'  doi:[10.1111/evo.14382](https://doi.org/10.1111/evo.14382).
#'
#' @importFrom stats cov2cor
#'
#' @export
#'
#' @seealso
#' \link{Exv.VXX} for moments of sample eigenvalue dispersion indices;
#' \link{VXXa} for bias-corrected versions.
#'
#' @examples
#' # For a population covariance matrix or population eigenvalues
#' set.seed(6835)
#' Lambda <- c(4, 2, 1, 1)
#' (Sigma <- GenCov(evalues = Lambda, evectors = "random"))
#' VE(L = Lambda)
#' VE(S = Sigma) # Same
#'
#' # For a random sample, sample covariance matrix or its eigenvalues
#' N <- 20
#' X <- eigvaldisp:::rmvn(N = N, Sigma = Sigma)
#' S <- cov(X)
#' L <- eigen(S)$values
#' VE(X = X)
#' VE(S = S) # Same
#' VE(L = L) # Same
#' # Thus, providing X is usually the most straightforward for a sample
#' # (and this is usally quicker for p > 30-50 or so).
#' # Also, observe bias in these quantities compared to population values.
#'
#' # Same for maximum likelihood estimator (divisor m = N)
#' VE(X = X, divisor = "ML")
#' VE(X = X, m = N)        # Same, but any divisor can be specified
#' VE(S = S * (N - 1) / N) # Same
#' # L and meanL are (N - 1) / N times the above,
#' # VE is ((N - 1) / N) ^ 2 times the above, whereas VR remains the same.
#'
#' # For a sample correlation matrix
#' R <- cor(X)
#' VE(S = R)
#' VE(X = X, scale. = TRUE) # Same, hence is usually quicker when p is large
#'
#' # Interested in eigenvectors?
#' VE(X = X, nv = 2)
#'
#' # Singular covariance matrix
#' Lambda2 <- c(4, 2, 1, 0)
#' (Sigma2 <- GenCov(evalues = Lambda2, evectors = "random"))
#' VE(S = Sigma2)                # Calculated in the full space
#' VE(S = Sigma2, sub = 1:3)     # In the subspace of the first 3 PCs
#' VE(S = Sigma2, drop_0 = TRUE) # Dropping zero eigenvalues (same in this case)
#'
#' # Sample from singular covariance
#' X2 <- eigvaldisp:::rmvn(N = N, Sigma = Sigma2, sqrt_method = "pivot")
#' VE(X = X2)                    # In the full space
#' VE(X = X2, sub = 1:3)         # In the subspace of the first 3 PCs
#' VE(X = X2, drop_0 = TRUE)     # Practically the same
#'
#' # Just to note, the null space is identical between the population and sample
#' # in this case (where N - 1 > p)
#' eigen(Sigma2)$vectors[, 4]
#' eigen(cov(X2))$vectors[, 4]
#'
#' # This is of course not the case when N - 1 < p, although
#' # a sample null space always encompasses the population null space.
#' Lambda3 <- 9:0
#' Sigma3 <- GenCov(evalues = Lambda3, evectors = "random")
#' X3 <- eigvaldisp:::rmvn(N = 6, Sigma = Sigma3, sqrt_method = "pivot")
#' eigS3 <- eigen(cov(X3))
#' (Popul_null <- eigen(Sigma3)$vectors[, Lambda3 < 1e-12])
#' (Sample_null <-eigS3$vectors[, eigS3$values < 1e-12])
#' crossprod(Popul_null, Sample_null)
#' # None of these vectors are identical, but
#' tcrossprod(crossprod(Popul_null, Sample_null))
#' # sum of squared cosines equals 1, as expected
#'
VE <- function(X, S, L, center = TRUE, scale. = FALSE,
               divisor = c("UB", "ML"), m = switch(divisor, UB = N - 1, ML = N),
               nv = 0, sub = seq_len(length(L)),
               drop_0 = FALSE, tol = .Machine$double.eps * 100) {
    divisor <- match.arg(divisor)
    if(!missing(X)) {
        X <- scale2(X, center = center, scale = scale.)
        N <- nrow(X)
        p <- ncol(X)
        svd.X <- svd(X, nu = 0, nv = nv)
        L <- svd.X$d
        L <- (L ^ 2)
        if(scale.) {
            L <- L / (N - 1)
        } else {
            L <- L / m
        }
        L <- c(L, rep_len(0, max(p - length(L), 0)))
    } else if(!missing(S)) {
        p <- ncol(S)
        if(scale.) S <- cov2cor(S)
        svd.X <- svd(S, nu = 0, nv = nv)
        L <- svd.X$d
    }
    L <- L[sub]
    if(drop_0) L <- L[L > tol]
    p <- length(L)
    mL <- mean(L)
    VE <- sum((L - mL) ^ 2) / p   # Faster than var(L) * (p - 1) / p,
    VR <- VE / ((p - 1) * mL ^ 2) # as mL is already available
    ans <- list(VE = VE, VR = VR, meanL = mL, L = L)
    if(nv > 0) ans <- c(ans, list(U = svd.X$v))
    return(ans)
}


##### VXXa #####
#' ``Bias-corrected'' eigenvalue dispersion indices
#'
#' These functions calculate bias-corrected or adjusted
#' eigenvalue dispersion indices
#'
#' These functions are to be used with sample data matrix,
#' covariance/correlation matrix, or eigenvalues, and not to be used with
#' population quantities. Unless \code{X} is provided, the degrees of freedom
#' \code{n} should be specified as this is required for adjustment.
#'
#' These functions internally call \code{VE(L = L, nv = 0, ...)} with
#' appropriately constructed (or user-specified) \code{L},
#' and adjusted eigenvalue dispersion indices are appended to the outcome list.
#' If \code{nv > 0}, eigenvectors are calculated before this function call
#' and returned as well.
#'
#' Bias correction is possible for eigenvalue variance of covariance matrices,
#' but not straightforward for relative eigenvalue variance (of either
#' covariance or correlation matrices), because the latter is a nonlinear
#' function of the population value (see Watanabe, 2021).
#' Although this could potentially be achievable for correlation matrices,
#' expression for its variance is not known to the author.
#'
#' The bias-corrected eigenvalue variance of a sample covariance matrix
#' also has smaller variance than the uncorrected version.
#' On the other hand, the adjusted relative eigenvalue variance of a
#' sample covariance/correlation matrix has larger variance than
#' the unadjusted version.
#'
#' @name VXXa
#'
#' @inheritParams VE
#'
#' @param n
#'   Degrees of freedom; required unless X is provided.
#' @param ...
#'   Arguments \code{sub}, \code{drop_0}, \code{tol} can be
#'   passed to \code{VE()}. Other arguments will not influence the result.
#'
#' @return
#' A list similar to the output of \code{VE()}, with an additional element:
#' \describe{
#'   \item{\code{VESa()}}{\code{$VESa}: Bias-corrected eigenvalue
#'     variance of covariance matrix.}
#'   \item{\code{VRSa()}}{\code{$VRSa}: Adjusted relative
#'     eigenvalue variance of covariance matrix}
#'   \item{\code{VRRa()}}{\code{$VRRa}: Adjusted relative
#'     eigenvalue variance of correlation matrix}
#' }
#'
#' @seealso
#' \link{VE} for the main function;
#' \link{Exv.VXX} for expectation (bias) in the ordinary estimators;
#' \link{Exv.VXXa} for the expectation/variance of the adjusted estimators.
#'
#' @references
#' Watanabe, J. (2021). Statistics of eigenvalue dispersion indices:
#'  quantifying the magnitude of phenotypic integration. *Evolution*,
#'  doi:[10.1111/evo.14382](https://doi.org/10.1111/evo.14382).
#'
#' @examples
#' # Spherical covariance matrix (VE = VR = 0)
#' Sigma <- diag(4)
#' VE(S = Sigma)
#'
#' N <- 20
#' set.seed(375)
#' X <- eigvaldisp:::rmvn(N = N, Sigma = Sigma)
#' S <- cov(X)
#' L <- eigen(S)$values
#'
#' # Bias-corrected eigenvalue variance of covariance matrix
#' VESa(X = X)
#' VESa(S = S, n = N - 1)$VESa # Same, but n is required
#' VESa(L = L, n = N - 1)$VESa # Same
#' # Note the overestimation in VE and VR,
#' # and (slightly) better performance of VESa
#' # (although it does not always work this well)
#'
#' # Adjusted relative eigenvalue variance of covariance matrix
#' VRSa(X = X)$VRSa
#'
#' # Population value for the correlation matrix (same to Sigma in this case)
#' VE(S = stats::cov2cor(Sigma))
#'
#' # Adjusted relative eigenvalue variance of correlation matrix
#' VRRa(X = X)
#'
#'
#' # Covariance/correlation matrix with strong covariation (VR = 0.8)
#' (Sigma2 <- GenCov(p = 4, VR = 0.8, evectors = "Givens"))
#' VE(S = Sigma2)
#'
#' N <- 20
#' X2 <- eigvaldisp:::rmvn(N = N, Sigma = Sigma2)
#'
#' VESa(X = X2)
#' # This is better than the un-corrected version
#'
#' VRSa(X = X2)
#' VRRa(X = X2)
#' # But these are slightly worse (as expected)
#'
NULL

##### VESa #####
#' Bias-corrected eigenvalue variance of covariance matrix
#'
#' \code{VESa}: bias-corrected eigenvalue variance of covariance matrix.
#'
#' @rdname VXXa
#'
#' @export
#'
VESa <- function(X, S, L, n = N - as.numeric(center), divisor = c("UB", "ML"),
                 m = switch(divisor, UB = N - 1, ML = N),
                 center = TRUE, nv = 0, ...) {
    divisor <- match.arg(divisor)
    if(!missing(X)) {
        X <- scale2(X, center = center, scale = FALSE)
        N <- nrow(X)
        p <- ncol(X)
        svd.X <- svd(X, nu = 0, nv = nv)
        L <- svd.X$d
        L0 <- (L ^ 2)
        L <- L0 / m
        L <- c(L, rep_len(0, max(p - length(L), 0)))
    } else {
        if(missing(n)) stop("Provide n (or X)")
        if(!missing(S)) {
            p <- ncol(S)
            # if(scale.) S <- cov2cor(S)
            svd.X <- svd(S, nu = 0, nv = nv)
            L <- svd.X$d
        }
        L0 <- L * n
    }
    ans <- VE(L = L, nv = 0, ...)
    p <- length(ans$L)
    t1 <- sum(L0)
    t2 <- sum(L0 ^ 2)
    ans$VESa <- ((p * n + 2) * t2 - (p + n + 1) * t1 ^ 2) /
                    (p ^ 2 * n * (n - 1) * (n + 2))
    if(nv > 0) ans <- c(ans, list(U = svd.X$v))
    return(ans)
}

##### VRSa #####
#' Adjusted relative eigenvalue variance of covariance matrix
#'
#' \code{VRSa}: adjusted relative eigenvalue variance of covariance matrix.
#'
#' @rdname VXXa
#'
#' @export
#'
VRSa <- function(X, S, L, n = N - as.numeric(center), divisor = c("UB", "ML"),
                 m = switch(divisor, UB = N - 1, ML = N),
                 center = TRUE, nv = 0, ...) {
    divisor <- match.arg(divisor)
    if(!missing(X)) {
        X <- scale2(X, center = center, scale = FALSE)
        N <- nrow(X)
        p <- ncol(X)
        svd.X <- svd(X, nu = 0, nv = nv)
        L <- svd.X$d
        L <- (L ^ 2) / m
        L <- c(L, rep_len(0, max(p - length(L), 0)))
    } else {
        if(missing(n)) stop("Provide n (or X)")
        if(!missing(S)) {
            p <- ncol(S)
            # if(scale.) S <- cov2cor(S)
            svd.X <- svd(S, nu = 0, nv = nv)
            L <- svd.X$d
        }
    }
    ans <- VE(L = L, nv = 0, ...)
    p <- length(ans$L)
    ans$VRSa <- (ans$VR * (p * n + 2) - (p + 2)) / p / (n - 1)
    if(nv > 0) ans <- c(ans, list(U = svd.X$v))
    return(ans)
}


##### VRRa #####
#' Adjusted relative eigenvalue variance of correlation matrix
#'
#' \code{VRRa}: adjusted relative eigenvalue variance of correlation matrix.
#'
#' @rdname VXXa
#'
#' @export
#'
VRRa <- function(X, S, L, n = N - as.numeric(center),
                 center = TRUE, nv = 0, ...) {
    if(!missing(X)) {
        X <- scale2(X, center = center, scale = TRUE)
        N <- nrow(X)
        p <- ncol(X)
        svd.X <- svd(X, nu = 0, nv = nv)
        L <- svd.X$d
        L <- (L ^ 2) / (N - 1)
        L <- c(L, rep_len(0, max(p - length(L), 0)))
    } else {
        if(missing(n)) stop("Provide n (or X)")
        if(!missing(S)) {
            p <- ncol(S)
            if(TRUE) S <- cov2cor(S)
            svd.X <- svd(S, nu = 0, nv = nv)
            L <- svd.X$d
        }
    }
    ans <- VE(L = L, nv = 0, ...)
    # p <- length(ans$L)
    ans$VRRa <- (ans$VR * n - 1) / (n - 1)
    if(nv > 0) ans <- c(ans, list(U = svd.X$v))
    return(ans)
}



##### simulateVE #####
#' Simulate eigenvalue dispersion indices
#'
#' \code{simulateVE()} iteratively generates multivariate normal variates,
#' from which eigenvalues, eigenvectors (optional, if \code{nv > 0}),
#' and eigenvalue dispersion indices
#' (both unstandardized \eqn{V} and standardized \eqn{V_{rel}}) of
#' sample covariance and correlation matrices are obtained.
#' These are returned as an invisible list.
#'
#' When \code{simulateVE()} is called, either data (\code{X}),
#' a covariance matrix (\code{Sigma}),
#' or its Cholesky factor (\code{cSigma}) should be given.
#' Specify the sample size \code{N} as well, unless \code{X} is provided.
#'
#' \code{simulateVE()} does not actually calculate sample covariance/correlation
#' matrices, but instead directly obtain eigenvalues and eigenvectors of these
#' from singular value decomposition \code{svd()} of normal variates
#' generated with \code{rmvn()}.
#'
#' In the output of \code{simulateVE()}, suffices \code{*v} and \code{*r}
#' denote covariance and correlation matrices, respectively.
#'
# #' @inheritParams rmvn
# #' @inheritParams VE
#' @param b
#'   Number of iterations.
#' @param X
#'   (Optional) Data matrix; when \code{Sigma} and \code{cSigma} are missing,
#'   the sample covariance matrix from \code{X} is used as
#'   the population covariance matrix in simulations (parametric bootstrapping).
# #' @param Sigma,cSigma,N   passed to rmvn
#' @param divisor,m,nv,drop_0,tol,center,scale.,sub
#'   These arguments are passed to \code{VE()}.
#'   \code{center} and \code{scale.} are also used construct \code{Sigma}
#'   when \code{X} is provided.
#'
#' @return
#'   \code{simulateVE()} invisiblly returns a list containing the following:
#'   \describe{
#'     \item{VEv}{Eigenvalue variance of cov matrix V(S) (b vector)}
#'     \item{VRv}{Relative eigenvalue variance of cov matrix Vrel(S) (b vector)}
#'     \item{Lv}{Eigenvalues of cov matrix (p * b matrix)}
#'     \item{VEr}{Eigenvalue variance of cor matrix V(R) (b vector)}
#'     \item{VRr}{Relative eigenvalue variance of cor matrix Vrel(R) (b vector)}
#'     \item{Lr}{Eigenvalues of cor matrix (p * b matrix)}
#'     (the rest are returned only when \code{nv > 0})
#'     \item{Uv}{Eigenvectors of cov matrix (p * nv * b array)}
#'     \item{Ur}{Eigenvectors of cor matrix (p * nv * b array)}
#'     \item{Uv.org}{Eigenvectors of population cov matrix (p * nv matrix)}
#'     \item{Ur.org}{Eigenvectors of population cor matrix (p * nv matrix)}
#'     (\code{Uv.org} and \code{Ur.org} are provided
#'     as references for signs of eigenvectors)
#'   }
#' @seealso \link{VE}, \link{sqrt_methods}
#'
#' @importFrom stats cov2cor
#'
#' @export
#'
#' @examples
#' Lambda <- c(4, 2, 1, 1)
#' (Sigma <- GenCov(evalues = Lambda, evectors = "random"))
#' N <- 10
#' set.seed(35638)
#' X1 <- eigvaldisp:::rmvn(N = N, Sigma = Sigma)
#' cSigma <- chol(Sigma)
#' set.seed(35638)
#' X2 <- eigvaldisp:::rmvn(N = N, cSigma = cSigma)
#' stopifnot(all.equal(X1, X2))
#' # These are identical. Providing cSigma is quicker as it skips sqrtfun(Sigma),
#' # so this may be useful when multiple simulations are to be run with
#' # the same population covariance matrix.
#'
#' # A small simulation
#' sim_result <- simulateVE(b = 50L, Sigma = Sigma, N = N)
#' str(sim_result)
#' # Results are returned as a list (see "Value" for details)
#'
#' # Syntax for parametric bootstrapping
#' param_boot <- simulateVE(b = 50L, X = X1)
#'
simulateVE <- function(b = 100L, X, Sigma, cSigma = sqrtfun(Sigma), N = nrow(X),
                   divisor = c("UB", "ML"),
                   m = switch(divisor, UB = N - 1, ML = N),
                   center = TRUE, scale. = FALSE, sub = seq_len(ncol(cSigma)),
                   nv = 0, sqrt_method = "chol",
                   drop_0 = FALSE, tol = .Machine$double.eps * 100) {
    sqrt_method <- match.arg(sqrt_method, c("default", "chol", "chol_piv",
                                   "chol_qr", "matsqrt", "pivot", "qr", "sqrt"))
    divisor <- match.arg(divisor)
    ## If X is provided (and Sigma and cSigma are missing), sample cov matrix is used
    if(missing(Sigma) && missing(cSigma) && !missing(X)) {
        X <- scale2(X, center = center, scale = scale.)
        Sigma <- crossprod(X) / m
    }
    ## Only cSigma is used for simulations (rather than X and Sigma)
    if(missing(cSigma)) {
        sqrtfun <- switch(sqrt_method, matsqrt = matsqrt, sqrt = matsqrt,
                          chol_qr = chol_qr, qr = chol_qr,
                          pivot = chol_piv, chol_piv = chol_piv, chol)
        cSigma <- sqrtfun(Sigma)
    } else if(!missing(Sigma) || !missing(X)) {
        warning("The arguments X and/or Sigma were ignored as cSigma was given")
    }
    p <- ncol(cSigma)
    ## Eigenvectors of the original cov matrix, taken from cSigma
    Uv.org <- svd(cSigma, nu = 0, nv = nv)$v
    Ur.org <- svd(cov2cor(crossprod(cSigma)), nu = 0, nv = nv)$v
    ## Objects to store simulation results
    ansv.L <- matrix(numeric(b * p), p, b)
    ansr.L <- matrix(numeric(b * p), p, b)
    ansv.VE <- numeric(b)
    ansv.VR <- numeric(b)
    ansr.VE <- numeric(b)
    ansr.VR <- numeric(b)
    ansv.U <- array(numeric(b * p * nv), dim = c(p, nv, b))
    ansr.U <- array(numeric(b * p * nv), dim = c(p, nv, b))
    dimnames(ansv.L)[[1]] <- paste0("L", seq_len(p))
    dimnames(ansr.L)[[1]] <- paste0("L", seq_len(p))
    dimnames(ansv.U)[[1]] <- paste0("v", seq_len(p))
    dimnames(ansr.U)[[1]] <- paste0("v", seq_len(p))
    i <- 1L
    ## Simulation runs; a while loop and tryCatch() are used
    ## as svd() (in VE()) occationally returns an error.
    while(i <= b) {
        Xi <- rmvn(cSigma = cSigma, N = N, sqrt_method = sqrt_method)
        ansv <- tryCatch(VE(Xi, center = center, scale. = FALSE, nv = nv,
                            sub = sub, divisor = divisor, m = m,
                            drop_0 = drop_0, tol = tol),
                         error = function(e) list(-1, numeric(p)))
        ansr <- tryCatch(VE(Xi, center = center, scale. = TRUE, nv = nv,
                            sub = sub, divisor = divisor, m = m,
                            drop_0 = drop_0, tol = tol),
                         error = function(e) list(-1, numeric(p)))
        if(ansv[[1]] < 0 || ansr[[1]] < 0) next
        ansv.VE[i] <- ansv$VE
        ansv.VR[i] <- ansv$VR
        ansr.VE[i] <- ansr$VE
        ansr.VR[i] <- ansr$VR
        ansv.L[, i] <- ansv$L
        ansr.L[, i] <- ansr$L
        ansv.U[, , i] <- ansv$U
        ansr.U[, , i] <- ansr$U
        i <- i + 1L
    }
    ans <- list(VEv = ansv.VE, VRv = ansv.VR, Lv = ansv.L, VEr = ansr.VE,
                VRr = ansr.VR, Lr = ansr.L)
    if(nv > 0) {
        ## Simulated eigenvectors are aligned with the original eigenvectors
        ## so that the inner products are positive.
        for(i in seq_len(nv)) {
            Indv <- crossprod(Uv.org[, i], ansv.U[, i, ]) < 0
            Indr <- crossprod(Ur.org[, i], ansr.U[, i, ]) < 0
            ansv.U[, i, Indv] <- -ansv.U[, i, Indv]
            ansr.U[, i, Indr] <- -ansr.U[, i, Indr]
        }
        ans <- c(ans, list(Uv = ansv.U, Ur = ansr.U,
                           Uv.org = Uv.org, Ur.org = Ur.org))
    }
    invisible(ans)
}

##### rmvn #####
#' Generate multivariate normal variates
#'
#' \code{rmvn()} is an internal function to generate
#' multivariate normal variates.
#'
#' @rdname simulateVE
#'
#' @param N
#'   Sample size in each iteration or each run of \code{rmvn()}.
#' @param p
#'   Number of variables; ignored when \code{Sigma} or \code{cSigma} is provided.
#' @param s2
#'   Population variance used for the spherical condition;
#'   ignored when \code{Sigma} or \code{cSigma} is provided.
#' @param Sigma
#'   Population covariance matrix, assumed validly constructed;
#'   by default a spherical covariance is used;
#'   ignored when \code{cSigma} is provided.
#' @param cSigma
#'   Cholesky factor of \code{Sigma}; this can be specified instead of \code{Sigma}.
#'   Potentially useful when multiple simulations are run for the same \code{Sigma}
#'   (although this will not substantially improve speed for a single call of
#'   \code{simulateVE(Sigma = ...)}, where \code{sqrtfun(Sigma)} is called only once).
#' @param mean
#'   Population mean vector; default is a \eqn{p} vector of 0's.
#' @param sqrt_method
#'   Method for matrix square root, which defines the internal function
#'   \code{sqrtfun()}. Choose one from the following:
#'   \code{"chol"} or \code{"default"} for \code{chol()},
#'   \code{"chol_piv"} or \code{"pivot"} for \code{chol_piv()},
#'   \code{"chol_qr"} or \code{"qr"} for \code{chol_qr()},
#'   \code{"matsqrt"} or \code{"sqrt"} for \code{matsqrt()}.
#'   See \link{sqrt_methods} for details. Ignored when \code{cSigma} is provided.
#'
#' @return
#'   \code{rmvn()} returns a \eqn{N * p} matrix with
#'   the specified population mean and covariance.
#'
#' @importFrom stats rnorm
#'
# #' @export
#'
rmvn <- function(N, p = 2L, s2 = 1, Sigma = s2 * diag(p), cSigma = sqrtfun(Sigma),
                 mean = rep_len(0, p), sqrt_method = "chol") {
    if(missing(cSigma)) {
        sqrt_method <- match.arg(sqrt_method, c("default", "chol", "chol_piv",
                                   "chol_qr", "matsqrt", "pivot", "qr", "sqrt"))
        sqrtfun <- switch(sqrt_method, matsqrt = matsqrt, sqrt = matsqrt,
                          chol_qr = chol_qr, qr = chol_qr,
                          pivot = chol_piv, chol_piv = chol_piv, chol)
        cSigma <- sqrtfun(Sigma)
    }
    p <- nrow(cSigma)
    X <- matrix(rnorm(p * N), N, p)
    X <- X %*% cSigma
    # if (missing(mean)) mean <- rep_len(0, p)
    X <- sweep(X, 2, mean, "+")
    return(X)
}
