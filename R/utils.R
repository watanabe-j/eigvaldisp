#### sqrt_methods (dummy) #####
#' Matrix factorization functions
#'
#' Assuming \code{A} is a symmetric, non-negative definite matrix:
#'
#' Cholesky factors are useful in simulating multivariate normal variates,
#' but \code{chol(A)} by default returns an error unless \code{A} is strictly
#' positive definite (i.e., singular matrices are not allowed).
#' One way to avoid this is to use \code{chol(A, pivot = TRUE)} and arrange
#' the columns appropriately (see the documentation of \link[base]{chol}).
#' \code{chol_piv(A)} is a utility to conduct this.
#'
#' Unexpectedly, however, this can still fail in certain situations,
#' e.g., when more than one eigenvalues are 0.
#' This is despite that a Cholesky factorization can be defined for any
#' nonnegative definite matrix, as can be confirmed via QR factorization,
#' (e.g., Schott, 2017; although this factor may not be unique).
#' \code{chol_piv()} returns a warning in this case.
#' \code{chol_qr()} handles this factorization via a naive implementation;
#' first taking a matrix square root via eigendecomposition (\code{matsqrt(A)})
#' and then conducting QR factorization of this (with pivoting).
#' This allows for a correct Cholesky factorization of
#' any positive semidefinite matrix, although for most simulation purposes
#' the matrix square root by \code{matsqrt()} usually suffices.
#'
#' @name sqrt_methods
#' @param A
#'   Numeric matrix, assumed to be symmetric.
#' @param method
#'   Either \code{"svd"} (default) or \code{"eigen"} to specify
#'   the function to be used for eigendecomposition. Results should be
#'   identical, but svd would be faster in most environments.
#' @references
#'   Schott, J. R. (2017). *Matrix Analysis for Statistics*, 3rd edition.
#'     John Wiley & Sons, Hoboken, New Jersey.
#' @seealso \code{\link[base]{chol}}, \code{\link[base]{qr}}
#' @examples
#' (A <- diag(c(2, 1, 1, 0)))
#' \dontrun{chol(A)} # This returns an error because of singularity
#' cA <- eigvaldisp:::chol_piv(A)
#' all.equal(A, crossprod(cA))
#' # TRUE, as expected
#'
#' B <- matrix(1, 4, 4)
#' B[1:2, 3:4] <- B[3:4, 1:2] <- 0
#' print(B)
#' cB1 <- eigvaldisp:::chol_piv(B)
#' all.equal(B, crossprod(cB1))
#' # not TRUE! (though perhaps environment-dependent)
#' crossprod(cB1)
#'
#' cB2 <- eigvaldisp:::chol_qr(B)
#' all.equal(B, crossprod(cB2))
#' # TRUE, as it should be
#' crossprod(cB2)
NULL

##### matsqrt #####
#' Matrix square root
#'
#' \code{matsqrt(A)} returns a matrix square root of \code{A}
#' with spectral or singular value decomposition.
#'
#' @rdname sqrt_methods
#' @return
#'   \code{matsqrt(A)}: Symmetric matrix which is a square root of \code{A}.
matsqrt <- function(A, method = c("svd", "eigen")) {
    method <- match.arg(method)
    if(method == "svd") {
       svdA <- svd(A, nu = 0)
       tv <- t(svdA$v)
       d <- svdA$d
    } else {
        eigenA <- eigen(A)
        tv <- t(eigenA$vectors)
        d <- eigenA$values
    }
    sqrA <- crossprod(tv * sqrt(pmax(d, 0)), tv)
    return(sqrA)
}

#### chol_piv #####
#' Cholesky factorization with pivoting
#'
#' \code{chol_piv(A)} conducts Cholesky factorization of \code{A}
#' with "pivoting".
#' \code{crossprod(chol_piv(A))} will be equal to \code{A},
#' unless multiple zero eigenvalues exist
#' (in which case the results can be spurious; a warning results).
#' @rdname sqrt_methods
#' @return
#'   \code{chol_piv(A)}/\code{chol_qr(A)}: Upper triangular matrix
#'   which is a Cholesky factor of \code{A}.
chol_piv <- function(A) {
    cA <- suppressWarnings(chol(A, pivot = TRUE))
    if(attr(cA, "rank") < ncol(A) - 1) {
        warning("The rank of matrix is ", attr(cA, "rank"),
        ", which is smaller than its dimension by 2 or more.\n  ",
        "chol_piv() may fail in this case. ",
        "Better to use chol_qr() or matsqrt().")
    }
    pivorder <- order(attr(cA, "pivot"))
    cA <- cA[, pivorder]
    return(cA)
}

##### chol_qr #####
#' Cholesky factorization with QR decomposition
#'
#' \code{chol_qr(A)} conducts Cholesky factorization of \code{A}
#' with QR factorization.
#' \code{crossprod(chol_qr(A))} is equal to \code{A},
#' even with multiple zero eigenvalues.
#' This is more reliable, but slower than \code{chol_piv()}.
#' \code{matsqrt()}, which is called internally, usually suffices.
#'
#' @rdname sqrt_methods
# #' @return
# #'   \code{chol_qr(A)}: Upper triangular matrix
# #'   which is the Cholesky factor of \code{A}.
chol_qr <- function(A, method = c("svd", "eigen")) {
    method <- match.arg(method)
    sqrA <- matsqrt(A, method = method)
    QRsqrA <- qr(sqrA)
    cA <- qr.R(QRsqrA)
    pivorder <- order(QRsqrA$pivot)
    cA <- cA[, pivorder]
    return(cA)
}

##### chol_qrsvd #####
#' Cholesky factorization with QR decomposition and SVD
#'
#' \code{chol_qrsvd(A)} is essentially \code{chol_qr(A, method = "svd")}.
#'
#' @rdname sqrt_methods
chol_qrsvd <- function(A) {
    svdA <- svd(A, nu = 0)
    tv <- t(svdA$v)
    sqrA <- crossprod(tv * sqrt(pmax(svdA$d, 0)), tv)
    QRsqrA <- qr(sqrA)
    cA <- qr.R(QRsqrA)
    pivorder <- order(QRsqrA$pivot)
    cA <- cA[, pivorder]
    return(cA)
}


##### cetering #####
#' Utilities for centering/scaling
#' @name centering
#' @param x Numeric matrix
#' @param center,scale Passed to \code{scale(x, ...)}
NULL

##### mc, scale2 #####
#' Utility for centering
#'
#' \code{mc(x)} is a short for \code{sweep(x, 2, colMeans(x))}.
#' @rdname centering
mc <- function(x) {
    sweep(x, 2, colMeans(x))
}

##### scale2 #####
#' Utility for scaling
#'
#' \code{scale2(x)} is essentially \code{scale(x)},
#' but returns 0 rather than NaN for constant columns.
#' @rdname centering
scale2 <- function(x, center = TRUE, scale = FALSE) {
    ans <- scale(x, center = center, scale = scale)
    ans[, attr(ans, "scaled:scale") == 0] <- 0
    return(ans)
}


##### digit #####
#' Obtain number of digits
#'
#' Internal utility function to obtain the number of digit of an integer.
#' For internal use in \code{AVar.VRr_pfx()} functions (and others).
#'
#' @param x
#'   The object whose number of digits is to be obtained.
#'   Length can be >1. Integer or integer-like numeric are assumed.
#'
#' @return The number of digits, as an integer vector
#'
digit <- function(x) {
    as.integer(floor(log(x, 10))) + 1L
}

##### divInd #####
#' Divide index vector
#'
#' Utility function to divide a vector of index to a list.
#' To be used internally in \code{AVar.VRr_pfd()}.
#'
#' The argument \code{Max} defines the maximum of \code{length(bd[[i]])} in
#' \code{AVar.VRr_pfd}, which roughly defines the amount of RAM required
#' (see description of that function).
# #' The number of steps (s) is defined as the maximum integer that satisfies
# #' the necessary condition: \code{(p + p - s) * (s + 1) / 2 < Max}
# #' Alternatively, Length could be specified to define the length of the list.
#'
#' @param b
#'   Vector to be divided into a list
#' @param Max
#'   Maximum acceptable length of vector that results from the divided list
# #' @param Length
# #'   Length of the list generated
#'
#' @return A list
divInd <- function(b, Max = 2e6) { # , Length = 2) {
    p <- length(b)
    # if(missing(Max) && !missing(Length)) Max <- ceiling(nc / (Length - 1) - 1)
    ans <- list()
    while(p > 0) {
        if((2 * p - 1) ^ 2 + 8 * (p - Max) < 0) {
            s <- p - 1
        } else {
            s <- pmax(
                floor((2 * p - 1 - sqrt((2 * p - 1) ^ 2 + 8 * (p - Max))) / 2),
                0)
        }
        ans <- c(ans, list(b[1:(s + 1)]))
        b <- b[-(1:(s + 1))]
        p <- length(b)
    }
    return(ans)
}

##### hgf #####
#' Hypergeometric function wrapper
#'
#' Wrapper for hypergeometric function with two upper and one lower arguments
#' (or the classic hypergeometric function).
#'
#' This function utilizes the function \code{genhypergeo} from
#' the package \code{hypergeo}; when it appears to have failed,
#' then the function \code{hypergeo} is used.
#'
#' Elements of \code{x} are assumed to be real numbers within \eqn{[-1, 1]}.
#' Complex arguments/returns are not assumed. No check is done.
#'
#' Tries to reduce computational time by dropping duplicated values from
#' \code{x} while trying to retain its original \code{names} and \code{dim}.
#'
#' This function primarily uses \code{hypergeo::genhypergeo()}, as this seems
#' to be more numerically stable than \code{hypergeo::hypergeo()} in some
#' conditions tested in development.
#' For the present purpose, \code{gsl::hyperg_2F1} may suffice,
#' but this has not been implemented.
#'
#' @param a1,a2 Parameters in the numerator
#' @param b1    Parameter in the denominator
#' @param x     Argument of the hypergeometric function
#' @param tol,maxiter  Passed to \code{genhypergeo} (and \code{hypergeo})
#'
#' @return A numeric vector corresponding to \code{x}
#'
#' @seealso \code{\link[hypergeo]{genhypergeo}}, 
#'   \code{\link[hypergeo]{hypergeo}}
#'
hgf <- function(a1, a2, b1, x, tol = 0, maxiter = 2000) {
    # require(hypergeo)
    if(any(length(a1) > 1, length(a2) > 1, length(b1) > 1)) {
        stop("Only scalars are acceptable for a1, a2, and b1")
    }
    ox <- x
    dx <- dim(ox)
    nx <- names(ox)
    dim(ox) <- NULL
    names(ox) <- ox
    x <- unique(ox)
    names(x) <- ox[!duplicated(ox)]
    F <- Re(hypergeo::genhypergeo(c(a1, a2), b1, x,
                                  tol = tol, maxiter = maxiter))
    if(any(is.na(F) & !is.na(x))) {
        F <- ifelse(is.na(F), Re(hypergeo::hypergeo(a1, a2, b1,
                                           x, tol = tol, maxiter = maxiter)), F)
    }
    F <- F[names(ox)]
    names(F) <- nx
    dim(F) <- dx
    return(F)
}
