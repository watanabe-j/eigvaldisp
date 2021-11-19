##### Exv.rx #####
#' Moments of correlation coefficients
#'
#' Internal functions to calculate exact and asymptotic moments of
#' sample correlation coefficients \eqn{r} from
#' population correlation coefficients.
#'
#' All these functions except \code{ACov.r1()} requires the arguments
#' \code{n} and \code{x}. It can be \code{length(x) > 1}
#' (at least supeficially vectorized), but it is assumed \code{length(n) = 1}.
#'
#' Covariance between two \eqn{r}'s is a function of
#' at least six population correlation coefficients.
#' \code{ACov.r1()} takes as an argument the entire population correlation
#' matrix \code{R} instead of individual correlation coefficients, and indices
#' for the focal correlation coefficients. Alternatively, the all six relevant
#' population correlation coefficients can directly be specified.
#'
#' The exact expressions are from Soper et al. (1917) or Ghosh (1966).
#' In particular:
#' \itemize{
#'   \item \code{Exv.r3()}: From Soper et al. (1917, eq. xxviii);
#'     that of Ghosh (1966) seems inaccurate for this moment.
#'   \item \code{Exv.r3d()}: Alternative to \code{Exv.r3()} from
#'     Soper et al. (1917, eq. xxvi). Equivalent but slightly slower.
#'   \item \code{Exv.r4()}: From Ghosh (1966, eq. 1).
#'   \item \code{Exv.r4d()}: Alternative to \code{Exv.r4()} from
#'     Soper et al. (1917, eq. xxi). Equivalent but slightly slower.
#' }
#'
#' Some of the exact expressions might be instable when the \code{x} is small.
#'
#' The asymptotic moments are from Ghosh (1966)
#' (note that the symbol \eqn{n} is used for sample size there).
#' *Note*: The validity of \code{ACm3.r1()} may need critical evaluation,
#' as Ghosh's (1966) equation 1 (from which the result is supposedly drawn)
#' seems inaccurate (see Soper et al. \[1917: eq. xxviii\]).
#' The following functions depending on this might be inaccurate:
#' \code{AExv.r3()}, \code{AExv.r4()}, and \code{AVar.r2()}.
#'
#' The asymptotic moment functions takes the argument \code{order.}.
#' This is the exponent in \eqn{O(n^r)}.
#' Higher orders do not always yield better approximations for ordinary
#' values of \eqn{n}. The allowed ranges and "good" choice
#' (empirically confirmed by comparison to exact expressions) are as follows:
#' \itemize{
#'   \item \code{AExv.r1()},  \code{AExv.r2()}: 1--7. 2 is good
#'     except when \eqn{n < 6} where 3 is better.
#'   \item \code{AVar.r1()}: 1--7. 2 is good for \eqn{n} <100 or so;
#'     for larger \eqn{n}, the difference tends to be negligible.
#'   \item \code{ACm3.r1()}: 2--5. 3 is best for \eqn{n} > 15,
#'     while 4 may work better for smaller \eqn{n}.
#'   \item \code{ACm4.r1()}: 2--5. 3 is good for \eqn{n} > 20 especially when
#'     \eqn{x} is small; 4 is better for small \eqn{n} and \eqn{x > 0.5}.
#'   \item \code{AExv.r3()}: 1--5. 2 is good.
#'   \item \code{AExv.r4()}: 1--5. 3 or 4 is good for small \eqn{x},
#'     while 2 is better for \eqn{x > 0.4}.
#'   \item \code{AVar.r2()}: 1--5. For \eqn{n < 100}, 3 or 4 is good for
#'     very small \eqn{x}, while 2 is better for \eqn{x > 0.2}.
#'     Slight overestimation seems to be common.
#' }
#'
#' \code{ACov.r1()} is as in, e.g., Olkin & Siotani (1976), Olkin & Finn (1990).
#' When \code{ind1 == ind2}, it calculates the variance for the coefficient.
#'
#' Multivariate normality is assumed for all moments.
#' Nevertheless, the distribution of \eqn{r} remains the same in many other
#' conditions (e.g., in all elliptically contoured distributions;
#' Anderson, 2003).
#'
#' @name Exv.rx
#'
#' @param n
#'   Degrees of freedom (not sample size), numeric of length 1.
#' @param x
#'   Population correlation coefficient, numeric of length 1 or more.
#'   For some functions, squared correlation coefficients can be passed
#'   as an argument when \code{do.square = FALSE}.
#' @param tol
#'   Tolerance used to judge wheter x is sufficiently close to 0 or +/- 1
#'   (where numerical instability can happen).
#' @param tol.hg,maxiter
#'   Passed to \code{hgf()} as tol and maxiter there.
#' @param do.square
#'   Whether \code{x} is to be squared before passed to \code{hgf()};
#'   \code{TRUE} by default; set this \code{FALSE}
#'   when \code{x} is already squared.
#' @param order.
#'   Order of asymptotic approximation (exponent \eqn{r} in \eqn{O(n^r)};
#'   for asymptotic expressions only).
#' @param return_terms
#'   When TRUE, returns individual terms of different orders
#'   (for asymptotic expressions only).
#'   Retained mainly for debugging purposes, but could be used to
#'   track the order in derived statistics.
#' @param Rho
#'   Correlation matrix; assumed to be validly constructed.
#' @param ind1,ind2
#'   Indices for the two focal coefficients (vectors of length 2).
#'   Superseded by \code{i, j, k, l}.
#' @param i,j,k,l
#'   Alternative indices to specify the focal coefficients.
#'   Superseded by \code{rij, rik, ril, rjk, rjl, rkl}.
#' @param rij,rik,ril,rjk,rjl,rkl
#'   The relevant population correlation coefficients.
#' @param ...
#'   Additional arguments in these functions are silently ignored.
#'
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#'  Analysis*, 3rd edition. John Wiley & Sons, Hoboken, New Jersey.
#'
#' Ghosh, B. K. (1966). Expansions for the moments of the distribution
#'  of correlation coefficients. *Biometrika* **53**, 258--262.
#'  doi:[10.2307/2334076](https://doi.org/10.2307/2334076).
#'
#' Olkin, I. & Finn, J. D. (1990). Testing correlated correlations.
#'  *Psychological Bulletin* **108**, 330--333.
#'  doi:[10.1037/0033-2909.108.2.330](https://doi.org/10.1037/0033-2909.108.2.330).
#'
#' Olkin, I. & Siotani, M. (1976). Asymptotic distribution of functions of
#'  a correlation matrix. Pp. 235--251 *in* Editorial Committee for Publication
#'  of Essays in Probability and Statistics, eds. *Essays in Probability and
#'  Statistics in Honor of Professor Junjiro Ogawa*. Shinko Tsusho, Tokyo.
#'
#' Soper, H. E., Young, A. W., Cave, B. M., Lee, A. & Pearson, K. (1917).
#'  On the distribution of the correlation coefficients in small samples.
#'  Appendix II to the papers of ``Student'' and R. A. Fisher.
#'  *Biometrika* **11**, 328--413.
#'  doi:[10.1093/biomet/11.4.328](https://doi.org/10.1093/biomet/11.4.328).
#'
#' @seealso
#' \link{Exv.VXX}, \link{AVar.VRR_xx} for outer functions.
#'
#' @examples
#' # Trivial examples
#' n <- 10    # Which typically corresponds to sample size of 11
#' eigvaldisp:::Exv.r1(n, 0)
#' eigvaldisp:::Exv.r1(n, 0.5) # Underestimate
#' eigvaldisp:::Exv.r2(n, 0)   # This is 1/n for x = 0
#'
#' eigvaldisp:::AExv.r1(n, 0)
#' eigvaldisp:::AExv.r1(n, 0.5)
#' eigvaldisp:::AExv.r2(n, 0)  # Fairly close approximations of the above
#'
#' # Comparison between exact vs. asymptotic expressions
#' n <- 3                 # Extremely small n
#' x <- seq(0, 1, 0.01)
#' cols <- rainbow(7)
#'
#' # Exv.r1() vs AExv.r1() with different orders
#' plot(x, sapply(x, eigvaldisp:::Exv.r1, n = n), type = "l", lwd = 2,
#'      xlim = c(0, 1), ylim = c(0, 1))
#' for(i in 1:7) {
#'     lines(x, sapply(x, eigvaldisp:::AExv.r1, n = n, order = i),
#'           col = cols[i])
#' }
#' legend("topleft", title = "Order", legend = c(1:7, "Exact"),
#'        col = c(cols, "black"), lty = 1)
#'
#' # Exv.r2() vs AExv.r2() with different orders
#' plot(x, sapply(x, eigvaldisp:::Exv.r2, n = n), type = "l", lwd = 2,
#'      xlim = c(0, 1), ylim = c(0, 1))
#' for(i in 1:7) {
#'     lines(x, sapply(x, eigvaldisp:::AExv.r2, n = n, order = i),
#'           col = cols[i])
#' }
#' legend("topleft", title = "Order", legend = c(1:7, "Exact"),
#'        col = c(cols, "black"), lty = 1)
#'
NULL

##### Exv.r1 #####
#' Expectation of correlation coefficient
#'
#' \code{Exv.r1()}: exact expectation of \eqn{r}.
#'
#' @rdname Exv.rx
#'
Exv.r1 <- function(n, x, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    C <- 2 / n * (gamma((n + 1) / 2) / gamma(n / 2)) ^ 2
    ifelse(abs(x) > 1 - tol, x, # When x ~ 1, the below may not converge
        ifelse(abs(x) < tol,
               0, # When x ~ 0
               C * x * hgf(1 / 2, 1 / 2, (n + 2) / 2, x ^ 2,
                           tol = tol.hg, maxiter = maxiter)))
}

##### Exv.r2 #####
#' Expectation of correlation coefficient squared
#'
#' \code{Exv.r2()}: exact expectation of \eqn{r^2}.
#'
#' @rdname Exv.rx
#'
Exv.r2 <- function(n, x, do.square = TRUE, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    if(do.square) x <- x ^ 2
    ifelse(abs(x) > 1 - tol, x,
        ifelse(abs(x) < tol,
               1 / n, # Same as below when x == 0
               {F <- hgf(1, 1, (n + 2) / 2, x,
                         tol = tol.hg, maxiter = maxiter);
                1 - (n - 1) * (1 - x) / n * F}))
}

##### Exv.r3 #####
#' Expectation of correlation coefficient cubed
#'
#' \code{Exv.r3()}: exact expectation of \eqn{r^3}.
#'
#' @rdname Exv.rx
#'
Exv.r3 <- function(n, x, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    mu1_n0 <- Exv.r1(n, x, tol = tol, tol.hg = tol.hg, maxiter = maxiter)
    mu1_n2 <- Exv.r1(n + 2, x, tol = tol, tol.hg = tol.hg, maxiter = maxiter)
    mu1_n0 - n * (n - 1) * (mu1_n2 - mu1_n0)
}

##### Exv.r3d #####
#' Expectation of correlation coefficient cubed, alternative
#'
#' \code{Exv.r3d()}: alternative to \code{Exv.r3()}.
#'
#' @rdname Exv.rx
#'
Exv.r3d <- function(n, x, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    x2 <- x ^ 2
    C <- 2 / n * (gamma((n + 1) / 2) / gamma(n / 2)) ^ 2
    ifelse(abs(x) > 1 - tol, x, # When x ~ 1, the below may not converge
       ifelse(abs(x) < tol,
              0, # When x ~ 0
              {F1 <- hgf(1 / 2, 1 / 2, (n + 2) / 2, x2,
                         tol = tol.hg, maxiter = maxiter);
               F2 <- hgf(3 / 2, 3 / 2, (n + 4) / 2, x2,
                         tol = tol.hg, maxiter = maxiter);
               C * x * (F1 - (n - 1) * (1 - x2) / (n + 2) * F2)}))
}

##### Exv.r4 #####
#' Expectation of correlation coefficient powered to 4
#'
#' \code{Exv.r4()}: exact expectation of \eqn{r^4}.
#'
#' @rdname Exv.rx
#'
Exv.r4 <- function(n, x, do.square = TRUE, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    if(do.square) x <- x ^ 2
    ifelse(abs(x) > 1 - tol, x,
        ifelse(abs(x) < tol,
            3 / n / (n + 2), # Same as below when x == 0
            {F <- hgf(1, 1, (n + 2) / 2, x, tol = tol.hg, maxiter = maxiter);
             Fd <- hgf(1, 2, (n + 4) / 2, x, tol = tol.hg, maxiter = maxiter);
             1 + (n - 1) * (n - 3) * (1 - x) / 2 / n * F -
             (n + 1) * (n - 1) * (1 - x) / 2 / (n + 2) * Fd}))
}

##### Exv.r4d #####
#' Expectation of correlation coefficient powered to 4, alternative
#'
#' \code{Exv.r4d()}: alternative to \code{Exv.r4()}.
#'
#' @rdname Exv.rx
#'
Exv.r4d <- function(n, x, do.square = TRUE, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    if(do.square) x <- x ^ 2
    ifelse(abs(x) > 1 - tol, x,
        ifelse(abs(x) < tol,
            3 / n / (n + 2), # Same as below when x == 0
            {F <- hgf(1, 1, (n + 2) / 2, x, tol = tol.hg, maxiter = maxiter);
             Fd <- hgf(2, 2, (n + 4) / 2, x, tol = tol.hg, maxiter = maxiter);
             1 - 2 * (n - 1) * (1 - x) / n * F +
             (n + 1) * (n - 1) * (1 - x)^2 / n / (n + 2) * Fd}))
}

##### Var.r1 #####
#' Variance of correlation coefficient
#'
#' \code{Var.r1()}: exact variance of \eqn{r}.
#'
#' @rdname Exv.rx
#'
Var.r1 <- function(n, x, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    mu2_r1 <- Exv.r2(n, x, do.square = TRUE, tol = tol,
                     tol.hg = tol.hg, maxiter = maxiter)
    mu1_r1 <- Exv.r1(n, x, tol = tol, tol.hg = tol.hg, maxiter = maxiter)
    mu2_r1 - mu1_r1 ^ 2
}

##### Var.r2 #####
#' Variance of correlation coefficient squared
#'
#' \code{Var.r2()}: exact variance of \eqn{r^2}.
#'
#' @rdname Exv.rx
#'
Var.r2 <- function(n, x, do.square = TRUE, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    if(do.square) x <- x ^ 2
    ifelse(abs(x) > 1 - tol, tol, # When x ~ 1, the below may not converge
        ifelse(abs(x) < tol,
            2 * (n - 1) / (n ^ 2 * (n + 2)), # When x ~ 0
            {F <- hgf(1, 1, (n + 2) / 2, x, tol = tol.hg, maxiter = maxiter);
             Fd <- hgf(1, 2, (n + 4) / 2, x, tol = tol.hg, maxiter = maxiter);
             (n - 1) * (n + 1) * (1 - x) / 2 / n *
             (F - n / (n + 2) * Fd - 2 * (n - 1) * (1 - x)
                                    / n / (n + 1) * F ^ 2)}))
}

##### Cm3.r1 #####
#' Central third moment of correlation coefficient
#'
#' \code{Cm3.r1()}: exact central third moment of \eqn{r}.
#'
#' @rdname Exv.rx
#'
Cm3.r1 <- function(n, x, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    mu3_r1 <- Exv.r3(n = n, x = x, tol = tol, tol.hg = tol.hg, maxiter = maxiter)
    mu2_r1 <- Exv.r2(n, x, do.square = TRUE, tol = tol,
                     tol.hg = tol.hg, maxiter = maxiter)
    mu1_r1 <- Exv.r1(n, x, tol = tol, tol.hg = tol.hg, maxiter = maxiter)
    mu3_r1 - 3 * mu2_r1 * mu1_r1 + 2 * mu1_r1 ^ 3
}

##### Cm4.r1 #####
#' Central fourth moment of correlation coefficient
#'
#' \code{Cm4.r1()}: exact central fourth moment of \eqn{r}.
#'
#' @rdname Exv.rx
#'
Cm4.r1 <- function(n, x, tol = .Machine$double.eps * 100,
                   tol.hg = 0, maxiter = 2000, ...) {
    mu_r4 <- Exv.r4(n = n, x = x, tol = tol, tol.hg = tol.hg, maxiter = maxiter)
    mu_r3 <- Exv.r3(n = n, x = x, tol = tol, tol.hg = tol.hg, maxiter = maxiter)
    mu_r2 <- Exv.r2(n = n, x = x, tol = tol, tol.hg = tol.hg, maxiter = maxiter)
    mu_r1 <- Exv.r1(n = n, x = x, tol = tol, tol.hg = tol.hg, maxiter = maxiter)
    mu_r4 - 4 * mu_r3 * mu_r1 + 6 * mu_r2 * mu_r1 ^ 2 - 3 * mu_r1 ^ 4
}

##### AExv.r1 #####
#' Asymptotic expectation of correlation coefficient
#'
#' \code{AExv.r1()}: asymptotic expectation of \eqn{r}.
#'
#' @rdname Exv.rx
#'
AExv.r1 <- function(n, x, order. = 2, return_terms = FALSE, ...) {
    if(order. > 7 || order. < 1) {
        order.old <- order.
        order. <- pmin(pmax(order., 1), 7)
        warning("order. = ", order.old, " is not allowed. order. = ",
                order., " was used.")
    }
    m <- n + 5
    lx <- length(x)
    rs <- sapply(x, function(z) z ^ (2 * seq.int(0, 6)))
    r2 <- rs[2, ]
    C <- matrix(c(
                3,        1,        0,        0,        0,        0,        0,
              121,       70,       25,        0,        0,        0,        0,
             6479,     4923,     2925,     1225,        0,        0,        0,
            86341,    77260,    58270,    38220,    19845,        0,        0,
          2290597,  2281659,  1972050,  1575350,  1157625,   800415,        0,
         30242033, 32443250, 30545375, 27171900, 23147775, 18523890, 19324305),
         nrow = 6, ncol = 7, byrow = TRUE)
    ms <- (m ^ seq.int(-1, -6) / c(4/3, 8, 64, 128, 512, 1024)) * 3
    tms <- ms[seq_len(order. - 1)] * C[seq_len(order. - 1), ] %*% rs
    tms <- rbind(1, tms)
    tms <- t(- x * (1 - r2) / 2 / m * t(tms))
    tms <- rbind(x, tms)
    rownames(tms) <- seq.int(0, order.)
    ans <- colSums(tms)
    if(return_terms){
        return(list(exp.value = ans, terms = tms))
    } else {
        return(ans)
    }
}

##### AVar.r1 #####
#' Asymptotic variance of correlation coefficient
#'
#' \code{AVar.r1()}: asymptotic variance of \eqn{r}.
#'
#' @rdname Exv.rx
#'
AVar.r1 <- function(n, x, order. = 2, return_terms = FALSE, ...) {
    if(order. > 7 || order. < 1) {
        order.old <- order.
        order. <- pmin(pmax(order., 1), 7)
        warning("order. = ", order.old, " is not allowed. order. = ",
                order., " was used.")
    }
    m <- n + 5
    lx <- length(x)
    rs <- sapply(x, function(z) z ^ (2 * seq.int(0, 6)))
    r2 <- rs[2, ]
    C <- matrix(c(
               14,       11,        0,        0,        0,        0,        0,
               98,      130,       75,        0,        0,        0,        0,
             2744,     4645,     4422,     2565,        0,        0,        0,
            19208,    37165,    44499,    40299,    26685,        0,        0,
           268912,   561743,   762882,   838182,   784770,   657135,        0,
          1882384,  4105472,  6005139,  7313100,  7924830,  7692660,  9365895),
          nrow = 6, ncol = 7, byrow = TRUE)
    ms <- (m ^ seq.int(-1, -6) / c(2, 2, 8, 8, 16, 16))
    tms <- ms[seq_len(order. - 1)] * C[seq_len(order. - 1), ] %*% rs
    tms <- rbind(1, tms)
    tms <- t((1 - r2) ^ 2 / m * t(tms))
    rownames(tms) <- seq.int(1, order.)
    ans <- colSums(tms)
    if(return_terms){
        return(list(exp.value = ans, terms = tms))
    } else {
        return(ans)
    }
}

##### ACm3.r1 #####
#' Asymptotic central third moment of correlation coefficient
#'
#' \code{ACm3.r1()}: asymptotic central third moment of \eqn{r}.
#'
#' @rdname Exv.rx
#'
ACm3.r1 <- function(n, x, order. = 3, return_terms = FALSE, ...) {
    if(order. > 5 || order. < 2) {
        order.old <- order.
        order. <- pmin(pmax(order., 2), 5)
        warning("order. = ", order.old, " is not allowed. order. = ",
                order., " was used.")
    }
    m <- n + 5
    lx <- length(x)
    rs <- sapply(x, function(z) z ^ (2 * seq.int(0, 3)))
    r2 <- rs[2, ]
    C <- matrix(c(
               69,    88,     0,     0,
              797,  1691,  1560,     0,
            12325, 33147, 48099, 44109),
          nrow = 3, ncol = 4, byrow = TRUE)
    ms <- (m ^ seq.int(-1, -3) / c(1, 4 / 3, 8 / 3))
    tms <- ms[seq_len(order. - 2)] * C[seq_len(order. - 2), ] %*% rs
    tms <- rbind(6, tms)
    tms <- t(- x * (1 - r2) ^ 3 / m ^ 2 * t(tms))
    rownames(tms) <- seq.int(2, order.)
    ans <- colSums(tms)
    if(return_terms){
        return(list(exp.value = ans, terms = tms))
    } else {
        return(ans)
    }
}

##### ACm4.r1 #####
#' Asymptotic central fourth moment of correlation coefficient
#'
#' \code{ACm4.r1()}: asymptotic central fourth moment of \eqn{r}
#'
#' @rdname Exv.rx
#'
ACm4.r1 <- function(n, x, order. = 3, return_terms = FALSE, ...) {
    if(order. > 5 || order. < 2) {
        order.old <- order.
        order. <- pmin(pmax(order., 2), 5)
        warning("order. = ", order.old, " is not allowed. order. = ",
                order., " was used.")
    }
    m <- n + 5
    lx <- length(x)
    rs <- sapply(x, function(z) z ^ (2 * seq.int(0, 3)))
    r2 <- rs[2, ]
    C <- matrix(c(
              12,    35,     0,     0,
             436,  2028,  3025,     0,
            3552, 20009, 46462, 59751),
          nrow = 3, ncol = 4, byrow = TRUE)
    ms <- (m ^ seq.int(-1, -3) / c(1, 4, 4))
    tms <- ms[seq_len(order. - 2)] * C[seq_len(order. - 2), ] %*% rs
    tms <- rbind(1, tms)
    tms <- t(3 * (1 - r2) ^ 4 / m ^ 2 * t(tms))
    rownames(tms) <- seq.int(2, order.)
    ans <- colSums(tms)
    if(return_terms){
        return(list(exp.value = ans, terms = tms))
    } else {
        return(ans)
    }
}

##### AExv.r2 #####
#' Asymptotic expectation of correlation coefficient squared
#'
#' \code{AExv.r2()}: asymptotic expectation of \eqn{r^2}
#'
#' @rdname Exv.rx
#'
AExv.r2 <- function(n, x, order. = 2, return_terms = FALSE, ...) {
    if(order. > 7 || order. < 1) {
        order.old <- order.
        order. <- pmin(pmax(order., 1), 7)
        warning("order. = ", order.old, " is not allowed. order. = ",
                order., " was used.")
    }
    m <- n + 5
    lx <- length(x)
    rs <- sapply(x, function(z) z ^ (2 * seq.int(0, 7)))
    r2 <- rs[2, ]
    C <- matrix(c(
               1,    -2,     0,      0,      0,      0,      0,       0,
               7,    -8,    -8,      0,      0,      0,      0,       0,
              49,   -26,   -56,    -48,      0,      0,      0,       0,
             343,   -32,  -272,   -384,   -384,      0,      0,       0,
            2401,   526,  -944,  -2016,  -2688,  -3840,      0,       0,
           16807,  7432,  -728,  -7680, -13440, -15360, -46080,       0,
          117649, 70774, 27544, -12048, -48000, -88320,  46080, -645120),
          nrow = 7, ncol = 8, byrow = TRUE)
    ms <- m ^ seq.int(-0, -6)
    tms <- ms[seq_len(order.)] * C[seq_len(order.), ] %*% rs
    tms <- t((1 - r2) / m * t(tms))
    tms <- rbind(r2, tms)
    rownames(tms) <- seq.int(0, order.)
    ans <- colSums(tms)
    if(return_terms){
        return(list(exp.value = ans, terms = tms))
    } else {
        return(ans)
    }
}


##### AExv.r3 #####
#' Asymptotic expectation of correlation coefficient cubed
#'
#' \code{AExv.r3()}: asymptotic expectation of \eqn{r^3}
#'
#' @rdname Exv.rx
#'
AExv.r3 <- function(n, x, order. = 2, return_terms = FALSE, ...) {
    if(order. > 5 || order. < 1) {
        order.old <- order.
        order. <- pmin(pmax(order., 1), 5)
        warning("order. = ", order.old, " is not allowed. order. = ",
                order., " was used.")
    }
    m <- n + 5
    lx <- length(x)
    rs <- sapply(x, function(z) z ^ (2 * seq.int(0, 5)))
    r2 <- rs[2, ]
    C <- matrix(c(
              2,     -3,      0,      0,      0,        0,
             12,      1,    -25,      0,      0,        0,
            306,    371,   -100,  -1225,      0,        0,
           3112,   7011,   6175,   1225, -33075,        0,
          54486, 220811, 309080, 302330, 674730, -2401245),
          nrow = 5, ncol = 6, byrow = TRUE)
    ms <- m ^ seq.int(-0, -4) / c(2, 8/3, 16, 128/3, 256)
    tms <- ms[seq_len(order.)] * C[seq_len(order.), ] %*% rs
    tms <- t(3 * x * (1 - r2) / m * t(tms))
    tms <- rbind(x ^ 3, tms)
    rownames(tms) <- seq.int(0, order.)
    ans <- colSums(tms)
    if(return_terms){
        return(list(exp.value = ans, terms = tms))
    } else {
        return(ans)
    }
}

##### AExv.r4 #####
#' Asymptotic expectation of correlation coefficient powered to 4
#'
#' \code{AExv.r4()}: asymptotic expectation of \eqn{r^4}
#'
#' @rdname Exv.rx
#'
AExv.r4 <- function(n, x, order. = 2, return_terms = FALSE, ...) {
    if(order. > 5 || order. < 1) {
        order.old <- order.
        order. <- pmin(pmax(order., 1), 5)
        warning("order. = ", order.old, " is not allowed. order. = ",
                order., " was used.")
    }
    m <- n + 5
    lx <- length(x)
    rs <- sapply(x, function(z) z ^ (2 * seq.int(0, 6)))
    r2 <- rs[2, ]
    C <- matrix(c(
            0,    3,   -4,    0,     0,     0,      0,
            1,    1,   16,  -24,     0,     0,      0,
            6,   -9,   16,   88,  -128,     0,      0,
          109, -131,  -96,  144,  2688, -3200,      0,
          444, -291, -596, -464, -3200, 24960, -23040),
          nrow = 5, ncol = 7, byrow = TRUE)
    ms <- m ^ seq.int(-0, -4) / c(1/2, 1/3, 1/6, 1/3, 1/6)
    tms <- ms[seq_len(order.)] * C[seq_len(order.), ] %*% rs
    tms <- t((1 - r2) / m * t(tms))
    tms <- rbind(x ^ 4, tms)
    rownames(tms) <- seq.int(0, order.)
    ans <- colSums(tms)
    if(return_terms){
        return(list(exp.value = ans, terms = tms))
    } else {
        return(ans)
    }
}

##### AVar.r2 #####
#' Asymptotic variance of correlation coefficient squared
#'
#' \code{AVar.r2()}: asymptotic variance of \eqn{r^2}
#'
#' @rdname Exv.rx
#'
AVar.r2 <- function(n, x, order. = 2, return_terms = FALSE, ...) {
    if(order. > 5 || order. < 1) {
        order.old <- order.
        order. <- pmin(pmax(order., 1), 5)
        warning("order. = ", order.old, " is not allowed. order. = ",
                order., " was used.")
    }
    m <- n + 5
    lx <- length(x)
    rs <- sapply(x, function(z) z ^ (2 * seq.int(0, 5)))
    r2 <- rs[2, ]
    C <- matrix(c(
            0,    1,     0,     0,     0,     0,
            1,   -2,    26,     0,     0,     0,
           11,  -36,     8,   320,     0,     0,
           45,  -98,  -230,   -64,  2144,     0,
          323, -325, -1736, -2592, -6752, 32064),
          nrow = 5, ncol = 6, byrow = TRUE)
    ms <- m ^ seq.int(-0, -4) / c(1/4, 1/2, 1/2, 1/4, 1/4)
    tms <- ms[seq_len(order.)] * C[seq_len(order.), ] %*% rs
    tms <- t((1 - r2) ^ 2 / m * t(tms))
    rownames(tms) <- seq.int(1, order.)
    ans <- colSums(tms)
    if(return_terms){
        return(list(exp.value = ans, terms = tms))
    } else {
        return(ans)
    }
}

##### ACov.r1 #####
#' Asymptotic covariance of correlation coefficients
#'
#' \code{ACov.r1()}: asymptotic covariance between two \eqn{r}'s.
#'
#' @rdname Exv.rx
#'
ACov.r1 <- function(n, Rho, ind1 = c(1, 2), ind2 = c(1, 2),
                    i = ind1[1], j = ind1[2], k = ind2[1], l = ind2[2],
                    rij = Rho[i, j], rik = Rho[i, k], ril = Rho[i, l],
                    rjk = Rho[j, k], rjl = Rho[j, l], rkl = Rho[k, l], ...) {
    (rij * rkl * (rik^2 + ril^2 + rjk^2 + rjl^2) / 2 + rik * rjl + ril * rjk -
    (rij * rik * ril + rij * rjk * rjl + rik * rjk * rkl + ril * rjl * rkl)) / n
}
