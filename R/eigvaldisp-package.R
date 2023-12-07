#' eigvaldisp: Statistics of Eigenvalue Dispersion Indices
#'
#' This package involves functions for analyzing eigenvalue dispersion
#' indices of covariance and correlation matrices, providing practical
#' implementations for theoretical results of Watanabe (2022).
#'
#' The DESCRIPTION file:
#' \packageDESCRIPTION{eigvaldisp}
#' \packageIndices{eigvaldisp}
#'
#' Run \code{vignettes("eigvaldisp")} for detailed descriptions and examples.
#'
#' @section Author/Maintainer:
#' Junya Watanabe <Junya.Watanabe@uab.cat>
#'
#' @references
#' Watanabe, J. (2022) Statistics of eigenvalue dispersion indices:
#'  quantifying the magnitude of phenotypic integration. *Evolution*,
#'  **76**, 4--28. \doi{10.1111/evo.14382}.
#'
#' @seealso
#'   \code{\link{VE}}: Calculate eigenvalue dispersion indices
#'
#'   \code{\link{Exv.VXX}}: Moments of eigenvalue dispersion indices
#'
#'   \code{\link{AVar.VRR_xx}}: Approximate variance of relative eigenvalue
#'                              variance of correlation matrix (internal functions)
#'
#'   \code{\link{VXXa}}: \dQuote{Bias-corrected} eigenvalue dispersion indices
#'
#'   \code{\link{Exv.VXXa}}: Moments of \dQuote{bias-corrected} eigenvalue dispersion
#'                           indices
#'
#'   \code{\link{simulateVE}}: Simulate eigenvalue dispersion indices
#'
#'   \code{\link{GenCov}}: Generate covariance/correlation matrix with
#'                         known structure
#'
#'   \code{\link{Exv.rx}}: Moments of correlation coefficients (internal functions)
#'
#'   \code{\link{hgf}}, \code{\link{sqrt_methods}}, \code{\link{centering}},
#'   \code{\link{digit}}: Other internal utility functions
#'
#' @examples
#' ## Generate a population covariance matrix with known eigenvalues
#' Lambda <- c(4, 2, 1, 1)
#' (Sigma <- GenCov(evalues = Lambda, evectors = "random"))
#' all.equal(Lambda, eigen(Sigma)$values)  # TRUE
#'
#' ## Calculate eigenvalue dispersion indices of this matrix
#' EDI_pop <- VE(S = Sigma)
#'
#' ## Population eigenvalue variance ("V(Sigma)") and
#' ## relative eigenvalue variance ("Vrel(Sigma)"):
#' EDI_pop$VE
#' EDI_pop$VR
#'
#' ## Simulate a multivariate normal sample
#' N <- 20
#' X <- rmvn(N = N, Sigma = Sigma)
#'
#' ## Calculating eigenvalue dispersion indices from the sample
#' EDI_sam <- VE(X = X)
#'
#' ## Sample eigenvalue variance ("V(S)") and
#' ## relative eigenvalue variance ("Vrel(S)")
#' EDI_sam$VE
#' EDI_sam$VR
#' ## These are typically biased upward
#'
#' ## Expectation and sampling variance of eigenvalue variance
#' ## ("E[V(S)]" and "Var[V(S)]")
#' ## The argument n is for the degree of freedom, hence N - 1 in this case
#' (E_V_Sigma <- Exv.VES(Sigma = Sigma, n = N - 1))
#' (Var_V_Sigma <- Var.VES(Sigma = Sigma, n = N - 1))
#'
#' ## Same for relative eigenvalue variance ("E[Vrel(S)]", "Var[Vrel(S)]")
#' (E_Vrel_Sigma <- Exv.VRS(Sigma, N - 1))
#' (Var_Vrel_Sigma <- Var.VRS(Sigma, N - 1))
#'
#' ## Usually sample estimates are within a few S.D. away from the expectation:
#' (EDI_sam$VE - E_V_Sigma) / sqrt(Var_V_Sigma)
#' (EDI_sam$VR - E_Vrel_Sigma) / sqrt(Var_Vrel_Sigma)
#'
#' @docType package
#' @name eigvaldisp-package
#' @aliases eigvaldisp
#'
NULL
