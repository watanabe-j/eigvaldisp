#include <Rcpp.h>
using namespace Rcpp;

//*** Cov.r2C ***
//' Sum of covariances of squared correlations
//'
//' \code{Cov.r2C()} calculates (twice) sum of covariances of
//' all non-redundant pairs of squared correlation coefficients.
//' This is an internal C++ function to be used in within \code{AVar.VRr_pfc()}
//' with the package \code{Rcpp}.
//'
//' Primarily for internal use.
//'
//' It is assumed that \code{E} has the form of
//' \code{matrix(sapply(n, Exv.r1, x = R[upper.tri(R)]), ncol = length(n))}.
//'
//' \code{R} and \code{E} are read as NumericVector rather than
//' NumericMatrix to save a modest amount of time
//'
//' @param n
//'   Degrees of freedom (not sample sizes); numeric of length 1 or more.
//' @param R
//'   Population correlation matrix; assumed to be validly constructed;
//'   numeric of length \code{p * p}
//' @param E
//'   Matrix of expectation of correlation coefficients corresponding to
//'   the upper triangular of \code{R}; numeric of length
//'   \code{p * (p - 1) / 2 * length(n)}.
//'
//' @return
//' A numeric vector of \eqn{\sum 2 Cov(r_{ij}^2, r_{kl}^2)},
//' corresponding to \code{n}.
//'
//' @seealso \link{Exv.VXX}, \link{AVar.VRR_xx}
//'
// [[Rcpp::export]]
NumericVector Cov_r2C(NumericVector n, NumericVector R, NumericVector E) {
    int ln = n.size();
    int p = sqrt(R.size()); // If R was read as NumericMatrix, int p = R.nrow()
    int nrE = p * (p - 1) / 2;
    double rij, rkl, rik, ril, rjk, rjl, C;
    double Eij, Ekl, Cijkl;
    NumericVector out(ln);
    for(int i = 0; i < p - 1; i++) {
        for(int j = i + 1; j < p; j++) {
            for(int k = i; k < p - 1; k++) {
                for(int l = ((i == k) ? (j + 1) : (k + 1)); l < p; l++) {
                    rij = R[i + j * p];  // If R was a NumericMatrix, R(i, j)
                    rkl = R[k + l * p];  // Similar for Eij below
                    rik = R[i + k * p];
                    rjl = R[j + l * p];
                    ril = R[i + l * p];
                    rjk = R[j + k * p];
                    C = rij * rkl * (rik * rik + ril * ril +
                                     rjk * rjk + rjl * rjl) / 2 +
                        rik * rjl + ril * rjk - (rij * rik * ril
                        + rij * rjk * rjl + rik * rjk * rkl + ril * rjl * rkl);
                    for(int m = 0; m < ln; m++) {
                        Eij = E[i + j * p - j * (2 * p - j + 1) / 2 + m * nrE];
                        Ekl = E[k + l * p - l * (2 * p - l + 1) / 2 + m * nrE];
                        Cijkl = C / n[m];
                        out[m] += (4 * Eij * Ekl + 2 * Cijkl) * Cijkl;
                    }
                }
            }
        }
    } // for speed, out is doubled after skipping redundant pairs of (rij, rkl)
    out = out * 2;
    return out;
}
