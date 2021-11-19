# eigvaldisp
R package for statistics of eigenvalue dispersion indices

This package involves functions for analyzing eigenvalue dispersion
indices of covariance and correlation matrices&mdash;common measures
of phenotypic integration in biometrics.
The primary feature of this package is to calculate expectation and
variance (i.e., sampling bias and error) of eigenvalue dispersion indices
for arbitrary population covariance structures
under multivariate normality.

This package was designed to supplement Watanabe (2021),
and built on the supplementary scripts associated with that paper.
See that paper for theoretical details.


## Installation
```
# install.packages("devtools")
devtools::install_github("watanabe-j/eigvaldisp")
```
If you have the packages `rmarkdown` and `knitr`, and have
[pandoc](https://pandoc.org) installed on your machine, you can build
a vignette by the option `build_vignettes = TRUE` in `install_github()`
(not necessary, but recommended).

This package has the following dependencies:
```
Imports:
    stats,
    hypergeo,
    Rcpp
Suggests:
    parallel,
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
LinkingTo:
    Rcpp
```

`stats` and `hypergeo` are strictly necessary, whereas
`Rcpp` is used only in a single function.
If you want to try this package without installing `Rcpp`,
simply `source()` the .R files in `R/`, except `eigvaldisp-package.R`,
`RcppExports.R`, and `zzz.R`.


## Examples

To get some ideas, let's create a simple population covariance matrix
with the function `GenCov()`, which constructs a covariance/correlation
matrix with known eigenvalues/vectors, and then calculate
eigenvalue dispersion indices of this matrix with the function `VE()`:
```
# Creating a population covariance matrix with known eigenvalues
Lambda <- c(4, 2, 1, 1)
Sigma <- GenCov(evalues = Lambda, evectors = "random")

# Calculating eigenvalue dispersion indices of this matrix
EDI_pop <- VE(S = Sigma)

# Eigenvalue variance ("V(Sigma)"): 1.5
EDI_pop$VE

# Eigenvalue variance ("Vrel(Sigma)"): 0.125
EDI_pop$VR
```

It is trivial to calculate the population eigenvalue dispersion indices.
The problem is the presence of sampling bias (and error),
which renders inferences from a sample rather difficult.

To see this, simulate a small multivariate normal sample from
the same population covariance matrix using
the internal function `rmvn()`:
```
# Simulate a multivariate normal sample
N <- 20
X <- eigvaldisp:::rmvn(N = N, Sigma = Sigma)

# Calculating eigenvalue dispersion indices from the sample
EDI_sam <- VE(X = X)
# Same as VE(S = cov(X)) but faster

# Sample eigenvalue variance ("V(S)")
EDI_sam$VE

# Sample relative eigenvalue variance ("Vrel(S)")
EDI_sam$VR
```

These are typically larger than the population values,
although there is always some random fluctuation. In other words,
sample eigenvalue dispersion indices tend to overestimate
the population values in this case
(although underestimation can happen when the population value is large).

The main functionality of this package is to calculate
expectation and variance (i.e., estimates of sampling bias and error)
of eigenvalue dispersion indices from arbitrary
population covariance/correlation matrices:
```
# Expectation of eigenvalue variance ("E[V(S)]"): 2.487
# The argument n is for the degree of freedom, hence N - 1
(E_V_Sigma <- Exv.VES(Sigma, N - 1))

# Expected bias: 0.987
E_V_Sigma - EDI_pop$VE

# Error (sampling variance) of eigenvalue variance ("Var[V(S)]"): 3.127
Var.VES(Sigma, N - 1)

# Same for relative eigenvalue variance ("E[Vrel(S)]", "Var[Vrel(S)]"):
# 0.184, 0.059, and 0.008
(E_Vrel_Sigma <- Exv.VRS(Sigma, N - 1))
E_Vrel_Sigma - EDI_pop$VR
Var.VRS(Sigma, N - 1)
```

"Bias-corrected" estimators are also implemented, although this is
not globally unbiased for the relative eigenvalue variance:
```
# Bias-corrected eigenvalue variance
VESa(X = X)$VESa

# Its expectation: 1.5 (as this is unbiased)
Exv.VESa(Sigma, N - 1)

# Its variance: 2.094 (smaller than that of the ordinary one)
Var.VESa(Sigma, N - 1)


# Adjusted relative eigenvalue variance,
VRSa(X = X)$VRSa

# Its expectation: 0.116 (tendency for underestimation)
Exv.VRSa(Sigma, N - 1)

# Its variance: 0.009 (larger than that of the ordinary one)
Var.VRSa(Sigma, N - 1)
```

The same functionalities are also available for correlation matrices,
but via different functions (`Exv.VRR()`, `Var.VRR()`, `VRRa()`, etc.)
since their distributions are different.

Also involved in this package are:
- Function for Monte Carlo simulation of eigenvalue dispersion indices
- Cholesky factorization functions for singular covariance matrices
- Calculation of exact/asymptotic moments of correlation coefficients

There may be better R implementations for some of the functionalities,
but this package is intended to be as much self-contained as possible.

For more detailed descriptions, use `vignette("eigvaldisp")`.


## Copyright notice
The function `GenCov()` in this package involves algorithms originally
adopted (with modifications) from the package `fungible` version 1.99
(Waller, 2021), which is under GPL (>= 2). See description of that function for details.


## References
Waller, N. G. (2021). fungible: psychometric functions from
the Waller Lab. R package version 1.99.
[https://CRAN.R-project.org/package=fungible](https://CRAN.R-project.org/package=fungible).

Watanabe, J. (2021). Statistics of eigenvalue dispersion indices: quantifying the magnitude of phenotypic integration. *Evolution*, doi:[10.1111/evo.14382](https://doi.org/10.1111/evo.14382).
