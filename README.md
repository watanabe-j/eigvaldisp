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
[pandoc](https://pandoc.org) (including pandoc-citeproc) installed
on your machine, you can build a vignette by the option
`build_vignettes = TRUE` in `install_github()` (recommended).

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

At present, one function uses the C++ API via `Rcpp`.
Windows users will require
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)
at installation to compile C++ codes.
If you want to try this package without installing them,
simply `source()` the .R files in `R/`, except `eigvaldisp-package.R`,
`RcppExports.R`, and `zzz.R`.


## Examples

To get some ideas, let's create a simple population covariance matrix
with the function `GenCov()`, which constructs a covariance/correlation
matrix with known eigenvalues/vectors, and then calculate
eigenvalue dispersion indices of this matrix with the function `VE()`:
```
set.seed(30)
## Generate a population covariance matrix with known eigenvalues
Lambda <- c(4, 2, 1, 1)
(Sigma <- GenCov(evalues = Lambda, evectors = "random"))
#>              [,1]         [,2]       [,3]       [,4]
#> [1,]  2.525707231 -0.002305745  0.4878674 -1.2925039
#> [2,] -0.002305745  1.893637883  0.2794045 -0.4629743
#> [3,]  0.487867373  0.279404464  1.2438233 -0.5590454
#> [4,] -1.292503917 -0.462974308 -0.5590454  2.3368316
eigen(Sigma)$values
#> [1] 4 2 1 1

## Calculate eigenvalue dispersion indices of this matrix
EDI_pop <- VE(S = Sigma)

## Eigenvalue variance ("V(Sigma)")
EDI_pop$VE
#> [1] 1.5

## Eigenvalue variance ("Vrel(Sigma)"):
EDI_pop$VR
#> [1] 0.125
```

It is trivial to calculate the population eigenvalue dispersion indices.
The problem is the presence of sampling bias (and error),
which renders inferences from a sample rather difficult.

To see this, simulate a small multivariate normal sample from
the same population covariance matrix using
the internal function `rmvn()`:
```
## Simulate a multivariate normal sample
N <- 20
X <- eigvaldisp:::rmvn(N = N, Sigma = Sigma)
cov(X)
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  2.8116163  0.5779609  0.9517708 -1.3106853
#> [2,]  0.5779609  2.6693532  0.2962871 -0.9844253
#> [3,]  0.9517708  0.2962871  1.2897573 -0.5258998
#> [4,] -1.3106853 -0.9844253 -0.5258998  2.2849262
## Reasonable estimate of Sigma

## Calculating eigenvalue dispersion indices from the sample
EDI_sam <- VE(X = X)
## Same as VE(S = cov(X)) but usually faster

## Sample eigenvalue variance ("V(S)")
EDI_sam$VE
#> [1] 2.499072

## Sample relative eigenvalue variance ("Vrel(S)")
EDI_sam$VR
#> [1] 0.1625316
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
## Expectation of eigenvalue variance ("E[V(S)]")
## The argument n is for the degree of freedom, hence N - 1 in this case
(E_V_Sigma <- Exv.VES(Sigma = Sigma, n = N - 1))
#> [1] 2.486842

## Expected bias
E_V_Sigma - EDI_pop$VE
#> [1] 0.9868421

## Error (sampling variance) of eigenvalue variance ("Var[V(S)]")
Var.VES(Sigma, N - 1)
#> [1] 3.126513

## Same for relative eigenvalue variance ("E[Vrel(S)]", "Var[Vrel(S)]")
(E_Vrel_Sigma <- Exv.VRS(Sigma, N - 1))
#> [1] 0.18438
E_Vrel_Sigma - EDI_pop$VR
#> [1] 0.05938
Var.VRS(Sigma, N - 1)
#> [1] 0.007989498
```

"Bias-corrected" estimators are also implemented, although this is
not globally unbiased for the relative eigenvalue variance:
```
## Bias-corrected eigenvalue variance
VESa(X = X)$VESa
#> [1] 1.290192

## Its expectation (equals the population value)
Exv.VESa(Sigma, N - 1)
#> [1] 1.5

## Its variance (smaller than that of the ordinary one)
Var.VESa(Sigma, N - 1)
#> [1] 2.094455


## Adjusted relative eigenvalue variance
VRSa(X = X)$VRSa
#> [1] 0.09274259

## Its expectation (underestimates the population value)
Exv.VRSa(Sigma, N - 1)
#> [1] 0.1164117

## Its variance
Var.VRSa(Sigma, N - 1)
#> [1] 0.009376563
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
The function `GenCov()` in this package involves an algorithm originally
adopted (with modifications) from the package `fungible` version 1.99
(Waller, 2021), which is under GPL (>= 2). See description of that function for details.


## References
Waller, N. G. (2021). fungible: psychometric functions from
the Waller Lab. R package version 1.99.
[https://CRAN.R-project.org/package=fungible](https://CRAN.R-project.org/package=fungible).

Watanabe, J. (2021). Statistics of eigenvalue dispersion indices: quantifying the magnitude of phenotypic integration. *Evolution*, doi:[10.1111/evo.14382](https://doi.org/10.1111/evo.14382).
