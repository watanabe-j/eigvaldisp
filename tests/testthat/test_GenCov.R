test_that("symmetry of GenCov()", {
    check_symmetry <- function(A) {
        expect_equal(A[lower.tri(A)], t(A)[lower.tri(A)])
    }
    p <- sample(3:20, 1)
    q <- sample(1:p, 1)
    VR <- stats::runif(1, 0, (p - q) / q / (p - 1))
    stopifnot(VR >= 0 && VR <= 1)
    Upsilon <- rnorm(p)
    Upsilon <- Upsilon / sqrt(sum(Upsilon ^ 2))

    A_spec   <- GenCov(p = p, VR = VR, evectors = Upsilon)
    A_plain  <- GenCov(p = p, VR = VR, evectors = "plain")
    A_random <- GenCov(p = p, VR = VR, evectors = "random")
    A_Givens <- GenCov(p = p, VR = VR, evectors = "Givens")
    A_MAP    <- GenCov(p = p, VR = VR, evectors = "MAP")

    check_symmetry(A_plain)
    check_symmetry(A_random)
    check_symmetry(A_Givens)
    check_symmetry(A_MAP)
    check_symmetry(A_spec)
})

test_that("elements of GenCov(evectors = c('plain', 'Givens', 'MAP'))", {
    p <- sample(3:20, 1)
    q <- sample(1:p, 1)
    VR <- stats::runif(1, 0, (p - q) / q / (p - 1))
    stopifnot(VR >= 0 && VR <= 1)

    A_plain  <- GenCov(p = p, VR = VR, evectors = "plain")
    A_Givens <- GenCov(p = p, VR = VR, evectors = "Givens")
    A_MAP    <- GenCov(p = p, VR = VR, evectors = "MAP")

    expect_equal(A_plain[lower.tri(A_plain)], rep.int(0, p * (p - 1) / 2))
    expect_equal(diag(A_Givens), rep.int(1, p))
    expect_equal(diag(A_MAP), rep.int(1, p))
})

test_that("structure of GenCov(q = 1, evectors = c('Givens', 'MAP'))", {
    p <- sample(3:20, 1)
    VR <- stats::runif(1, 0, 1)

    A_Givens <- GenCov(p = p, q = 1, VR = VR, evectors = "Givens")
    A_MAP    <- GenCov(p = p, q = 1, VR = VR, evectors = "MAP")
    evec_A_Givens <- svd(A_Givens, nu = 0, nv = 1)$v[, 1]
    evec_A_MAP    <- svd(A_MAP, nu = 0, nv = 1)$v[, 1]

    expect_equal(A_Givens[lower.tri(A_Givens)]^2, rep.int(VR, p * (p - 1) / 2))
    expect_equal(A_MAP[lower.tri(A_MAP)]^2, rep.int(VR, p * (p - 1) / 2))
    expect_equal(evec_A_Givens ^ 2, rep.int(1 / p, p))
    expect_equal(evec_A_MAP ^ 2, rep.int(1 / p, p))
})

test_that("eigenvalues of GenCov(evalues = X)", {
    p <- sample(3:20, 1)
    Lambda <- sort(stats::rchisq(p, 5), decreasing = TRUE)
    Lambda2 <- Lambda / sum(Lambda) * p
    Upsilon <- rnorm(p)
    Upsilon <- Upsilon / sqrt(sum(Upsilon ^ 2))

    A_spec   <- GenCov(evalues = Lambda, evectors = Upsilon)
    A_plain  <- GenCov(evalues = Lambda, evectors = "plain")
    A_random <- GenCov(evalues = Lambda, evectors = "random")
    A_Givens <- GenCov(evalues = Lambda2, evectors = "Givens")
    A_MAP    <- GenCov(evalues = Lambda2, evectors = "MAP")

    expect_equal(svd(A_spec,   nu = 0, nv = 0)$d, Lambda)
    expect_equal(svd(A_plain,  nu = 0, nv = 0)$d, Lambda)
    expect_equal(svd(A_random, nu = 0, nv = 0)$d, Lambda)
    expect_equal(svd(A_Givens, nu = 0, nv = 0)$d, Lambda2)
    expect_equal(svd(A_MAP,    nu = 0, nv = 0)$d, Lambda2)
})

test_that("eigenvectors of GenCov(evectors = X) when X is full", {
    p <- sample(3:20, 1)
    # q <- sample(1:p, 1)
    # VR <- stats::runif(1, 0, (p - q) / q / (p - 1))
    # stopifnot(VR >= 0 && VR <= 1)
    Lambda <- sort(stats::rchisq(p, 5), decreasing = TRUE)
    Upsilon <- matrix(rnorm(p * p), p, p)
    Upsilon <- qr.Q(qr(Upsilon))

    A1 <- GenCov(evalues = Lambda, evectors = Upsilon)
    A2 <- GenCov(p = p, shape = "linearly_decreasing", evectors = Upsilon)

    expect_equal(diag(crossprod(svd(A1, nu = 0)$v, Upsilon)) ^ 2, rep.int(1, p))
    expect_equal(diag(crossprod(svd(A2, nu = 0)$v, Upsilon)) ^ 2, rep.int(1, p))
})

test_that("eigenvectors of GenCov(evectors = X) when X is partial", {
    p <- sample(3:20, 1)
    q <- sample(1:(p - 1), 1)
    # VR <- stats::runif(1, 0, (p - q) / q / (p - 1))
    # stopifnot(VR >= 0 && VR <= 1)
    Lambda <- sort(stats::rchisq(p, 5), decreasing = TRUE)
    Upsilon <- matrix(rnorm(p * p), p, p)
    Upsilon <- qr.Q(qr(Upsilon))[, 1:q]

    A1 <- GenCov(evalues = Lambda, evectors = Upsilon)
    A2 <- GenCov(p = p, shape = "linearly_decreasing", evectors = Upsilon)

    expect_equal(diag(crossprod(svd(A1, nu = 0, nv = q)$v, Upsilon)) ^ 2, rep.int(1, q))
    expect_equal(diag(crossprod(svd(A2, nu = 0, nv = q)$v, Upsilon)) ^ 2, rep.int(1, q))

    expect_equal(crossprod(svd(A1, nu = 0)$v), diag(p))
    expect_equal(crossprod(svd(A2, nu = 0)$v), diag(p))
})
