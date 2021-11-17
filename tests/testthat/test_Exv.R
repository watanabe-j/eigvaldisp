test_that("behavior of Exv.VXx under null", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    s2 <- stats::rchisq(1, 5)
    I <- diag(p)
    J <- s2 * I

    expect_equal(Exv.VES(I, n), (p - 1) * (p + 2) / (p * n))
    expect_equal(Exv.VES(J, n), (p - 1) * (p + 2) / (p * n) * s2 ^ 2)
    expect_equal(Exv.VRS(I, n), (p + 2) / (p * n + 2))
    expect_equal(Exv.VRS(J, n), (p + 2) / (p * n + 2))
    expect_equal(Exv.VER(I, n), (p - 1) / n)
    expect_equal(Exv.VRR(I, n), 1 / n)

    expect_equal(Exv.VESa(I, n), 0)
    expect_equal(Exv.VESa(J, n), 0)
    expect_equal(Exv.VRSa(I, n), 0)
    expect_equal(Exv.VRSa(J, n), 0)
    expect_equal(Exv.VRRa(I, n), 0)
})

test_that("behavior of Var.VXx under null", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    s2 <- stats::rchisq(1, 5)
    I <- diag(p)
    J <- s2 * I

    v1 <- 4 * (p - 1) * (p + 2) * (2 * p ^ 2 + p * n + 3 * p - 6) /
        (p ^ 3 * n ^ 3)
    v2 <- 4 * p ^ 2 * (p + 2) * (n - 1) * (n + 2) /
        ((p - 1) * (p * n + 2) ^ 2 * (p * n + 4) * (p * n + 6))
    v3 <- 4 * (n - 1) / (p * (p - 1) * n ^ 2 * (n + 2))

    v4 <- 4 * (p * n + 2) * (p - 1) * (p + 2) / (p ^ 3 * n * (n - 1) * (n + 2))

    expect_equal(Var.VES(I, n), v1)
    expect_equal(Var.VES(J, n), v1 * s2 ^ 4)
    expect_equal(Var.VRS(I, n), v2)
    expect_equal(Var.VRS(J, n), v2)
    expect_equal(Var.VER(I, n), v3 * (p - 1) ^ 2)
    expect_equal(Var.VRR(I, n), v3)

    expect_equal(Var.VESa(I, n), v4)
    expect_equal(Var.VESa(J, n), v4 * s2 ^ 4)
    expect_equal(Var.VRSa(I, n), v2 * ((p * n + 2) / (p * (n - 1))) ^ 2)
    expect_equal(Var.VRSa(J, n), v2 * ((p * n + 2) / (p * (n - 1))) ^ 2)
    expect_equal(Var.VRRa(I, n), v3 * (n / (n - 1)) ^ 2)

})

test_that("behavior of Exv/Var.VRx under complete integration", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    s2 <- stats::rchisq(1, 5)
    Rho <- matrix(1, p, p)
    Sigma <- Rho * s2

    expect_equal(Exv.VRS(Rho, n), 1)
    expect_equal(Exv.VRS(Sigma, n), 1)
    expect_equal(Exv.VRR(Rho, n), 1)

    expect_equal(Var.VRS(Rho, n), 0)
    expect_equal(Var.VRS(Sigma, n), 0)
    expect_equal(Var.VRR(Rho, n), 0)

    expect_equal(AVar.VRR_pfd(Rho, n), 0)
    expect_equal(AVar.VRR_klv(Rho, n), 0)

    expect_equal(Exv.VESa(Rho, n), p - 1)
    expect_equal(Exv.VRSa(Rho, n), 1)
    expect_equal(Exv.VRSa(Sigma, n), 1)
    expect_equal(Exv.VRRa(Rho, n), 1)

    expect_equal(Var.VRSa(Rho, n), 0)
    expect_equal(Var.VRSa(Sigma, n), 0)
    expect_equal(Var.VRRa(Rho, n), 0)
})

test_that("behavior of AVar.VRR_kxx under null", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    I <- diag(p)

    expect_equal(AVar.VRR_klv(I, n), 0)
    expect_equal(AVar.VRR_kl(I, n), 0)
    expect_equal(AVar.VRR_kr(I, n), 0)
    expect_equal(AVar.VRR_krv(I, n), 0)
})

test_that("behavior of AVar.VRR_kxx under arbitrary conditions", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    Lambda <- sort(stats::rchisq(p, 5), decreasing = TRUE)
    Lambda <- Lambda / sum(Lambda) * p
    A <- GenCov(evalues = Lambda, evectors = "Givens")

    klv <- AVar.VRR_klv(A, n)
    kl <- AVar.VRR_kl(A, n)
    kr <- AVar.VRR_kr(A, n)
    krv <- AVar.VRR_krv(A, n)
    expect_equal(klv, kl)
    expect_equal(klv, kr)
    expect_equal(klv, krv)
})


test_that("behavior of AVar.VRR_pfx under arbitrary conditions", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    Lambda <- sort(stats::rchisq(p, 5), decreasing = TRUE)
    Lambda <- Lambda / sum(Lambda) * p
    A <- GenCov(evalues = Lambda, evectors = "Givens")

    pfd <- AVar.VRR_pfd(A, n)
    pfv <- AVar.VRR_pfv(A, n)
    pf <- AVar.VRR_pf(A, n)
    pfc <- AVar.VRR_pfc(A, n)
    expect_equal(pfd, pfv)
    expect_equal(pfd, pf)
    expect_equal(pfd, pfc)
})
