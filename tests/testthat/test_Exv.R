test_that("behavior of Exv.VXx under null", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    s2 <- stats::rchisq(1, 5)
    I <- diag(p)
    J <- s2 * I

    expect_equal(Exv.VEv(I, n), (p - 1) * (p + 2) / (p * n))
    expect_equal(Exv.VEv(J, n), (p - 1) * (p + 2) / (p * n) * s2 ^ 2)
    expect_equal(Exv.VRv(I, n), (p + 2) / (p * n + 2))
    expect_equal(Exv.VRv(J, n), (p + 2) / (p * n + 2))
    expect_equal(Exv.VEr(I, n), (p - 1) / n)
    expect_equal(Exv.VRr(I, n), 1 / n)

    expect_equal(Exv.VEav(I, n), 0)
    expect_equal(Exv.VEav(J, n), 0)
    expect_equal(Exv.VRav(I, n), 0)
    expect_equal(Exv.VRav(J, n), 0)
    expect_equal(Exv.VRar(I, n), 0)
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

    expect_equal(Var.VEv(I, n), v1)
    expect_equal(Var.VEv(J, n), v1 * s2 ^ 4)
    expect_equal(Var.VRv(I, n), v2)
    expect_equal(Var.VRv(J, n), v2)
    expect_equal(Var.VEr(I, n), v3 * (p - 1) ^ 2)
    expect_equal(Var.VRr(I, n), v3)

    expect_equal(Var.VEav(I, n), v4)
    expect_equal(Var.VEav(J, n), v4 * s2 ^ 4)
    expect_equal(Var.VRav(I, n), v2 * ((p * n + 2) / (p * (n - 1))) ^ 2)
    expect_equal(Var.VRav(J, n), v2 * ((p * n + 2) / (p * (n - 1))) ^ 2)
    expect_equal(Var.VRar(I, n), v3 * (n / (n - 1)) ^ 2)

})

test_that("behavior of Exv/Var.VRx under complete integration", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    s2 <- stats::rchisq(1, 5)
    Rho <- matrix(1, p, p)
    Sigma <- Rho * s2

    expect_equal(Exv.VRv(Rho, n), 1)
    expect_equal(Exv.VRv(Sigma, n), 1)
    expect_equal(Exv.VRr(Rho, n), 1)

    expect_equal(Var.VRv(Rho, n), 0)
    expect_equal(Var.VRv(Sigma, n), 0)
    expect_equal(Var.VRr(Rho, n), 0)

    expect_equal(AVar.VRr_pfd(Rho, n), 0)
    expect_equal(AVar.VRr_klv(Rho, n), 0)

    expect_equal(Exv.VEav(Rho, n), p - 1)
    expect_equal(Exv.VRav(Rho, n), 1)
    expect_equal(Exv.VRav(Sigma, n), 1)
    expect_equal(Exv.VRar(Rho, n), 1)

    expect_equal(Var.VRav(Rho, n), 0)
    expect_equal(Var.VRav(Sigma, n), 0)
    expect_equal(Var.VRar(Rho, n), 0)
})

test_that("behavior of AVar.VRr_kxx under null", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    I <- diag(p)

    expect_equal(AVar.VRr_klv(I, n), 0)
    expect_equal(AVar.VRr_kl(I, n), 0)
    expect_equal(AVar.VRr_kr(I, n), 0)
    expect_equal(AVar.VRr_krv(I, n), 0)
})

test_that("behavior of AVar.VRr_kxx under arbitrary conditions", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    Lambda <- sort(stats::rchisq(p, 5), decreasing = TRUE)
    Lambda <- Lambda / sum(Lambda) * p
    A <- GenCov(evalues = Lambda, evectors = "Givens")

    klv <- AVar.VRr_klv(A, n)
    kl <- AVar.VRr_kl(A, n)
    kr <- AVar.VRr_kr(A, n)
    krv <- AVar.VRr_krv(A, n)
    expect_equal(klv, kl)
    expect_equal(klv, kr)
    expect_equal(klv, krv)
})


test_that("behavior of AVar.VRr_pfx under arbitrary conditions", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    Lambda <- sort(stats::rchisq(p, 5), decreasing = TRUE)
    Lambda <- Lambda / sum(Lambda) * p
    A <- GenCov(evalues = Lambda, evectors = "Givens")

    pfd <- AVar.VRr_pfd(A, n)
    pfv <- AVar.VRr_pfv(A, n)
    pf <- AVar.VRr_pf(A, n)
    pfc <- AVar.VRr_pfc(A, n)
    expect_equal(pfd, pfv)
    expect_equal(pfd, pf)
    expect_equal(pfd, pfc)
})
