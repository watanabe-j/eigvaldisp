test_that("value of VE(V = GenCov(VR = VR))", {
    p <- sample(3:20, 1)
    q <- sample(1:p, 1)
    VR <- stats::runif(1, 0, (p - q) / q / (p - 1))
    Upsilon <- matrix(rnorm(p * p), p, p)
    Upsilon <- qr.Q(qr(Upsilon))

    A_spec   <- GenCov(p = p, q = q, VR = VR, evectors = Upsilon)
    A_plain  <- GenCov(p = p, q = q, VR = VR, evectors = "plain")
    A_random <- GenCov(p = p, q = q, VR = VR, evectors = "random")
    A_Givens <- GenCov(p = p, q = q, VR = VR, evectors = "Givens")
    A_MAP    <- GenCov(p = p, q = q, VR = VR, evectors = "MAP")

    VE_A_spec   <- VE(V = A_spec  )
    VE_A_plain  <- VE(V = A_plain )
    VE_A_random <- VE(V = A_random)
    VE_A_Givens <- VE(V = A_Givens)
    VE_A_MAP    <- VE(V = A_MAP   )

    expect_equal(VE_A_spec  $VR, VR)
    expect_equal(VE_A_plain $VR, VR)
    expect_equal(VE_A_random$VR, VR)
    expect_equal(VE_A_Givens$VR, VR)
    expect_equal(VE_A_MAP   $VR, VR)

    expect_equal(VE_A_spec  $VE, VR * (p - 1))
    expect_equal(VE_A_plain $VE, VR * (p - 1))
    expect_equal(VE_A_random$VE, VR * (p - 1))
    expect_equal(VE_A_Givens$VE, VR * (p - 1))
    expect_equal(VE_A_MAP   $VE, VR * (p - 1))
})

test_that("value of VE(V = GenCov(evalues = Lambda))", {
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

    VE_A_spec   <- VE(V = A_spec  )
    VE_A_plain  <- VE(V = A_plain )
    VE_A_random <- VE(V = A_random)
    VE_A_Givens <- VE(V = A_Givens)
    VE_A_MAP    <- VE(V = A_MAP   )

    expect_equal(VE_A_spec  $L, Lambda)
    expect_equal(VE_A_plain $L, Lambda)
    expect_equal(VE_A_random$L, Lambda)
    expect_equal(VE_A_Givens$L, Lambda2)
    expect_equal(VE_A_MAP   $L, Lambda2)

    expect_equal(VE_A_spec  $meanL, mean(Lambda))
    expect_equal(VE_A_plain $meanL, mean(Lambda))
    expect_equal(VE_A_random$meanL, mean(Lambda))
    expect_equal(VE_A_Givens$meanL, mean(Lambda2))
    expect_equal(VE_A_MAP   $meanL, mean(Lambda2))

    expect_equal(VE_A_spec  $VR, var(Lambda) / p / mean(Lambda) ^ 2)
    expect_equal(VE_A_plain $VR, var(Lambda) / p / mean(Lambda) ^ 2)
    expect_equal(VE_A_random$VR, var(Lambda) / p / mean(Lambda) ^ 2)
    expect_equal(VE_A_Givens$VR, var(Lambda2) / p / mean(Lambda2) ^ 2)
    expect_equal(VE_A_MAP   $VR, var(Lambda2) / p / mean(Lambda2) ^ 2)

    expect_equal(VE_A_spec  $VE, var(Lambda) * (p - 1) / p)
    expect_equal(VE_A_plain $VE, var(Lambda) * (p - 1) / p)
    expect_equal(VE_A_random$VE, var(Lambda) * (p - 1) / p)
    expect_equal(VE_A_Givens$VE, var(Lambda2) * (p - 1) / p)
    expect_equal(VE_A_MAP   $VE, var(Lambda2) * (p - 1) / p)
})

test_that("value of VE(V = GenCov(evalues = Lambda))", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    Lambda <- sort(stats::rchisq(p, 5), decreasing = TRUE)
    Sigma <- GenCov(evalues = Lambda, evectors = "random")

    X <- rmvn(N, V = Sigma)
    S <- stats::cov(X)
    R <- stats::cor(X)
    L <- svd(S, nu = 0, nv = 0)$d
    K <- svd(R, nu = 0, nv = 0)$d

    VE_X  <- VE(X = X, scale. = FALSE)
    VE_S  <- VE(V = S, scale. = FALSE)
    VE_L  <- VE(L = L)
    VE_Xs <- VE(X = X, scale. = TRUE)
    VE_Ss <- VE(V = S, scale. = TRUE)
    VE_R  <- VE(V = R)
    VE_K  <- VE(L = K)

    expect_equal(VE_X, VE_S)
    expect_equal(VE_X, VE_L)
    expect_equal(VE_Xs, VE_Ss)
    expect_equal(VE_Xs, VE_R)
    expect_equal(VE_Xs, VE_K)
})
