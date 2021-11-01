test_that("behavior of crossprod(chol_xx)", {
    A1 <- matrix(0, 4, 4)
    diag(A1) <- 1
    A1[c(1, 2), c(2, 1)] <- 1
    A2 <- A1
    A2[c(3, 4), c(4, 3)] <- 1
    A3 <- matrix(0, 4, 4)
    diag(A3) <- 1
    A3[c(1, 4), c(4, 1)] <- 1
    A4 <- A3
    A4[c(2, 3), c(3, 2)] <- 1
    A5 <- A6 <- matrix(1, 4, 4)
    A5[1:3, 4] <- A5[4, 1:3] <- 0

    expect_equal(crossprod(matsqrt(A1)), A1)
    expect_equal(crossprod(chol_qr(A1)), A1)
    expect_equal(crossprod(chol_piv(A1)), A1)

    expect_equal(crossprod(matsqrt(A2)), A2)
    expect_equal(crossprod(chol_qr(A2)), A2)
    expect_warning(chol_piv(A2))

    expect_equal(crossprod(matsqrt(A3)), A3)
    expect_equal(crossprod(chol_piv(A3)), A3)
    expect_equal(crossprod(chol_qr(A3)), A3)

    expect_equal(crossprod(matsqrt(A4)), A4)
    expect_equal(crossprod(chol_qr(A4)), A4)
    expect_warning(chol_piv(A4))

    expect_equal(crossprod(matsqrt(A5)), A5)
    expect_equal(crossprod(chol_qr(A5)), A5)
    expect_warning(chol_piv(A5))

    expect_equal(crossprod(matsqrt(A6)), A6)
    expect_equal(crossprod(chol_qr(A6)), A6)
    expect_warning(chol_piv(A6))
})
