test_that("sample size is rounded up", {
  expect_identical(
    hte_power(
      n = NULL, alpha = 0.05, kappa = 0.8, p = 0.5,
      sigma2 = 10, var_s0 = 0.25
    ),
    1256L
  )
})

test_that("each scalar quantity can be recovered", {
  n <- 1308
  alpha <- 0.05
  kappa <- 0.8
  p <- 0.4
  sigma2 <- 10
  var_s0 <- 0.25
  z_sum <- qnorm(1 - alpha / 2) + qnorm(kappa)
  standardized_n <- sqrt(n * p * (1 - p) * var_s0 / sigma2)

  expect_equal(
    hte_power(n, NULL, kappa, p, sigma2, var_s0),
    2 * pnorm(standardized_n - qnorm(kappa), lower.tail = FALSE)
  )
  expect_equal(
    hte_power(n, alpha, NULL, p, sigma2, var_s0),
    pnorm(standardized_n - qnorm(1 - alpha / 2))
  )
  expect_equal(
    hte_power(n, alpha, kappa, p, NULL, var_s0),
    n * p * (1 - p) * var_s0 / z_sum^2
  )
  expect_equal(
    hte_power(n, alpha, kappa, p, sigma2, NULL),
    sigma2 * z_sum^2 / (n * p * (1 - p))
  )
})

test_that("solving for treatment probability returns symmetric roots", {
  z_sum <- qnorm(1 - 0.05 / 2) + qnorm(0.8)
  n <- 1308
  roots <- hte_power(
    n = n, alpha = 0.05, kappa = 0.8, p = NULL,
    sigma2 = 10, var_s0 = 0.25
  )

  expect_named(roots, c("lower", "upper"))
  expect_equal(sum(roots), 1)
  expect_equal(
    roots[["lower"]] * roots[["upper"]],
    10 * z_sum^2 / (n * 0.25)
  )
})

test_that("invalid leave-one-unknown calls fail clearly", {
  expect_error(hte_power(), "Exactly one")
  expect_error(
    hte_power(NULL, NULL, 0.8, 0.5, 10, 0.25),
    "Exactly one"
  )
  expect_error(
    hte_power(NULL, 0.05, 0.8, 1, 10, 0.25),
    "strictly between"
  )
  expect_error(
    hte_power(1256.5, 0.05, 0.8, 0.5, NULL, 0.25),
    "positive integer"
  )
  expect_error(
    hte_power(0, 0.05, 0.8, 0.5, NULL, 0.25),
    "positive integer"
  )
})
