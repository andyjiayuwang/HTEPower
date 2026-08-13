#' Solve the simplified HTE power formula
#'
#' Solve for exactly one missing quantity in the simplified power formula
#' \deqn{n = \frac{\sigma^2(z_{1-\alpha/2}+z_\kappa)^2}
#' {p(1-p)\mathrm{Var}\{s_0(X)\}}.}
#' All supplied quantities must be numeric scalars. Set exactly one argument to
#' `NULL`; the function returns the equality boundary implied by the other five
#' arguments. When `n` is unknown, the result is rounded up to the next integer.
#'
#' Because `p * (1 - p)` is symmetric around `p = 0.5`, solving for `p` can
#' produce two equivalent treatment probabilities. In that case, both are
#' returned as a named numeric vector with elements `lower` and `upper`.
#'
#' @param n Total sample size. Must be a positive integer when supplied.
#' @param alpha Two-sided significance level in `(0, 1)`.
#' @param kappa Target power in `[0.5, 1)`.
#' @param p Treatment probability in `(0, 1)`.
#' @param sigma2 Outcome disturbance variance, which must be positive.
#' @param var_s0 Variance of the CATE, `Var{s0(X)}`, which must be positive.
#'
#' @return A numeric equality boundary for the missing quantity. For an unknown
#'   `n`, the boundary is rounded up and returned as an integer. For an unknown
#'   `p`, one or two feasible treatment probabilities are returned.
#' @export
#'
#' @examples
#' hte_power(
#'   n = NULL, alpha = 0.05, kappa = 0.8, p = 0.5,
#'   sigma2 = 10, var_s0 = 0.25
#' )
#'
#' hte_power(
#'   n = 1256, alpha = 0.05, kappa = 0.8, p = NULL,
#'   sigma2 = 10, var_s0 = 0.25
#' )
hte_power <- function(n = NULL, alpha = NULL, kappa = NULL,
                      p = NULL, sigma2 = NULL, var_s0 = NULL) {
  inputs <- list(
    n = n,
    alpha = alpha,
    kappa = kappa,
    p = p,
    sigma2 = sigma2,
    var_s0 = var_s0
  )
  missing <- vapply(inputs, is.null, logical(1))

  if (sum(missing) != 1L) {
    stop("Exactly one argument must be NULL.", call. = FALSE)
  }

  supplied <- inputs[!missing]
  bad_scalar <- !vapply(
    supplied,
    function(x) is.numeric(x) && length(x) == 1L && is.finite(x),
    logical(1)
  )
  if (any(bad_scalar)) {
    stop("All supplied arguments must be finite numeric scalars.", call. = FALSE)
  }

  if (!is.null(n) && (n <= 0 || n != floor(n))) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  if (!is.null(alpha) && (alpha <= 0 || alpha >= 1)) {
    stop("`alpha` must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.null(kappa) && (kappa < 0.5 || kappa >= 1)) {
    stop("`kappa` must lie in [0.5, 1).", call. = FALSE)
  }
  if (!is.null(p) && (p <= 0 || p >= 1)) {
    stop("`p` must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.null(sigma2) && sigma2 <= 0) {
    stop("`sigma2` must be positive.", call. = FALSE)
  }
  if (!is.null(var_s0) && var_s0 <= 0) {
    stop("`var_s0` must be positive.", call. = FALSE)
  }

  unknown <- names(inputs)[missing]

  if (unknown == "n") {
    z_sum <- stats::qnorm(1 - alpha / 2) + stats::qnorm(kappa)
    required_n <- sigma2 * z_sum^2 / (p * (1 - p) * var_s0)
    return(as.integer(ceiling(required_n)))
  }

  standardized_n <- sqrt(n * p * (1 - p) * var_s0 / sigma2)

  if (unknown == "alpha") {
    z_alpha <- standardized_n - stats::qnorm(kappa)
    if (z_alpha <= 0) {
      stop(
        "The supplied values imply no feasible `alpha` in (0, 1).",
        call. = FALSE
      )
    }
    return(2 * stats::pnorm(z_alpha, lower.tail = FALSE))
  }

  if (unknown == "kappa") {
    implied_kappa <- stats::pnorm(
      standardized_n - stats::qnorm(1 - alpha / 2)
    )
    if (implied_kappa < 0.5) {
      stop(
        "The supplied values imply target power below the supported range [0.5, 1).",
        call. = FALSE
      )
    }
    return(implied_kappa)
  }

  z_sum <- stats::qnorm(1 - alpha / 2) + stats::qnorm(kappa)

  if (unknown == "p") {
    required_allocation_variance <- sigma2 * z_sum^2 / (n * var_s0)
    discriminant <- 1 - 4 * required_allocation_variance
    tolerance <- sqrt(.Machine$double.eps)
    if (discriminant < -tolerance) {
      stop(
        "No treatment probability in (0, 1) satisfies the formula for the supplied values.",
        call. = FALSE
      )
    }
    discriminant <- max(discriminant, 0)
    roots <- c(
      lower = (1 - sqrt(discriminant)) / 2,
      upper = (1 + sqrt(discriminant)) / 2
    )
    if (discriminant == 0) {
      return(unname(roots[1]))
    }
    return(roots)
  }

  if (unknown == "sigma2") {
    return(n * p * (1 - p) * var_s0 / z_sum^2)
  }

  sigma2 * z_sum^2 / (n * p * (1 - p))
}
