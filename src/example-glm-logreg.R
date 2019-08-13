reprex::reprex({
  # Arnold, T., Kane, M., & Lewis, B. W. (2019). A Computational Approach to Statistical Learning. CRC Press.
  # GLM with Newton Raphson and QR but fixed to logistic regression
  # I made some changes compared to original
  # License LGPL I believe: https://github.com/statsmaths/casl/blob/master/R/ch05-glm.R
  library(armacmp)
  glm_fit <- function(X, y = type_colvec(), maxit = type_scalar_int(), tol = type_scalar_numeric()) {
    s <- rep.int(0, ncol(X))
    eta <- rep.int(0, nrow(X))
    QR <- qr(X)
    Q <- qr.Q(QR)
    R <- qr.R(QR)
    linkinv <- function(eta = type_colvec()) {
      return(exp(eta) / (1 + exp(eta)))
    }
    mu_eta <- function(eta = type_colvec()) {
      return(exp(eta) / (1 + exp(eta))^2)
    }
    variance <- function(mu = type_colvec()) {
      return(mu * (1 - mu))
    }
    for (j in seq_len(maxit)) {
      s_old <- s
      mu <- linkinv(eta)
      mu_p <- mu_eta(eta)
      z <- eta + (y - mu) / mu_p
      W <- mu_p^2 / variance(mu)
      C <- chol(crossprod(Q, diag(W) %*% Q))
      s1 <- forwardsolve(t(C), crossprod(Q, W * z))
      s <- backsolve(C, s1)
      eta <- Q %*% s
      if (sum(sqrt(crossprod(s - s_old))) < tol) break
    }
    beta <- backsolve(R, crossprod(Q, eta))
    return(beta)
  }

  glm_fit_cpp <- compile(glm_fit, verbose = TRUE)

  set.seed(50)
  n <- 1e6
  p <- 10
  true_beta <- rnorm(p)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), ncol = p - 1))
  y <- runif(n) < plogis(X %*% true_beta)

  coefs <- glm_fit_cpp(X, y, 25, 1e-6)
  coefs2 <- coef(glm.fit(X, y, family = binomial()))

  all.equal(
    as.numeric(coefs),
    as.numeric(coefs2)
  )

  # Note: glm.fit does more than simply computing the coefficients
  system.time(glm_fit_cpp(X, y, 25, 1e-6))
  system.time(glm.fit(X, y, family = binomial()))
}, outfile = "example-glm-logreg", venue = "R", show = FALSE)
