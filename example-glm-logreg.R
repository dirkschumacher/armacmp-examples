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
#> R function
#> 
#> function (X, y = type_colvec(), maxit = type_scalar_int(), tol = type_scalar_numeric()) 
#> {
#>     s <- rep.int(0, ncol(X))
#>     eta <- rep.int(0, nrow(X))
#>     QR <- qr(X)
#>     Q <- qr.Q(QR)
#>     R <- qr.R(QR)
#>     linkinv <- function(eta = type_colvec()) {
#>         return(exp(eta)/(1 + exp(eta)))
#>     }
#>     mu_eta <- function(eta = type_colvec()) {
#>         return(exp(eta)/(1 + exp(eta))^2)
#>     }
#>     variance <- function(mu = type_colvec()) {
#>         return(mu * (1 - mu))
#>     }
#>     for (j in seq_len(maxit)) {
#>         s_old <- s
#>         mu <- linkinv(eta)
#>         mu_p <- mu_eta(eta)
#>         z <- eta + (y - mu)/mu_p
#>         W <- mu_p^2/variance(mu)
#>         C <- chol(crossprod(Q, diag(W) %*% Q))
#>         s1 <- forwardsolve(t(C), crossprod(Q, W * z))
#>         s <- backsolve(C, s1)
#>         eta <- Q %*% s
#>         if (sum(sqrt(crossprod(s - s_old))) < tol) 
#>             break
#>     }
#>     beta <- backsolve(R, crossprod(Q, eta))
#>     return(beta)
#> }
#> 
#> C++ function translation
#> 
#> arma::colvec armacmp_fun(const arma::mat& X, const arma::colvec& y, int maxit, double tol)
#> {
#> arma::colvec s = arma::colvec(std::vector<double>(X.n_cols, 0.0));
#> arma::colvec eta = arma::colvec(std::vector<double>(X.n_rows, 0.0));
#> arma::mat QR__Q, QR__R;
#> arma::qr_econ( QR__Q, QR__R, X );
#> 
#> arma::mat Q = QR__Q;
#> arma::mat R = QR__R;
#> auto linkinv = [&](const arma::colvec& eta) mutable -> arma::colvec
#> {
#> return arma::exp(eta) / (1.0 + arma::exp(eta));
#> }
#> ;
#> auto mu_eta = [&](const arma::colvec& eta) mutable -> arma::colvec
#> {
#> return arma::exp(eta) / arma::pow((1.0 + arma::exp(eta)), 2.0);
#> }
#> ;
#> auto variance = [&](const arma::colvec& mu) mutable -> arma::colvec
#> {
#> return mu % (1.0 - mu);
#> }
#> ;
#> 
#> auto&& arma___tmp_00 = arma::linspace<arma::colvec>(1, maxit, maxit);
#> for (auto arma___tmp_11 = arma___tmp_00.begin();arma___tmp_11 != arma___tmp_00.end();++arma___tmp_11)
#> {
#> arma::colvec s_old = s;
#> arma::colvec mu = linkinv(eta);
#> arma::colvec mu_p = mu_eta(eta);
#> arma::colvec z = eta + (y - mu) / mu_p;
#> arma::colvec W = arma::pow(mu_p, 2.0) / variance(mu);
#> arma::mat C = arma::chol(arma::trans(Q) * (arma::diagmat(W) * Q));
#> arma::colvec s1 = arma::solve(arma::trimatl(arma::trans(C)), arma::trans(Q) * (W % z) );
#> s = arma::solve(arma::trimatu(C), s1 );
#> eta = Q * s;
#> 
#> if ( arma::accu(arma::sqrt(arma::trans(s - s_old) * (s - s_old))) < tol ) 
#> {
#> break;
#> 
#> }
#> 
#> }
#> 
#> 
#> arma::colvec beta = arma::solve(arma::trimatu(R), arma::trans(Q) * (eta) );
#> return beta;
#> }

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
#> [1] TRUE

# Note: glm.fit does more than simply computing the coefficients
system.time(glm_fit_cpp(X, y, 25, 1e-6))
#>    user  system elapsed 
#>   2.305   0.218   2.844
system.time(glm.fit(X, y, family = binomial()))
#>    user  system elapsed 
#>   4.767   1.004   6.424
