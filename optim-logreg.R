library(armacmp)
optimize_lbfgs <- compile_optimization_problem(
  data = list(design_matrix = type_matrix(), response = type_colvec()),
  evaluate = function(beta) {
    return(norm(response - design_matrix %*% beta)^2)
  },
  gradient = function(beta) {
    return(-2 %*% t(design_matrix) %*% (response - design_matrix %*% beta))
  },
  optimizer = optimizer_L_BFGS()
)

# this example is taken from the RcppEnsmallen package
# https://github.com/coatless/rcppensmallen/blob/master/src/example-linear-regression-lbfgs.cpp
n <- 1e6
beta <- c(-2, 1.5, 3, 8.2, 6.6)
p <- length(beta)
X <- cbind(1, matrix(rnorm(n), ncol = p - 1))
y <- X %*% beta + rnorm(n / (p - 1))

# Run optimization with lbfgs fullly in C++
optimize_lbfgs(
  design_matrix = X,
  response = y,
  beta = matrix(runif(p), ncol = 1)
)
#>           [,1]
#> [1,] -2.002257
#> [2,]  1.500687
#> [3,]  3.000419
#> [4,]  8.200791
#> [5,]  6.599552
