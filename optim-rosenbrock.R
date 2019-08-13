# Rosenbrock optimization problem
library(armacmp)

# using simulated annealing
eval_fun <- function(x) {
  return(100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2)
}
rosenbrock <- compile_optimization_problem(
  data = list(),
  evaluate = eval_fun,
  optimizer = optimizer_SA()
)

system.time(result <- rosenbrock(matrix(c(-23, 50), nrow = 1)))
#>    user  system elapsed 
#>   0.026   0.000   0.028
result
#>         [,1]     [,2]
#> [1,] 1.03939 1.080484

system.time(result <- stats::optim(c(-23, 50), eval_fun, method = "SANN"))
#>    user  system elapsed 
#>   0.033   0.002   0.037
result
#> $par
#> [1] -6.840219 46.826866
#> 
#> $value
#> [1] 61.61553
#> 
#> $counts
#> function gradient 
#>    10000       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL

# using the gradient and L-BFGS
gradient_fun <- function(x) {
  grad <- x # make a copy of x and use that as the gradient
  grad[1] <- -400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1])
  grad[2] <- 200 * (x[2] - x[1]^2)
  return(grad)
}
rosenbrock <- compile_optimization_problem(
  data = list(),
  evaluate = eval_fun,
  gradient = gradient_fun,
  optimizer = optimizer_L_BFGS()
)

system.time(result <- rosenbrock(matrix(c(-23, 50), nrow = 1)))
#>    user  system elapsed 
#>       0       0       0
result
#>           [,1]      [,2]
#> [1,] 0.9999997 0.9999995

# now with optim
system.time(result <- stats::optim(c(-23, 50), eval_fun, gradient_fun, method = "L-BFGS-B"))
#>    user  system elapsed 
#>   0.014   0.000   0.015
result
#> $par
#> [1] 0.9999999 0.9999999
#> 
#> $value
#> [1] 3.403617e-15
#> 
#> $counts
#> function gradient 
#>      112      112 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
