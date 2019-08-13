reprex::reprex({
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
  result

  system.time(result <- stats::optim(c(-23, 50), eval_fun, method = "SANN"))
  result

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
  result

  # now with optim
  system.time(result <- stats::optim(c(-23, 50), eval_fun, gradient_fun, method = "L-BFGS-B"))
  result
}, outfile = "optim-rosenbrock", venue = "R", show = FALSE)
