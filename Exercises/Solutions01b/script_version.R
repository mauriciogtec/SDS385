library(tidyverse)

# Get the Data -----------------------------------------------------
wdbc <- read_csv(
  file = url("https://raw.githubusercontent.com/jgscott/SDS385/master/data/wdbc.csv"),
  col_names = FALSE
)
X <- wdbc[ ,3:10] %>% # select the first columns
  scale() %>% # subtract mean and divide by standard deviation for each column
  cbind(1, .) %>% # add columns of one 
  as.matrix()
y <- (wdbc[[2]] == "M") %>% 
  as.numeric() # conver to 0's and 1's

# Define the likelihood ----------------------------------------------

#' @title Sigmoid function
#' @params
#'   u a numeric vector
#' @value Applies the transformation ui -> (1 / (1 + exp(-ui))) to each entry ui of u
sigmoid <- function(u) {
  # Data type validation
  stopifnot(is.numeric(u))
  # OutputS
  1 / (1 + exp(-u))
}

#' @title Logloss function
#' @params
#'   beta a numeric vector of coefficients
#'   X a matrix with rows as observations and columns as numeric variables
#'   y a logical vector of responses
#' @value Computes the negative of the loglikelihood
logloss <- function(beta, X, y) {
  # Data type validation
  stopifnot(
    is.numeric(beta),
    inherits(X, c("matrix", "Matrix")),
    is.numeric(y) # Should I require logical here?
  )
  # Compue predictions
  reg <- X %*% beta 
  w <- sigmoid(as.numeric(reg))
  # Logloss
  ll <- -sum(log(1 - w)) - crossprod(y, reg)
  # Output likelihood
  drop(ll)
}

#' @title Logloss gradient function
#' @params
#'   beta a numeric vector of coefficients
#'   X a matrix with rows as observations and columns as numeric variables
#'   y a logical vector of responses
#' @value Computes the gradient of the negative of the loglikelihood. Returns a numeric vector
logloss_grad <- function(beta, X, y) {
  # Data type validation
  stopifnot(
    is.numeric(beta),
    inherits(X, c("matrix", "Matrix")),
    is.numeric(y) # Should I require logical here?
  )
  # Compute predictions
  reg <- X %*% beta 
  w <- sigmoid(drop(reg))
  # Output gradient
  drop(crossprod(X, w - y))
}

# LineSearch -------------------------------------
#' @title Linesearch with backtracking
#' @params
#'   x the current point at which the function will be evaluated
#'   direction a search direction in which to search  
#'   f the function f we are trying to minimize
#'   grad the gradient of f
#'   control additional options to control the line serach algorithm (see details)
#' @value the step size that satistfies the armijo condition
bt_linesearch <- function(x, direction, f, grad, control = NULL) {
  # Data type validation
  stopifnot(
    is.numeric(x),
    is.numeric(direction) && length(direction) == length(x),
    is.function(f) && is.function(grad),
    is.list(control) || is.null(control)
  )
  # Default values and user controls
  bt_rate <- 0.9
  min_decrease_const <- 0.1
  init_step_size <- 1.0
  maxit <- 1e4
  for (i in seq_along(names(control))) { 
    assign(names(control)[i], control[[i]]) 
  }
  # Compute minimum dicrease for the Armijo condiion
  graddir_prod <- crossprod(grad(x), direction)
  if (graddir_prod >= 0) {
    warning("Not a valid descent direction, returning step_size zero")
    return(0)
  }
  sufficient_decrease <- f(x) - min_decrease_const * step_size * graddir_prod
  # Backtrack to find optimal step size
  i  <- 1
  step_size <- init_step_size
  repeat {
    x_new <- x + step_size * direction
    f_new <- f(x_new) 
    if (f_new < sufficient_decrease) {
      break
    } 
    if (i == maxit) {
      warning("Max number of iterations in linesearch, returning step size zero")
      step_size <- 0
      break
    }
    step_size <- step_size * bt_rate
    i <- i + 1
  }
  # Output
  list(
    x_new = x_new,
    f_new = f_new,
    step_size = step_size,
    iter = i
  )
}

# Optimization -----------------------------------------------------
