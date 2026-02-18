## ID: Recruits.R, last updated 2026-02-18, F.Osorio

Recruits <- function(x, coef) 
{ # nonlinear mean function (Ricker, 1954)
  theta <- c(exp(coef[1]), -coef[2])
  y <- theta[1] * x * exp(-theta[2] * x)
  y
}

lin2orig <- function(coef)
{ # linearized estimates to original parameterization
  theta <- c(exp(coef[1]), -coef[2])
  theta
}
