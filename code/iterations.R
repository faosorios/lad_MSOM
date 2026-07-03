## ID: iterations.R, last updated 2026-06-28, F.Osorio

iter.LAD <- function(model) 
{ ## computing iterations 
  if (!inherits(model, "lad"))
    stop("Use only with 'lad' objects")
  obj <- model
  y <- model.response(obj$model, "numeric")
  x <- model.matrix(obj$terms, obj$model, obj$contrast)
  n <- obj$dims[1]
  p <- obj$dims[2]

  # checking
  if (obj$method != "BR")
    stop("Use only with objects fitted by 'BR' option")

  # BR estimation
  ITER <- rep(0, n)
  now <- proc.time()
  for (i in 1:n) {
    z <- lad.fit.BR(x[-i,], y[-i])
    ITER[i] <- z$iterations
  }
  BR.it <- ITER
  BR.speed <- proc.time() - now

  # EM estimation
  ITER <- rep(0, n)
  now <- proc.time()
  for (i in 1:n) {
    z <- lad.fit.EM(x[-i,], y[-i], maxiter = 500)
    ITER[i] <- z$iterations
  }
  EM.it <- ITER
  EM.speed <- proc.time() - now

  out <- list(BR = list(iterations = BR.it, speed = BR.speed), EM = list(iterations = EM.it, speed = EM.speed))
  out
}
