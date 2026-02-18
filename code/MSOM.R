## ID: MSOM.R, last updated 2026-01-28, F.Osorio

msom.lad <- function(model) 
{ ## Cook's distance for LAD regression
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

  # extracting elements estimates using the full data
  lambda <- obj$lambda
  basic <- obj$basic
  resid <- obj$resid

  # computing leverages
  rs <- svd(x, nv = 0)
  levs <- rowSums(rs$u^2)

  # deleted statistics
  cf <- rel <- matrix(0, nrow = n, ncol = p)
  SAD <- dev <- rep(0, n)
  for (i in 1:n) {
    z <- lad.fit.BR(x[-i,], y[-i])
    SAD[i] <- z$SAD
    dev[i] <- y[i] - sum(x[i,] * z$coef)
    cf[i,] <- z$coef
    rel[i,] <- (z$coef - obj$coef) / obj$coef
  }

  # computing test statistics
  LR <- 2 * (obj$SAD - SAD) / lambda
  Wald <- (1 - levs) * (dev / lambda)^2
  R <- 1 / (1 - levs)
  R[basic] <- 0 
  BF <- sign(resid) * dev / lambda

  # standardized residuals
  std <- sqrt(1 - levs) * dev / lambda
  qr  <- residuals(fm, type = "quantile")

  names(SAD) <- names(dev) <- as.character(1:n)
  rownames(cf) <- rownames(rel) <- as.character(1:n)
  colnames(cf) <- colnames(rel) <- colnames(obj$coef)
  names(LR) <- names(R) <- names(Wald) <- as.character(1:n)
  out <- list(SAD = SAD, dev = dev, coef = cf, rel = 100 * rel, levs = levs, LR = LR, Rao = R, Wald = Wald, gradient = BF, std = std, qr = qr)
  out
}
