## ID: envelope.R, last updated 2026-02-02, F.Osorio

envelope.std <- function(object, reps = 50, conf = 0.95, plot.it = TRUE)
{ ## simulated envelope for LAD regression
  envel <- function(n, x, mu, dispersion, lambda, reps, conf) {
    conf <- 1 - conf
    # initialize progress bar
    cat(" Progress:\n")
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    elims <- matrix(0, nrow = n, ncol = reps)
    for (i in 1:reps) {
      rsp <- rlaplace(n, location = mu, scale = dispersion)
      # computing 'studentized' residual
      res <- rstudent.lad(x, rsp, lambda) # rsp
      elims[,i] <- sort(res)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    band <- matrix(0, nrow = n, ncol = 2)
    for (i in 1:n)
      band[i,] <- quantile(elims[i,], probs = c(conf / 2, 1 - conf / 2))
    band
  }

  if (!inherits(object, "lad"))
    stop("Use only with 'lad' objects")
  if (object$method != "BR")
    stop("Use only with objects fitted by 'BR' option")
  x <- model.matrix(object$terms, object$model, object$contrast)
  y <- model.response(object$model)
  mu <- fitted(object)
  n <- object$dims[1]
  scale <- object$scale
  lambda <- object$lambda

  # selected statistic and envelope
  z <- rstudent.lad(x, y, lambda)
  band <- envel(n, x, mu, scale, lambda, reps, conf)

  if (plot.it) {
    ylim <- range(z, band)
    qqnorm(z, ylim = ylim, main = paste("Q-Q plot of studentized residuals"))
    par(new = TRUE)
    qqnorm(band[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
    par(new = TRUE)
    qqnorm(band[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
  }

  output <- list(residuals = z, envelope = band)
  invisible(output)
}

rstudent.lad <- function(x, y, lambda)
{ ## 'studentized' residual for lad
  n <- nrow(x)

  # computing leverages
  rs <- svd(x, nv = 0)
  levs <- rowSums(rs$u^2)

  dev <- rep(0, n)
  # computing deleted residuals
  for (i in 1:n) {
    cf <- lad.fit.BR(x[-i,], y[-i])$coef
    dev[i] <- y[i] - sum(x[i,] * cf)
  }

  # standardized residuals
  std <- sqrt(1 - levs) * dev / lambda
  std
}
