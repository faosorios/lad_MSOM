## loading packages and reading sources
library(india)
source("../code/MSOM.R")

## fitting linear models and computing gradient statistic 
fm <- lad(stack.loss ~ ., data = stackloss)
z1 <- msom.lad(fm)$gradient

data(skeena)
skeena$ratio <- skeena$recruits / skeena$spawners
fm <- lad(log(ratio) ~ spawners, data = skeena)
z2 <- msom.lad(fm)$gradient

data(ereturns)
fm <- lad(m.marietta ~ CRSP, data = ereturns)
z3 <- msom.lad(fm)$gradient

data(aircraft)
fm <- lad(cost ~ ., data = aircraft)
z4 <- msom.lad(fm)$gradient

## Fig 5(a)
obs <- c(4,21)
par(pty = "s")
plot(z1, type = "h", ylab = "Gradient statistic", ylim = c(0,3), lwd = 2, cex.lab = 1.3)
text(obs, z1[obs], as.character(obs), pos = 3, cex = 1.2)

## Fig 5(b)
par(pty = "s")
plot(z2, type = "h", ylab = "Gradient statistic", ylim = c(0,2.5), lwd = 2, cex.lab = 1.3)
text(12, z2[12], as.character(12), pos = 3, cex = 1.2)

## Fig 5(c)
cutoff <- qchisq(0.95, df = 1)
par(pty = "s")
plot(z3, type = "h", ylab = "Gradient statistic", ylim = c(0,6), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(8, z3[8], as.character(8), pos = 3, cex = 1.2)

## Fig 5(d)
par(pty = "s")
plot(z4, type = "h", ylab = "Gradient statistic", ylim = c(0,5), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(22, z4[22], as.character(22), pos = 3, cex = 1.2)
