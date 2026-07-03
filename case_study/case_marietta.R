## loading packages and reading sources
library(india) # also attach 'fastmatrix' and 'L1pack'
source("../code/MSOM.R")
source("../code/envelope.R")

## fitting linear model
data(ereturns)
fm <- lad(m.marietta ~ CRSP, data = ereturns)

## get info for Fig.1
x <- ereturns$CRSP
y <- ereturns$m.marietta
id <- fm$basic
cf <- fm$coef

## Fig.1
par(pty = "s")
plot(x, y, xlab = "CRSP", ylab = "Martin Marietta", xlim = c(-.1, .12), ylim = c(-.15,.8), lwd = 2, cex.lab = 1.3)
par(new = TRUE)
plot(x[id], y[id], xlim = c(-.1, .12), ylim = c(-.15,.8), xlab = "", ylab = "", lwd = 2, pch = 21, bg = "red", axes = FALSE)
abline(cf, lwd = 2, lty = 2, col = "red")
obs <- c(8,15,34)
text(x[obs], y[obs], as.character(obs), pos = 3, cex = 1.2)
text(x[58], y[58], as.character(58), pos = 1, cex = 1.2)

## Table 1: removing selected observations
rm08 <- lad(m.marietta ~ CRSP, data = ereturns, subset = -8)
rm15 <- lad(m.marietta ~ CRSP, data = ereturns, subset = -15)
rm34 <- lad(m.marietta ~ CRSP, data = ereturns, subset = -34)
rm58 <- lad(m.marietta ~ CRSP, data = ereturns, subset = -58)
rm00 <- lad(m.marietta ~ CRSP, data = ereturns, subset = -c(8,15,34)) # remove obs. 8, 15 and 34

summary(fm)$coef[,1:2]
#                Estimate  Std.Error
#(Intercept) -0.003761392 0.01345317
#CRSP         1.252373418 0.31323871
 
summary(rm08)$coef[,1:2]
#               Estimate  Std.Error
#(Intercept) -0.01202718 0.01326754
#CRSP         1.11083591 0.32747437
 
summary(rm15)$coef[,1:2]
#                Estimate  Std.Error
#(Intercept) -0.003761392 0.01346953
#CRSP         1.252373418 0.31201685
 
summary(rm34)$coef[,1:2]
#                Estimate  Std.Error
#(Intercept) -0.003761392 0.01350527
#CRSP         1.252373418 0.31188365
 
summary(rm58)$coef[,1:2]
#                Estimate  Std.Error
#(Intercept) -0.003761392 0.01516169
#CRSP         1.252373418 0.35378253

summary(rm00)$coef[,1:2]
#               Estimate  Std.Error
#(Intercept) -0.01202718 0.01464172
#CRSP         1.11083591 0.35663212

c(fm$scale, rm08$scale, rm15$scale, rm34$scale, rm58$scale, rm00$scale) 
#[1] 0.08363633 0.07174567 0.08021267 0.08008132 0.08157402 0.06362382

c(fm$lambda, rm08$lambda, rm15$lambda, rm34$lambda, rm58$lambda, rm00$lambda) 
#[1] 0.1020079 0.1004063 0.1014499 0.1014499 0.1143539 0.1090054

## Table 2: percentage of relative change
rel <- matrix(0, nrow = 5, ncol = 4)
rel[1,1:2] <- (rm08$coef - fm$coef) / fm$coef
rel[2,1:2] <- (rm15$coef - fm$coef) / fm$coef
rel[3,1:2] <- (rm34$coef - fm$coef) / fm$coef
rel[4,1:2] <- (rm58$coef - fm$coef) / fm$coef
rel[5,1:2] <- (rm00$coef - fm$coef) / fm$coef
rel[,3] <- (c(rm08$scale, rm15$scale, rm34$scale, rm58$scale, rm00$scale) - fm$scale) / fm$scale
rel[,4] <- (c(rm08$lambda, rm15$lambda, rm34$lambda, rm58$lambda, rm00$lambda) - fm$lambda) / fm$lambda
rel <- 100 * rel
rnames <- c(names(fm$coef), "scale", "lambda")
colnames(rel) <- rnames
rownames(rel) <- c("8","15","34","58","all")

rel
#    (Intercept)      CRSP      scale     lambda
#8      219.7535 -11.30154 -14.217102 -1.5700369
#15       0.0000   0.00000  -4.093509 -0.5470343
#34       0.0000   0.00000  -4.250556 -0.5470343
#58       0.0000   0.00000  -2.465801 12.1030044
#all    219.7535 -11.30154 -23.928011  6.8598166

## computing MSOM diagnostic measures
z <- msom.lad(fm)

## some definitions
nobs <- fm$dims[1]
cutoff <- qchisq(0.95, df = 1)

## Fig 2(a): Likelihood ratio test statistic
obs <- (1:nobs)[z$LR > cutoff]
par(pty = "s")
plot(z$LR, type = "h", ylab = "LR statistic", ylim = c(0,12), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(obs, z$LR[obs], as.character(obs), pos = 3, cex = 1.2)

## Fig 2(b): Wald statistic
par(pty = "s")
plot(z$Wald, type = "h", ylab = "Wald statistic", ylim = c(0,30), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(obs, z$Wald[obs], as.character(obs), pos = 3, cex = 1.2)

## Fig 3(a): quantile residuals
obs <- c(8,15,34)
par(pty = "s")
plot(z$qr, ylim = c(-2,6), ylab = "quantile residuals", lwd = 2, cex.lab = 1.3)
abline(h = 0, lwd = 2, col = "gray75")
abline(h = c(-2,2), lwd = 2, lty = 2, col = "red")
text(obs, z$qr[obs], as.character(obs), pos = 3, cex = 1.2)

## Fig 3(b): standardized residuals
obs <- c(8,34)
par(pty = "s")
plot(z$std, ylim = c(-2,6), ylab = "standardized residuals", lwd = 2, cex.lab = 1.3)
abline(h = 0, lwd = 2, col = "gray75")
abline(h = c(-2,2), lwd = 2, lty = 2, col = "red")
text(obs, z$std[obs], as.character(obs), pos = 3, cex = 1.2)

## Fig 3(c): simulated envelope based on quantile residuals
z <- envelope(fm, reps = 3000)
ylim <- c(-4.1, 5.3)
par(pty = "s")
qqnorm(z$resid, ylim = ylim, lwd = 2, main = "", cex.lab = 1.3)
par(new = TRUE)
qqnorm(z$envel[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(z$envel[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")

## Fig 3(d): simulated envelope based on 'standarized' residuals
z <- envelope.std(fm, reps = 3000)
ylim <- c(-4.1, 5.3)
par(pty = "s")
qqnorm(z$resid, ylim = ylim, lwd = 2, main = "", cex.lab = 1.3)
par(new = TRUE)
qqnorm(z$envel[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(z$envel[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
