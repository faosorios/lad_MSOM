## loading packages and reading sources
library(india)
source("../code/MSOM.R")
source("../code/envelope.R")

## fitting linear model
data(aircraft)
fm <- lad(cost ~ ., data = aircraft)

## Table 4: removing selected observations
rm08 <- lad(cost ~ ., data = aircraft, subset = -8)
rm14 <- lad(cost ~ ., data = aircraft, subset = -14)
rm22 <- lad(cost ~ ., data = aircraft, subset = -22)

summary(fm)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  1.9510964363 12.682807120
#aspect      -3.0466304002  2.210406065
#lift2drag    1.4735937299  1.487921344
#weight       0.0022304441  0.000600496
#thrust      -0.0009646206  0.000625113

summary(rm08)$coef[,1:2]
#                Estimate     Std.Error
#(Intercept)  1.671983020 16.3236829105
#aspect      -4.749064250  2.7712211853
#lift2drag    2.652676768  1.9408435921
#weight       0.003265141  0.0008083801
#thrust      -0.001952660  0.0008281230

summary(rm14)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  2.4923401252 14.434034130
#aspect      -3.9593658135  2.627340459
#lift2drag    4.5386928750  4.430381944
#weight       0.0019543745  0.000734989
#thrust      -0.0009448061  0.000725396

summary(rm22)$coef[,1:2]
#                Estimate     Std.Error
#(Intercept)  9.376461426 13.3903723769
#aspect      -4.110943106  2.2711120945
#lift2drag    1.635141030  1.5405008579
#weight       0.002038673  0.0008230097
#thrust      -0.001067926  0.0007302342

c(fm$scale, rm08$scale, rm14$scale, rm22$scale) 
#[1] 7.361162 7.475893 7.450202 5.050452

c(fm$lambda, rm08$lambda, rm14$lambda, rm22$lambda)
#[1] 10.53901 13.19639 11.99349 10.77186

## computing MSOM diagnostic measures
z <- msom.lad(fm)
cutoff <- qchisq(0.95, df = 1)

## Fig 9(a): Likelihood ratio test statistic
par(pty = "s")
plot(z$LR, type = "h", ylab = "LR statistic", ylim = c(0,9), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(22, z$LR[22], as.character(22), pos = 3, cex = 1.2)

## Fig 9(b): Wald statistic
par(pty = "s")
plot(z$Wald, type = "h", ylab = "Wald statistic", ylim = c(0,9), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(22, z$Wald[22], as.character(22), pos = 3, cex = 1.2)

## gradient statistic (Fig 4(d), Supplementary material)
par(pty = "s")
plot(z$gradient, type = "h", ylab = "Gradient statistic", ylim = c(0,5), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(22, z$gradient[22], as.character(22), pos = 3, cex = 1.2)

## Fig 10(a): quantile residuals
par(pty = "s")
plot(z$qr, ylab = "quantile residuals", ylim = c(-2,4), lwd = 2, cex.lab = 1.3)
abline(h = 0, lwd = 2, col = "gray75")
abline(h = c(-2,2), lwd = 2, lty = 2, col = "red")
text(22, z$qr[22], as.character(22), pos = 3, cex = 1.2)

## Fig 10(b): standardized residuals
par(pty = "s")
plot(z$std, ylab = "standardized residuals", ylim = c(-2,4), lwd = 2, cex.lab = 1.3)
abline(h = 0, lwd = 2, col = "gray75")
abline(h = c(-2,2), lwd = 2, lty = 2, col = "red")
text(22, z$std[22], as.character(22), pos = 3, cex = 1.2)

## Fig 10(c): simulated envelope based on quantile residuals
z <- envelope(fm, reps = 3000)
ylim <- c(-3.2, 3.4)
par(pty = "s")
qqnorm(z$resid, ylim = ylim, lwd = 2, main = "", cex.lab = 1.3)
par(new = TRUE)
qqnorm(z$envel[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(z$envel[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")

## Fig 10(d): simulated envelope based on 'standarized' residuals
z <- envelope.std(fm, reps = 3000)
ylim <- c(-3.2, 3.4)
par(pty = "s")
qqnorm(z$resid, ylim = ylim, lwd = 2, main = "", cex.lab = 1.3)
par(new = TRUE)
qqnorm(z$envel[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(z$envel[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
