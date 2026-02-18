## loading packages and reading sources
library(india)
source("../code/MSOM.R")
source("../code/envelope.R")

## fitting linear model
fm <- lad(stack.loss ~ ., data = stackloss)

## Table 1: removing selected observations
rm02 <- lad(stack.loss ~ ., data = stackloss, subset = -2)
rm03 <- lad(stack.loss ~ ., data = stackloss, subset = -3)
rm04 <- lad(stack.loss ~ ., data = stackloss, subset = -4)
rm17 <- lad(stack.loss ~ ., data = stackloss, subset = -17)
rm21 <- lad(stack.loss ~ ., data = stackloss, subset = -21)

summary(fm)$coef[,1:2]
#                Estimate  Std.Error
#(Intercept) -39.68985507 13.0314316
#Air.Flow      0.83188406  0.1477300
#Water.Temp    0.57391304  0.4031510
#Acid.Conc.   -0.06086957  0.1712118

summary(rm02)$coef[,1:2]
#                Estimate  Std.Error
#(Intercept) -39.69396552 12.6316608
#Air.Flow      0.82974138  0.1519089
#Water.Temp    0.57758621  0.3921722
#Acid.Conc.   -0.06034483  0.1704040

summary(rm03)$coef[,1:2]
#                Estimate  Std.Error
#(Intercept) -39.65188470 10.8313433
#Air.Flow      0.83037694  0.1254675
#Water.Temp    0.58093126  0.3328416
#Acid.Conc.   -0.06208426  0.1414111

summary(rm04)$coef[,1:2]
#                Estimate  Std.Error
#(Intercept) -39.98644986 10.6543713
#Air.Flow      0.83468835  0.1230113
#Water.Temp    0.56368564  0.3436712
#Acid.Conc.   -0.05691057  0.1397474

summary(rm17)$coef[,1:2]
#               Estimate  Std.Error
#(Intercept) -36.4245283 14.9724806
#Air.Flow      0.7990566  0.1344800
#Water.Temp    0.6962264  0.3688379
#Acid.Conc.   -0.1056604  0.1907009

summary(rm21)$coef[,1:2]
#                Estimate  Std.Error
#(Intercept) -39.98644986 11.8448327
#Air.Flow      0.83468835  0.1483138
#Water.Temp    0.56368564  0.4056148
#Acid.Conc.   -0.05691057  0.1554192

c(fm$scale, rm02$scale, rm03$scale, rm04$scale, rm17$scale, rm21$scale) 
#[1] 2.833893 2.975182 2.591366 2.435015 2.942298 2.302217

c(fm$lambda, rm02$lambda, rm03$lambda, rm04$lambda, rm17$lambda, rm21$lambda) 
#[1] 3.552933 3.443342 2.933291 2.899740 3.232495 3.206190

## computing MSOM diagnostic measures
z <- msom.lad(fm)
nobs <- fm$dims[1]
cutoff <- qchisq(0.95, df = 1)

## Fig 1(a): Likelihood ratio test statistic
obs <- (1:nobs)[z$LR > cutoff]
par(pty = "s")
plot(z$LR, type = "h", ylab = "LR statistic", ylim = c(0,6), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(obs, z$LR[obs], as.character(obs), pos = 3, cex = 1.2)

## Fig 1(b): Wald statistic
obs <- (1:nobs)[z$Wald > cutoff]
par(pty = "s")
plot(z$Wald, type = "h", ylab = "Wald statistic", ylim = c(0,6), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(obs, z$Wald[obs], as.character(obs), pos = 3, cex = 1.2)

## gradient statistic (Fig 4(a), Supplementary material)
par(pty = "s")
plot(z$gradient, type = "h", ylab = "Gradient statistic", ylim = c(0,3), lwd = 2, cex.lab = 1.3)
text(obs, z$gradient[obs], as.character(obs), pos = 3, cex = 1.2)

## Fig 2(a): quantile residuals
par(pty = "s")
plot(z$qr, ylim = c(-3,3), ylab = "quantile residuals", lwd = 2, cex.lab = 1.3)
abline(h = 0, lwd = 2, col = "gray75")
abline(h = c(-2,2), lwd = 2, lty = 2, col = "red")
text(4, z$qr[4], as.character(4), pos = 3, cex = 1.2)
text(21, z$qr[21], as.character(21), pos = 1, cex = 1.2)

## Fig 2(b): standardized residuals
par(pty = "s")
plot(z$std, ylim = c(-3,3), ylab = "standardized residuals", lwd = 2, cex.lab = 1.3)
abline(h = 0, lwd = 2, col = "gray75")
abline(h = c(-2,2), lwd = 2, lty = 2, col = "red")
text(4, z$std[4], as.character(4), pos = 3, cex = 1.2)
text(21, z$std[21], as.character(21), pos = 1, cex = 1.2)

## Fig 2(c): simulated envelope based on quantile residuals
z <- envelope(fm, reps = 3000)
ylim <- c(-3.2, 3.4)
par(pty = "s")
qqnorm(z$resid, ylim = ylim, lwd = 2, main = "", cex.lab = 1.3)
par(new = TRUE)
qqnorm(z$envel[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(z$envel[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")

## Fig 2(d): simulated envelope based on 'standarized' residuals
z <- envelope.std(fm, reps = 3000)
ylim <- c(-3.2, 3.4)
par(pty = "s")
qqnorm(z$resid, ylim = ylim, lwd = 2, main = "", cex.lab = 1.3)
par(new = TRUE)
qqnorm(z$envel[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(z$envel[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
