## loading packages and reading sources
library(india)
source("../code/Recruits.R")
source("../code/MSOM.R")
source("../code/envelope.R")

## fitting linear model
data(skeena)
skeena$ratio <- skeena$recruits / skeena$spawners
fm <- lad(log(ratio) ~ spawners, data = skeena)

## get info for Fig 3
x <- skeena$spawners
y <- skeena$recruits
id <- fm$basic
cf <- fm$coef

## Fig 3
par(pty = "s", mai = c(1,1,.35,.35))
plot(recruits ~ spawners, data = skeena, xlab = "Spawners", ylab = "Recruits", xlim = c(80,1100), ylim = c(0,3500), lwd = 2, cex.lab = 1.3)
par(new = TRUE)
plot(x[id], y[id], xlim = c(80,1100), ylim = c(0,3500), xlab = "", ylab = "", lwd = 2, pch = 21, bg = "red", axes = FALSE)
xgrid <- seq(70, 1100, length = 300)
ygrid <- Recruits(xgrid, cf)
lines(xgrid, ygrid, lwd = 2, lty = 2, col = "red")
obs <- c(1,4,5,9,12,16,18,20)
text(x[obs], y[obs], as.character(obs), pos =3, cex = 1.2)

## deleting case 12
rm12 <- lad(log(ratio) ~ spawners, data = skeena, subset = -12)
cf12 <- rm12$coef
ygrid12 <- Recruits(xgrid, cf12)
lines(xgrid, ygrid12, lwd = 2, lty = 4, col = "gray35")

## Table 2: removing selected observations
rm01 <- lad(log(ratio) ~ spawners, data = skeena, subset = -1)
rm04 <- lad(log(ratio) ~ spawners, data = skeena, subset = -4)
rm05 <- lad(log(ratio) ~ spawners, data = skeena, subset = -5)
rm09 <- lad(log(ratio) ~ spawners, data = skeena, subset = -9)
rm16 <- lad(log(ratio) ~ spawners, data = skeena, subset = -16)
rm18 <- lad(log(ratio) ~ spawners, data = skeena, subset = -18)
rm20 <- lad(log(ratio) ~ spawners, data = skeena, subset = -20)

summary(fm)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  1.1273270127 0.2882327717
#spawners    -0.0004569882 0.0004687638

summary(rm01)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  1.2501384155 0.3418382808
#spawners    -0.0007051402 0.0005715341
 
summary(rm04)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  1.2385782710 0.3415251805
#spawners    -0.0006563632 0.0005473425
 
summary(rm05)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  1.2501384155 0.3368656878
#spawners    -0.0007051402 0.0005561148
 
summary(rm09)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  1.2501384155 0.3481236128
#spawners    -0.0007051402 0.0005884437

summary(rm12)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  1.2385782710 0.3486655327
#spawners    -0.0006563632 0.0005576464
 
summary(rm16)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  1.0807863230 0.3396291201
#spawners    -0.0003735819 0.0005425925
 
summary(rm18)$coef[,1:2]
#                 Estimate   Std.Error
#(Intercept)  1.1273270127 0.298677191
#spawners    -0.0004569882 0.000481584
 
summary(rm20)$coef[,1:2]
#                 Estimate    Std.Error
#(Intercept)  1.2501384155 0.3444849256
#spawners    -0.0007051402 0.0005675321

c(fm$scale, rm01$scale, rm04$scale, rm05$scale, rm09$scale, rm12$scale, rm16$scale, rm18$scale, rm20$scale) 
#[1] 0.5082636 0.5157247 0.4984088 0.4955790 0.5162590 0.4530636 
#[7] 0.5088729 0.4953226 0.5244354

c(fm$lambda, rm01$lambda, rm04$lambda, rm05$lambda, rm09$lambda, rm12$lambda, rm16$lambda, rm18$lambda, rm20$lambda) 
#[1] 0.6454322 0.7506382 0.7366803 0.7506382 0.7506382 0.7366803 
#[7] 0.7004734 0.6609216 0.7688171

## computing MSOM diagnostic measures
z <- msom.lad(fm)
cutoff <- qchisq(0.95, df = 1)

## Fig 4(a): Likelihood ratio test statistic
par(pty = "s")
plot(z$LR, type = "h", ylab = "LR statistic", ylim = c(0,5), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(12, z$LR[12], as.character(12), pos = 3, cex = 1.2)

## Fig 4(b): Wald statistic
par(pty = "s")
plot(z$Wald, type = "h", ylab = "Wald statistic", ylim = c(0,5), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(12, z$Wald[12], as.character(12), pos = 3, cex = 1.2)

## gradient statistic (Fig 4(b), Supplementary material)
par(pty = "s")
plot(z$gradient, type = "h", ylab = "Gradient statistic", ylim = c(0,2.5), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(12, z$gradient[12], as.character(12), pos = 3, cex = 1.2)

## Fig 5(a): quantile residuals
par(pty = "s")
plot(z$qr, ylab = "quantile residuals", ylim = c(-3,2), lwd = 2, cex.lab = 1.3)
abline(h = 0, lwd = 2, col = "gray75")
abline(h = c(-2,2), lwd = 2, lty = 2, col = "red")
text(12, z$qr[12], as.character(12), pos = 1, cex = 1.2)

## Fig 5(b): standardized residuals
par(pty = "s")
plot(z$std, ylab = "standardized residuals", ylim = c(-3,2), lwd = 2, cex.lab = 1.3)
abline(h = 0, lwd = 2, col = "gray75")
abline(h = c(-2,2), lwd = 2, lty = 2, col = "red")
text(12, z$std[12], as.character(12), pos = 1, cex = 1.2)

## Fig 5(c): simulated envelope based on quantile residuals
z <- envelope(fm, reps = 3000)
ylim <- c(-3.5, 3.4)
par(pty = "s")
qqnorm(z$resid, ylim = ylim, lwd = 2, main = "", cex.lab = 1.3)
par(new = TRUE)
qqnorm(z$envel[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(z$envel[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")

## Fig 5(d): simulated envelope based on 'standarized' residuals
z <- envelope.std(fm, reps = 3000)
ylim <- c(-3.5, 3.4)
par(pty = "s")
qqnorm(z$resid, ylim = ylim, lwd = 2, main = "", cex.lab = 1.3)
par(new = TRUE)
qqnorm(z$envel[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(z$envel[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
