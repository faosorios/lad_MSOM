## experiment presented in Section 5

## loading packages and reading sources
library(india)
source("../code/iterations.R")

## fitting linear models (see Section 4 of the manuscript and Section 1 from 
## supplement)
data(ereturns)
fm <- lad(m.marietta ~ CRSP, data = ereturns)
marietta.IT <- iter.LAD(fm)

data(skeena)
skeena$ratio <- skeena$recruits / skeena$spawners
fm <- lad(log(ratio) ~ spawners, data = skeena)
skeena.IT <- iter.LAD(fm)

fm <- lad(stack.loss ~ ., data = stackloss)
stackloss.IT <- iter.LAD(fm)

data(aircraft)
fm <- lad(cost ~ ., data = aircraft)
aircraft.IT <- iter.LAD(fm)

## Elapsed total time of BR and IRLS algorithms
speed <- matrix(0, nrow = 8, ncol = 6)
speed[1,2:6] <- marietta.IT$BR$speed
speed[2,2:6] <- marietta.IT$EM$speed
speed[3,2:6] <- skeena.IT$BR$speed
speed[4,2:6] <- skeena.IT$EM$speed
speed[5,2:6] <- stackloss.IT$BR$speed
speed[6,2:6] <- stackloss.IT$EM$speed
speed[7,2:6] <- aircraft.IT$BR$speed
speed[8,2:6] <- aircraft.IT$EM$speed
tnames <- c("method",names(aircraft.IT$EM$speed))
colnames(speed) <- tnames
rownames(speed) <- c("marietta","marietta ","skeena","skeena ","stackloss","stackloss ","aircraft","aircraft ")
speed <- as.data.frame(speed)
speed$method <- rep(c("BR","IRLS"), 4)

# output (times may vary slightly). "elapsed" time is reported in Table 5
speed[,c(1,2,4)]
#           method user.self elapsed
#marietta       BR     0.003   0.003
#marietta     IRLS     0.018   0.018
#skeena         BR     0.006   0.006
#skeena       IRLS     0.029   0.029
#stackloss      BR     0.001   0.001
#stackloss    IRLS     0.007   0.007
#aircraft       BR     0.002   0.004
#aircraft     IRLS     0.014   0.014

## Elapsed total time of BR and IRLS algorithms
ITER <- matrix(0, nrow = 8, ncol = 5)
ITER[1,c(2,3,5)] <- quantile(marietta.IT$BR$iter, probs = c(0, .5, 1))
ITER[1,4] <- mean(marietta.IT$BR$iter)
ITER[2,c(2,3,5)] <- quantile(marietta.IT$EM$iter, probs = c(0, .5, 1))
ITER[2,4] <- mean(marietta.IT$EM$iter)
ITER[3,c(2,3,5)] <- quantile(skeena.IT$BR$iter, probs = c(0, .5, 1))
ITER[3,4] <- mean(skeena.IT$BR$iter)
ITER[4,c(2,3,5)] <- quantile(skeena.IT$EM$iter, probs = c(0, .5, 1))
ITER[4,4] <- mean(skeena.IT$EM$iter)
ITER[5,c(2,3,5)] <- quantile(stackloss.IT$BR$iter, probs = c(0, .5, 1))
ITER[5,4] <- mean(stackloss.IT$BR$iter)
ITER[6,c(2,3,5)] <- quantile(stackloss.IT$EM$iter, probs = c(0, .5, 1))
ITER[6,4] <- mean(stackloss.IT$EM$iter)
ITER[7,c(2,3,5)] <- quantile(aircraft.IT$BR$iter, probs = c(0, .5, 1))
ITER[7,4] <- mean(aircraft.IT$BR$iter)
ITER[8,c(2,3,5)] <- quantile(aircraft.IT$EM$iter, probs = c(0, .5, 1))
ITER[8,4] <- mean(aircraft.IT$EM$iter)
colnames(ITER) <- c("method","min","median","average","max")
rownames(ITER) <- c("marietta","marietta ","skeena","skeena ","stackloss","stackloss ","aircraft","aircraft ")
ITER <- as.data.frame(ITER)
ITER$method <- rep(c("BR","IRLS"), 4)

ITER
#           method min median    average max
#marietta       BR   2    2.0   2.583333   5
#marietta     IRLS  17   43.0  72.633333 500
#skeena         BR   2    3.5   3.821429   6
#skeena       IRLS  12   36.5  58.964286 267
#stackloss      BR   5    6.0   6.523810  10
#stackloss    IRLS  17   56.0  64.761905 190
#aircraft       BR   6    7.0   7.608696  11
#aircraft     IRLS  55  161.0 159.260870 419

## Fig.7
BR <- list(marietta = log(marietta.IT$BR$iter), skeena = log(skeena.IT$BR$iter), stackloss = log(stackloss.IT$BR$iter), aircraft = log(aircraft.IT$BR$iter))
EM <- list(marietta = log(marietta.IT$EM$iter), skeena = log(skeena.IT$EM$iter), stackloss = log(stackloss.IT$EM$iter), aircraft = log(aircraft.IT$EM$iter))
bnames <- c("Martin-Marietta","Skeena","Stackloss","Aircraft")
xticks <- c(3,8,20,55,150,400)
ticks <- log(xticks)

par(pty = "s")
boxplot(BR, ylim = c(.65, 6.25), ylab = "Number of iterations", names = bnames, boxwex = 0.6, lwd = 2, cex.lab = 1.3, yaxt = "n")
axis(side = 2, at = ticks, labels = as.character(xticks))

par(pty = "s")
boxplot(EM, ylim = c(.65, 6.25), ylab = "Number of iterations", names = bnames, boxwex = 0.6, lwd = 2, cex.lab = 1.3, yaxt = "n")
axis(side = 2, at = ticks, labels = as.character(xticks))
