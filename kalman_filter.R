# The purpose of this script is to build a kalman filter. 
library(KFAS)
# old ---------------------------------------------------------------------


data("alcohol")

deaths <- window(alcohol[, 2], end = 2007)
population <- window(alcohol[, 6], end = 2007)

Zt <- matrix(c(1, 0), 1, 2)
Ht <- matrix(NA)
Tt <- matrix(c(1, 0, 1, 1), 2, 2)
Rt <- matrix(c(1, 0), 2, 1)
Qt <- matrix(NA)
a1 <- matrix(c(1, 0), 2, 1)
P1 <- matrix(0, 2, 2)
P1inf <- diag(2)

input_data <- deaths / population
# Insert some missing values :D 
input_data[c(10, 20, 30)] <- NA


model_gaussian <- SSModel(input_data ~ -1 + 
														SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, 
																			P1 = P1, P1inf = P1inf,),
													H = Ht)

plot(input_data, type="o")

fit_gaussian <- fitSSM(model_gaussian, inits = c(0, 0), method = "BFGS")
out_gaussian <- KFS(fit_gaussian$model)


plot.ts(deaths / population, type="o")

# Smoothed estimates of state
a_smoth <- out_gaussian$alphahat

# Filtered estimates of state
a_fil <- out_gaussian$att

lines(a_fil[, 1] + a_fil[, 2], col="blue")
lines(a_smoth[, 1] + a_smoth[, 2], col="red")

library(astsa)
library(xts)
library(fGarch)

djiar = diff(log(djia$Close))[-1]
plot(djiar, main="DJIA Returns", type="n")
lines(djiar)

djia.g <- garchFit(~garch(1,1), data=djiar,cond.dist='std')

plot(djia.g)

source("stocks_xts.Rda.Rda")


str(tbl)

test_stock <- as.data.frame(diff(log(tbl[,1])))


num <- nrow(test_stock)

A <- array(rep(1, 2), dim = c(1, 2, num))
mu0 <- mean(test_stock[,1])
sigma0 <- 1

Linn <- function(para) {
	delta <- para[1]
	beta <- para[2]
	Phi <- diag(c(delta, beta))
	kf <- Kfilter1(num, test_stock, A, mu0, sigma0, Phi, 0, 0, 1, 1, 0)
}

Linn(c(1, 1))




load("stocks_raw.Rda")
source("Kfilter1.R")

s1 <- cbind(tbl[,1])
s1 <- diff(log(s1))

s2 <- log(s1^2)

s2[which(!is.finite(s2))] <- 0

num <- nrow(s2)
A <- array(rep(1, 2), dim = c(1,  2, num))
mu0 <- matrix(0.1, ncol = 1, nrow = 2)
sigma0 <- diag(2)
Q <- diag(2)
R <- diag(1)

Linn <- function(para) {
	delta <- para[1]
	beta <- para[2]
	phi <- diag(c(delta, beta))
	kf = Kfilter1(num, s2, A, mu0, sigma0, phi, 0, 0, Q, R, 0)
	return(kf$like)
}
# Small sanity check
out <- Kfilter1(num, s2, A, mu0, sigma0, diag(2), 0, 0, Q, R, 0)


init.par <- c(.1, .1)
optim(init.par, Linn, NULL, method = c("BFGS"), hessian = TRUE)
optim(init.par, Linn)

Linn(init.par)
#			control = list(trace=1, REPORT=1))
