# Simulate garch process
library(TSA)
# Alpha = c(intercept, arch coef, arch coef 2..)
# Beta = garch coef 

plot(s1, type="l")

# Simulate garch process
# Remove missing elements
# Specifify model in state space representation 
num = 1000
x <- garch.sim(alpha = c(0.0001, 0.1), beta = 0.45, n = num)

plot(1:num, x, type="l")
plot(cumsum(x), type="l")

# 10% uniform distrbution (note we exclude endpoints)
mis1 <- ceiling(runif(num * 0.7, min = 2, max = (num - 1)))

x_mis1 <- x 
x_mis1[mis1] <- NA

plot(1:num, x_mis1, type="l")
plot(cumsum(x_mis1[!is.na(x_mis1)]), type="l")

# State representation 
A <- array(c(1, 1), dim = c(1, 2, num))
for(k in 1:num){
	if(!is.na(x_mis1[k])) A[,,k] <- diag(1)
}

Phi <- array(c(.14 ,0,0, .45), dim = c(2, 2, num)) # No need for array
Phi <- matrix(c(.1, 0, 0, .45), nrow = 2, ncol = 2)
mu0 <- matrix(c(0, 0), nrow = 2, ncol = 1)
sigma0 <- diag(2)
cQ  = matrix(c(.1, 0, 0, .1), nrow = 2, ncol = 2)

out <- Kfilter1(num = num, y = x, A = A, mu0 = mu0, Sigma0 = sigma0, Phi = Phi, Ups = 0, 
				 				Gam = 0, cQ = cQ, cR = .1, input = 0)

debug(Kfilter1)
undebug(Kfilter1)

plot(cumprod(1 + out$xp), col = "red", type="l")
lines(cumprod(1 + x), col = "blue", type="l")

plot(cumprod(1 + x), col = "blue", type = "l")
dim(out$xp)
plot(exp(out$xp[2,,]) - 1, type="l")
lines(x, type="l")
