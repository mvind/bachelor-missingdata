library(xts)

## AR(1) process pr. Shumway Stoffer Chp. 6
## NO Missing data (yet)

# Consts
N = 200
psi = 0.5
sig_w = 1
sig_v = 1 

# random AR(1) normal noise 
dates <- as.Date(Sys.Date()) - (0:(N-1))
x <- rep(NA, N)
x[1] <- rnorm(1, 0, 1)
for(i in 2:(N-1)) {
	x[i] <- psi*x[i - 1] + rnorm(1, 0, sig_w)	
}
r1 <- xts(x, dates)
plot(r1)
plot(cumsum(r1[-1]))

hist(r1, breaks = 40)

y <- x + rnorm(length(x), 0, sig_v)

# R1 : signal state
# y  : observation signal (noisy)
ranwalk <- dlm(FF = 1, V = 1, GG = 1, W = 1, m0 = 0, C0 = 1)

r1 <- as.ts(r1)

buildFun <- function(x) {
	dlmModPoly(1, dV=x[1], dW = x[2])
}
fit <- dlmMLE(r1, parm=c(10, 10), build = buildFun)

dlmNile <- buildFun(fit$par)

dlmNile$GG
out <- dlmFilter(r1, dlmNile)

plot(cumsum(out$y[-1]), type="l")

lines(cumsum(out$m), type="o")
qqnorm(out$m)
out$mod$V
