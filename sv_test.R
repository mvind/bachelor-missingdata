library(astsa)

y = log(nyse ^ 2)
num = length(y)
# Initial Parameters
phi0 = 0
phi1 = .95
sQ = .2
alpha = mean(y)
sR0 = 1
mu1 = -3
sR1 = 2
init.par = c(phi0, phi1, sQ, alpha, sR0, mu1, sR1) # Innovations Likelihood
Linn = function(para) {
	phi0 = para[1]
	phi1 = para[2]
	sQ  = para[3]
	alpha = para[4]
	sR0  = para[5]
	mu1  = para[6]
	sR1 = para[7]
	sv = SVfilter(num, y, phi0, phi1, sQ, alpha, sR0, mu1, sR1)
	return(sv$like)
}
# Estimation
(est = optim(
	init.par,
	Linn,
	NULL,
	method = 'BFGS',
	hessian = TRUE,
	control = list(trace = 1, REPORT = 1)
))
SE = sqrt(diag(solve(est$hessian)))
u = cbind(estimates = est$par, SE)
rownames(u) = c('phi0', 'phi1', 'sQ', 'alpha', 'sigv0', 'mu1', 'sigv1')
u
# Graphics (need filters at the estimated parameters)
phi0 = est$par[1]
phi1 = est$par[2]
sQ = est$par[3]
alpha = est$par[4] sR0 = est$par[5]
mu1 = est$par[6]
sR1 = est$par[7]
sv = SVfilter(num, y, phi0, phi1, sQ, alpha, sR0, mu1, sR1)
# densities plot (f is chi-sq, fm is fitted mixture)
x = seq(-15, 6, by = .01)
f = exp(-.5 * (exp(x) - x)) / (sqrt(2 * pi))
f0 = exp(-.5 * (x ^ 2) / sR0 ^ 2) / (sR0 * sqrt(2 * pi))
f1 = exp(-.5 * (x - mu1) ^ 2 / sR1 ^ 2) / (sR1 * sqrt(2 * pi))
fm = (f0 + f1) / 2
plot(x, f, type = 'l')
lines(x, fm, lty = 2, lwd = 2)
dev.new()
Time = 701:1100
plot (
	Time,
	nyse[Time],
	type = 'l',
	col = 4,
	lwd = 2,
	ylab = '',
	xlab = '',
	ylim = c(-.18, .12)
)
lines(Time, sv$xp[Time] / 20, lwd = 2, col = 6)

par(mfrow=c(2, 1))
plot(exp(sv$xp), type = "l")
plot(nyse)
