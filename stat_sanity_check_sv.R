source("sanity_check_lib.R")

# Optim func
get_optim_func <- function(rtn, A = NULL) {
	out <- function(par) {
		rtn <- rtn
		phi <- par[1]
		sig_eta <- par[2]
		ups <- par[3]
		sigma0 <- sig_eta^2 / (1 - phi^2)	
		mu0 <- ups / (1 - phi)
		
		# We need special case for missing A 	
		if(is.null(A)) {
			a <- array(1, dim = c(1, 1, length(rtn)))
		} else {
			a <- A
		}
		est <- Kfilter1(length(rtn), y = rtn, A = a, Sigma0 = sigma0, mu0 = mu0, Phi = phi, 
										Ups = ups, Gam = -1.27, cQ = sqrt(sig_eta), cR = sqrt((pi^2) / 2), input = 1)
		return(est$like)
	}
}

simulateAverageError_SV <- function(mis_e = .1, n) {
	# garch(1, 1)
	x <- TSA::garch.sim(alpha = c(.1, .3), beta = .45, n = n + 1)
	p <- 100 + cumsum(x)	
	p_rtn <- diff(log(p))
	rtn <- p_rtn - mean(p_rtn) # pr. Harvey et al 1994
	lrtn <- log(rtn^2)	# Base case
	
	# Build two testing series
	mis <- ceiling(runif(n * mis_e, min = 2, max = n - 1))
	lrtn_mis <- lrtn
	lrtn_mis[mis] <- NA
	lrtn_naive <- naive1(lrtn)
	
	init.par = c(.1, 1, 0)
	
	# Fit base case
	base_func <- get_optim_func(lrtn)
	base_qmle <- optim(init.par, base_func, gr = NULL, method = "BFGS", hessian = TRUE)
	
	# State space setup	
	A <- array(1, dim = c(1,1,n))
	base_sigma0 <- base_qmle$par[2]^2 / (1 - base_qmle$par[1]^2)
	base_mu0 <- base_qmle$par[3] / (1 - base_qmle$par[1])
	R <- sqrt((pi^2) / 2)
	
	# Fit base case	
	base_fit <- Kfilter1(num = length(lrtn), y = lrtn, A = a,  Sigma0 = base_sigma0,
												mu0 = base_mu0, Phi = base_qmle$par[1], Ups = 
													base_qmle$par[3], Gam = -1.27, cQ = sqrt(base_qmle$par[2]),
												cR = R, input = 1)

	# Fit naive
	naive_func <- get_optim_func(lrtn_naive)
	naive_qmle <- optim(init.par, naive_func, gr = NULL, method = "BFGS", hessian = TRUE)
	
	naive_sigma0 <- base_qmle$par[2]^2 / (1 - base_qmle$par[1]^2)
	naive_mu0 <- base_qmle$par[3] / (1 - base_qmle$par[1])
	
	naive_fit <- Kfilter1(num = length(lrtn_naive), y = lrtn_naive, A = a, Sigma0 = naive_sigma0,
												mu0 = naive_sigma0, Phi = naive_qmle$par[1], Ups = 
													naive_qmle$par[3], Gam = -1.27, cQ = sqrt(naive_qmle$par[2]),
												cR = R, input = 1)
	
	# Apply EM to lrtn_mis
	# Build observation matrix
	A_mis <- array(NA, dim = c(1, 1, length(lrtn_mis)))
	for(k in 1:(length(lrtn_mis))) {
		if(!is.na(lrtn_mis[k])) A_mis[,,k] <- diag(1)
	}
	
}


# Model specification / 1 test-----------------------------------------------------
raw <- read.csv("DEXUSUK.csv")
# Real data example
p <- raw[, 2]
p <- p[!(p == ".")]
p <- as.numeric(as.vector(p))
# p <- cumsum(x) + 100
# p <- x + 100
rtn <- diff(log(p))
rtn <- rtn - mean(rtn)
rtn <- log((rtn)^2)
library(astsa)
source("stat_sanity_check.R")
# log(y^t) = h_t + eps_t 

fun1 <- function(par) {
	phi <- par[1]
	sig_eta <- abs(par[2])
	ups <- par[3]
	sigma0 <- sig_eta^2 / (1 - phi^2)	
	mu0 <- ups / (1 - phi)
	a <- array(1, dim = c(1, 1, length(sv1[, 2])))
	a[m] <- 0
	est <- Kfilter1(length(sv2), y = sv2, A = a, Sigma0 = sigma0, mu0 = mu0, Phi = phi, 
									Ups = ups, Gam = -1.27, cQ = sig_eta, cR = sqrt((pi^2) / 2), input = 1)
#	cat("Phi: ", phi, " Ups: ", ups, " eta ",  sig_eta, "\n")
	return(est$like)
}
init.par <- c(0.5, 1, 1)

mis <- ceiling(runif(length(sv1[, 2]) * .2, min = 2, max = length(sv1[,2])))
rtn_mis <- rtn
rtn_mis[mis] <- NA

sv2 <- sv1[, 2]
sv2[m] <- NA
a

A <- array(0, dim = c(1, 1, length(rtn)))
for(k in 1:length(rtn)){
	if(!is.na(rtn_mis[k])) A[,,k] <- diag(1)
}

em1 <- EM1(num = length(rtn), y = rtn_mis, A = A, mu0 = mu0, Sigma0 =  sigma0, Phi = out$par[1],
					 cQ = sqrt(out$par[2]), cR = R)

# out <- optim(init.par, fun1, gr = NULL, method = "Nelder-Mead", hessian = FALSE,
# 						 control = list(trace = 1, REPORT = 1))
out <- optim(init.par, fun1, gr = NULL, method = "BFGS", hessian = TRUE,
						 control = list(trace = 1, REPORT = 1))

c("phi" = out$par[1], "sig_eta" = out$par[2]^2, "ups" = out$par[3])


## Estimate with correct Em2 specification
init.par = c(0.9, .1, -.5)

n <- length(rtn)
A <- array(1, c(1,1,n))
sigma0 <- init.par[2]^2 / (1- init.par[1]^2)
mu0 <- init.par[3] / (1 - init.par[1])
ups <- init.par[3]
phi <- init.par[1]
cR <- sqrt((pi^2)/2)
cQ <- sqrt(init.par[2])

out2 <- EM2(num = length(sv2),y = sv2, A = a, Sigma0 = sigma0, mu0 = mu0, Phi = phi, cQ = cQ, cR = cR,
		Ups = ups, Gam = -1.27, input = 1, max.iter = 1000)

debug(EM2)
c("phi" = out2$Phi, "sig_eta" = sqrt(out2$Q), "ups" = out2$Ups)
out2$like
out$value

# constants
A <- array(1, dim = c(1,1,num))
sigma0 <- out$par[2]^2 / (1 - out$par[1]^2)
mu0 <- out$par[3] / (1 - out$par[1])
R <- sqrt((pi^2) / 2)

res <- Kfilter1(length(rtn), y=rtn, A = A, Sigma0 = sigma0, mu0 = mu0, Phi = out$par[1], 
				 Ups = out$par[3], Gam = -1.27, cQ = sqrt(out$par[2]), cR = R, input = 1)

res1 <- Kfilter1(n, y=rtn, A = A, Sigma0 = out2$Sigma0, mu0 = out2$mu0, Phi = out2$Phi, 
				 Ups = out2$Ups, Gam = -1.27, cQ = sqrt(out2$Q), cR = R, input = 1)

sres <- Ksmooth1(length(rtn), y=rtn, A = A, Sigma0 = sigma0, mu0 = mu0, Phi = out$par[1], 
				 Ups = out$par[3], Gam = -1.27, cQ = sqrt(out$par[2]), cR = R, input = 1)

plot(res$xp, type="l", col="red")
lines(sres$xs, type="l", col="blue")
lines(rtn, type="l", col="black")

rtn <- diff(log(p))
plot(abs(rtn), type="l", main="absolute returns and smoothed volatility")
lines(exp(sres$xs / 2), type="l", col = "red")
lines(exp(res$xf / 2), type="l", col = "blue")
# > out$par
# [1]  0.99929225 -0.36026304  0.01363386

# Low p-value i.e we can reject null hypothesis i.e its not stationary
Box.test(rtn, lag = 25, type = "Ljung-Box")
library(tseries)
adf.test(rtn)




