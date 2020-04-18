 
randomwalk <- function() {
	n = 500 						# length
	
	# x <- vector("numeric", n)	
	# #x[1] <- rnorm(1, mean = 0, sd = 1)
	# x[1] <- 0
	# for(i in 2:n) x[i] <-  0.7*x[i - 1] + rnorm(1, mean = 0, sd = 1)
	x <- arima.sim(n = n +1, list(ar=.8), sd = 1)
	cbind("state" = x[-1], "measurement" = x[-1] + rnorm(n, mean = 0, sd = 1))
}

mis <- function(x, p = 5) {
	# x: random walk 
	# p: % missing data 
	mis <- ceiling(runif(nrow(x) * p / 100, min = 2, max = nrow(x) - 1))
	mis
}

# Optim func
get_optim_func <- function(rtn, A = NULL) {
	out <- function(par) {
		phi <- par[1]
		sig_q <- par[2] # state
		sig_r <- par[3] # obs
		sig0 <- (sig_q)^2 / (1 - phi^2)
		sig0[sig0 < 0] = 0
				
		est <- Kfilter1(length(rtn), y = rtn, A = A, Sigma0 = sig0, mu0 = 0, Phi = phi, 
										Ups = 0, Gam = 0, cQ = sig_q, cR = sig_r, input = 0)
		return(est$like)
	}
}
library(astsa)
main <- function() {
	# Generate random walk 
	x <- randomwalk()
	
	# Generate missing data points 
	m <- mis(x, p = 25)
	
	y <- x[, 2]  # observation
	y_mis <- y # missing observations
	y_mis[m] <- NA
	
	A <- array(1, dim = c(1, 1, length(y)))	

	# Naive 1
	y_naive1 <- naive1(y_mis)
	n1_optim_func <- get_optim_func(y_naive1, A)
	init.par = c(0.8, 1, 1)
	
	# n1_mle <-	optim(init.par, n1_optim_func, gr = NULL, method = "BFGS", hessian = T)
	# n1_par <- n1_mle$par
	# n1_fit <- Kfilter1(num = 500, y = y_naive1, A = A,  mu0 = 0, Sigma0 = n1_par[2]^2 / (1 - n1_par[1]^2), 
	# 											Phi = n1_par[1], Ups = 0, Gam = 0, cQ = n1_par[2], cR = n1_par[3],
	# 											input = 0)
	# n1_f <- n1_fit$xf[m]
	
	# Naive 2
	y_naive2 <- naive2(y_mis)
	n2_optim_func <- get_optim_func(y_naive2, A)
	# n2_mle <- optim(init.par, n2_optim_func, gr = NULL, method = "BFGS", hessian = T)
	# n2_par <- n2_mle$par
	# n2_fit <- Kfilter1(num = 500, y = y_naive2, A = A, mu0 = 0, Sigma0 = n2_par[2]^2 / (1 - n2_par[1]^2),
	# 									 Phi = n2_par[1], Ups = 0, Gam = 0, cQ = n2_par[2], cR = n2_par[3],
	# 									 input = 0)
	# n2_f <- n2_fit$xf[m]
	
	# MLE - no handling
	A_mis <- A
	A_mis[,,m] <- 0
	mis_optim_func <- get_optim_func(y_mis, A_mis)
	#init.par = c(0.8, 1, sqrt(var(log(diff(y_mis^2)), na.rm = T)))
	# mis_fit <- optim(init.par, mis_optim_func, method = "BFGS", gr = NULL, hessian = T)
	# mis_par <- mis_fit$par
	# 
	# # optim fails when missing p is different from 0
	# mle_fit <- Kfilter1(num = 500, y = y_mis, A = A_mis, mu0 = 0, Sigma0 = mis_par[2]^2 / (1 - mis_par[1]^2),
	# 										Phi = mis_par[1], Ups = 0, Gam = 0, cQ = mis_par[2], cR = mis_par[3],
	# 										input = 0)
	# 
	# mle_f <- mle_fit$xf[m]
	# 
	# # MLE RESULT VEC
	# cbind("base" = (x[m, 1] - mle_f),
	# 	"naive1" = (x[m, 1] - n1_f),
	# 	"naive2" = (x[m, 1] - n2_f))	
	### EM
	y_em <- y
	y_em[m] <- 0 # We use 0's instead of NA's

	# Base case no handling
	em_base <- astsa::EM1(num = 500, y = y_em, A = A_mis, mu0 = 0, Sigma0 = 1,
						 Phi = 0.85, cQ = 1.1, cR = 1.2)

	# Note Kfilter usees NA's instead of 0's
	base_fit <- Kfilter1(num = 500, y = y_em, A = A_mis, mu0 = em_base$mu0, Sigma0 = em_base$Sigma0,
											 Phi = em_base$Phi, Ups = 0, Gam = 0, cQ = sqrt(em_base$Q), cR = sqrt(em_base$R), input = 0)
	base_f <- base_fit$xf[m]

	# naive 1
	em_n1 <- astsa::EM1(num = 500, y = y_naive1, A = A, mu0 = 0, Sigma0 = 1,
											Phi = 0.85, cQ = 1.1, cR = 1.2)

	n1_fit <- Kfilter1(num = 500, y = y_naive1, A = A, mu0 = em_n1$mu0, Sigma0 = em_n1$Sigma0,
										 Phi = em_n1$Phi, Ups = 0, Gam = 0, cQ = sqrt(em_n1$Q), cR = sqrt(em_n1$R), input = 0)
	n1_f <- n1_fit$xf[m]

	# naive 2
	em_n2 <- astsa::EM1(num = 500, y = y_naive2, A = A, mu0 = 0, Sigma = 1,
											Phi = 0.85, cQ = 1.1, cR = 1.2)

	n2_fit <- Kfilter1(num = 500, y = y_naive2, A = A, mu0 = em_n2$mu0, Sigma0 = em_n2$Sigma0,
										 Phi = em_n2$Phi, Ups = 0, Gam = 0, cQ = sqrt(em_n2$Q), cR = sqrt(em_n2$R), input = 0)

	n2_f <- n2_fit$xf[m]

	# EM RESULTS
	cbind("base" = (x[m, 1] - base_f),
		"naive1" = (x[m, 1] - n1_f),
		"naive2" = (x[m, 1] - n2_f))
}
main()
result <- lapply(1:50, function(i){ cat("Itertion: ", i, "\n") 
								main()})

tbl <- do.call(rbind, result)

summarize(tbl) 

summarize <- function(tbl) {
	cat("Naive - EM", "\n")
	cat("Mean: ", round(apply(tbl, 2, mean), digits = 4), "\n")
	cat("Sd: ", round(apply(tbl, 2, sd), digits = 6))
}


result[1]
do.call(rbind, result)

#0.26332473  0.175876619  0.30695062
result[2]

125 * 1:50
apply(tbl[1:125, ], 2, mean)
par(mfrow=c(3, 1))
hist(tbl[, 1], breaks = 50, main = "No Imputation - 25% Missing Data", xlab = "Average Error")
hist(tbl[, 2], breaks = 50, main = "Naive1 - 25% Missing Data", xlab = "Average Error")
hist(tbl[,3], breaks = 50, main = "Naive2 - 25% Missing Data", xlab = "Average Error")
