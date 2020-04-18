
sv <- function() {
	# sigma_t = log h_t
	# h_t = ups + phi h_t-1 + sig_eta 
	n = 500
	ups = -.2
	phi = 0.95
	sig_eta = 0.1
	
	h <- vector("numeric", n + 1)
	h[1] <- rnorm(1, ups, sd = sqrt(sig_eta^2 / (1 - phi^2)))
	for(i in 2:(n + 1)) {
		h[i] <- ups + phi * h[i - 1] + rnorm(1, 0, sig_eta)
	}
	y <- h[2:(n+1)] + log(rnorm(n, 0, 1)^2)
	cbind("state" = h[-1], "measurement" = y)
}

get_optim_func <- function(rtn, A = NULL) {
	out <- function(par) {
		phi <- par[1]
		sig_eta <- abs(par[2])
		ups <- par[3]
		sigma0 <- sig_eta^2 / (1 - phi^2)	
		mu0 <- ups / (1 - phi)
		est <- Kfilter1(length(rtn), y = rtn, A = A, Sigma0 = sigma0, mu0 = mu0, Phi = phi, 
											Ups = ups, Gam = -1.27, cQ = sig_eta, cR = sqrt((pi^2) / 2), input = 1)
			#	cat("Phi: ", phi, " Ups: ", ups, " eta ",  sig_eta, "\n")
		return(est$like)
	}
}

main <- function() {
	# Generate stoc proc
	x <- sv()

	# Generate missing data points
	m <- mis(x, p = 25)

	y <- x[, 2] # Observations
	y_mis <- y # Missing observations
	y_mis[m] <- NA
	
	A <- array(1, dim = c(1, 1, length(y)))		

	### MLE 
	# Naive 1
	y_naive1 <- naive1(y_mis)
	n1_optim_func <- get_optim_func(y_naive1, A)
	init.par <- c(0.9, -.1, -.1)
	# 
	# n1_mle <- optim(init.par, n1_optim_func, gr = NULL, method = "BFGS", hessian = T)
	# n1_par <- n1_mle$par
	# 
	# n1_fit <- Kfilter1(num = 500, y = y_naive1, A = A,  mu0 = n1_par[3] / (1 - n1_par[1]), Sigma0 = n1_par[2]^2 / (1 - n1_par[1]^2), 
	# 									 Phi = n1_par[1], Ups = n1_par[3], Gam = -1.27, cQ = abs(n1_par[2]), cR = sqrt((pi^2)/2),
	# 									 input = 1)
	# n1_f <- n1_fit$xf[m]
	
	# Naive 2
	y_naive2 <- naive2(y_mis)
	n2_optim_func <- get_optim_func(y_naive2, A)
	# 
	# n2_mle <- optim(init.par, n2_optim_func, gr = NULL, method = "BFGS", hessian = T)
	# n2_par <- n2_mle$par 
	# 
	# n2_fit <- Kfilter1(num = 500, y = y_naive2, A = A,  mu0 = n2_par[3] / (1 - n2_par[1]), Sigma0 = n2_par[2]^2 / (1 - n2_par[1]^2), 
	# 									 Phi = n2_par[1], Ups = n2_par[3], Gam = -1.27, cQ = abs(n2_par[2]), cR = sqrt((pi^2)/2),
	# 									 input = 1)
	# n2_f <- n2_fit$xf[m]
	# 
	# MLE - no handling 
	A_mis <- A
	A_mis[m] <- 0
	mis_optim_func <- get_optim_func(y_mis, A_mis)
	# mis_fit <- optim(c(0.95, -.1, -.1), mis_optim_func, gr = NULL, method = "BFGS", hessian = T) # Some problem with initpars and NaN likelihood???
	# mis_par <- mis_fit$par
	# 
	# mis_fit <- Kfilter1(num = 500, y = y_mis, A = A_mis,  mu0 = mis_par[3] / (1 - mis_par[1]), Sigma0 = mis_par[2]^2 / (1 - mis_par[1]^2), 
	# 									 Phi = mis_par[1], Ups = mis_par[3], Gam = -1.27, cQ = abs(mis_par[2]), cR = sqrt((pi^2)/2),
	# 									 input = 1)
	# 
	# mle_f <- mis_fit$xf[m]
	# # MLE RESULT VEC
	# cbind("base" = (x[m, 1] - mle_f),
	# 	"naive1" = (x[m, 1] - n1_f),
	# 	"naive2" = (x[m, 1] - n2_f))	
		
	### EM
	init.par = c(0.9, .1, -.1)

	sigma0 <- init.par[2]^2 / (1- init.par[1]^2)
	mu0 <- init.par[3] / (1 - init.par[1])
	ups <- init.par[3]
	phi <- init.par[1]
	cR <- sqrt((pi^2)/2)
	cQ <- sqrt(init.par[2])


	y_em <- y
	y_em[m] <- 0

	# Base case no handling
	em_base <- EM2(num = 500, y = y_em, A = A_mis, Sigma0 = sigma0, mu0 = mu0, Phi = phi, cQ = cQ, cR = cR,
							Ups = ups, Gam = -1.27, input = 1, max.iter = 500)

	base_fit <- Kfilter1(num = 500, y = y_em, A = A_mis, mu0 = em_base$mu0, Sigma0 = em_base$Sigma0,
											 Phi = em_base$Phi, Ups = em_base$Ups, Gam = -1.27, cQ = sqrt(em_base$Q),
											 cR = cR, input = 1)
	base_f <- base_fit$xf[m]

	# naive 1
	em_n1 <- EM2(num = 500, y = y_naive1, A = A, Sigma0 = sigma0, mu0 = mu0, Phi = phi, cQ = cQ, cR = cR,
							 Ups = ups, Gam = -1.27, input = 1, max.iter = 500)

	n1_fit <- Kfilter1(num = 500, y = y_naive1, A = A, mu0 = em_n1$mu0, Sigma0 = em_n1$Sigma0,
											 Phi = em_n1$Phi, Ups = em_n1$Ups, Gam = -1.27, cQ = sqrt(em_n1$Q),
											 cR = cR, input = 1)
	n1_f <- n1_fit$xf[m]

	# naive 2
	em_n2 <- EM2(num = 500, y = y_naive2, A = A, Sigma0 = sigma0, mu0 = mu0, Phi = phi, cQ = cQ, cR = cR,
							 Ups = ups, Gam = -1.27, input = 1, max.iter = 500)

	n2_fit <- Kfilter1(num = 500, y = y_naive2, A = A, mu0 = em_n2$mu0, Sigma0 = em_n2$Sigma0,
										 Phi = em_n2$Phi, Ups = em_n2$Ups, Gam = -1.27, cQ = sqrt(em_n2$Q),
										 cR = cR, input = 1)

	n2_f <- n2_fit$xf[m]

	# # EM RESULTS

	cbind("base" = (x[m, 1] - base_f),
				"naive1" = (x[m, 1] - n1_f),
				"naive2" = (x[m, 1] - n2_f))	
			
}
result <- lapply(1:50, function(i){ cat("Itertion: ", i, "\n") 
	main()})

tbl <- do.call(rbind, result)

summarize(tbl) 

summarize <- function(tbl) {
	cat("Naive - EM", "\n")
	cat("Mean: ", round(apply(tbl, 2, mean), digits = 4), "\n")
	cat("Sd: ", round(apply(tbl, 2, sd), digits = 4))
}	


hist(tbl[,1])
mean(tbl[, 1])
