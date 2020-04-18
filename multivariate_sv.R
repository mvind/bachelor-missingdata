# get hist stock prices in tbl2
source("prices.R")

stocks <- tbl2

rtns <- coredata(stocks)
rtns1 <- diff(log(rtns))
rtn_mu <- apply(rtns1, 2, mean)
#rtns1 <- rtns - apply(rtns, 2, mean)
rtns1[, 1] <- rtns1[, 1] - rtn_mu[1]
rtns1[, 2] <- rtns1[, 2] - rtn_mu[2]
rtns1[, 3] <- rtns1[, 3] - rtn_mu[3]

rtns3 <- log(rtns1^2)

# Figure out how to build observtion matrix for multivariate data
m <- 5
rtns1[ceiling(runif(m/100 * nrow(rtns1), min = 2, max = nrow(rtns1) - 1)), 1] <- 0
rtns1[ceiling(runif(m/100 * nrow(rtns1), min = 2, max = nrow(rtns1) - 1)), 2] <- 0
rtns1[ceiling(runif(m/100 * nrow(rtns1), min = 2, max = nrow(rtns1) - 1)), 3] <- 0

build_mv_obs_matrix <- function(tbl) {
	n <- nrow(tbl)
	m <- ncol(tbl)	
	out <- array(diag(1, nrow = m), dim = c(m, m, n))	
	
	for(i in 1:n) {
		mis <- which(tbl[i, ] == 0)
		if(is_empty(mis)) {
			next
		} else {
			diag(out[,, i])[mis] <- 0 
		}
	}
	out
}
A <- build_mv_obs_matrix(rtns1)

init.par <- c(0.9, -.1, -.1)
sigma0 <- init.par[2]^2 / (1- init.par[1]^2)
mu0 <- init.par[3] / (1 - init.par[1])

sigma0 <- diag(sigma0, nrow = 3)
mu0 <- matrix(mu0, ncol = 1, nrow = 3)

psi_sig <- diag(pi^2 / 2, nrow = 3)
eta_sig <- diag(1, nrow = 3)
n <- nrow(rtns1)
A <- array(diag(1, nrow = 3), dim = c(3,3,n))
Phi <- diag(1, nrow = 3)
Gam <- matrix(0, nrow = 3, ncol = 3)
diag(Gam) <- -1.27
Ups <- matrix(0, nrow = 3, ncol = 3)
ut <- matrix(1, ncol = 3, nrow = n)
# Base case no handling
library(astsa)
em_base <- EM2(n, y = rtns1, A = A, Sigma0 = sigma0, mu0 = mu0, Phi = Phi, cQ = eta_sig, cR = chol(psi_sig),
							 Ups = Ups, Gam = Gam, input = ut, max.iter = 100)
undebug(EM2)
debug(EM2)

em_base$Q 
em_base$R
em_base$R / (pi^2/2) 
em_base$Sigma0
em_base$Phi
em_base$mu0
em_base$like


em_base$R %*% t(em_base$R) / (pi^2 / 2)
em_base$Q %*% t(em_base$Q) / 10e-3

Kfilter1(num = n, y = rtns2, A = A, Sigma0 = sigma0, mu0 = mu0, Phi = Phi,
				 cQ = eta_sig, cR = chol(psi_sig), Ups = 0, Gam = Gam, input = ut)
debug(Kfilter1)
undebug(Kfilter1)
