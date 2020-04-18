load("stocks_raw.Rda")

rtns <- diff(log(tbl)) # Maybe convert to ts obj

y <- log(rtns ^ 2)
y[which(!is.finite(y))] <- -7
num <- nrow(y)

# Initial parameters
phi0 <- 0
phi1 <- 0.95
sQ <- .2
alpha <- apply(y, 2, mean, na.rm = T)
sR0 <- 1
mu1 <- -3
sR1 <- 2

init.par <- c(phi0, phi1, sQ, alpha, sR0, mu1, sR1)

Linn <- function(para) {
	phi0 <- para[1]
	phi1 <- para[2]
	sQ <- para[3]
	alpha <- para[4]
	sR0 <- para[5]
	mu1 <- para[6]
	sR1 <- para[7]
	
	sv <- astsa::SVfilter(num, y, phi0, phi1, sQ, alpha, sR0, mu1, sR1)
	return(sv$like)
}

est <- optim(init.par, Linn, NULL, method = 'BFGS', hessian = TRUE)

phi0 = est$par[1]
phi1 = est$par[2]
sQ = est$par[3]
alpha = est$par[4] 
sR0 = est$par[5]
mu1 = est$par[6]
sR1 = est$par[7]

sv = SVfilter(num, y, phi0, phi1, sQ, alpha, sR0, mu1, sR1)


dim(sv$Pp)

plot(1:length(sv$xp), exp(sv$xp), type="l")
plot(101:length(sv$Pp), sv$Pp[-c(1:100)], type="l")

