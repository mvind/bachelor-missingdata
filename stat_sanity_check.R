source("sanity_check_lib.R")
# Simulate random walk with measurement error 

num = 10000
x <- vector("numeric", num)

for(i in 2:num) x[i] <- x[i - 1] + rnorm(1, 0, 0.1)

plot(1:num, x, type="l")

# 10% uniform distrbution (note we exclude endpoints)
mis1 <- ceiling(runif(num * 0.7, min = 2, max = (num - 1)))

x_mis1 <- x 
x_mis1[mis1] <- NA

plot(1:num, x_mis1, type="l")

# TODO implement EM Algorithm on local level model / (Random walk with measurement error)
A <- array(0, dim = c(1, 1, num))
for(k in 1:num){
	if(!is.na(x_mis1[k])) A[,,k] <- diag(1)
}
mu0 <- matrix(0, 1, 1)
sigma0 <- matrix(1, 1, 1)
phi <- matrix(1, 1, 1)

# Høger method 
# Fill using the middle point between the closest two data points

t1 <- naive1(x_mis1)

# Now we filter the filled up timeseries


# out_n1 <- Kfilter0(num, x, A=1, mu0 = matrix(0, 1, 1), Sigma0 = matrix(.1, 1, 1),
# 									 Phi = matrix(1, 1, 1),cQ = matrix(0.5, 1, 1), cR = matrix(var(t1), 1, 1))

# Estimate filtered values 

res <- Kfilter1(num = num, y = x_mis1, Phi = phi, A = A, mu0 = mu0, Sigma0 = sigma0,
				 cQ = diag(1), cR = diag(1), input = 0, Ups = 0, Gam = 0)


em1 <- EM1(num = num, y = x_mis1, A = A, mu0 = mu0, Sigma0 =  sigma0, Phi = phi,
									cQ = diag(1), cR = diag(1))

em_mis1 <- Kfilter1(num = num, y = x_mis1, Phi = em1$Phi, A = A, mu0 = em1$mu0, Sigma0 = em1$Sigma0,
				 cQ = em1$Q, cR = em1$R, input = 0, Ups = 0, Gam = 0)

plot(em_mis1$xf, type="l")
lines(x, type="l", col = "red")
plot(x_mis1, type="l")
lines(em_mis1$xf, col ="red")
mean(x[mis1] - em_mis1$xf[mis1])
mean(x[mis1] - t1[mis1])

mean(x - t1)
mean(x - em_mis1$xf)
mean(x - em_mis1$xp)

out2 <- Ksmooth1(num = num, y = x_mis1, A = A, mu0 = em1$mu0, Sigma0 = em1$Sigma0, Phi = em1$Phi, 
				 Ups = 0, Gam = 0, cQ = em1$Q, cR = em1$R, input = 0)
plot(out2$xs, type="l")

mean(x - out2$xs) 

# Generate Data
set.seed(999); num = 100
x = arima.sim(n=num+1, list(ar=.8), sd=1)
y = ts(x[-1] + rnorm(num,0,1))
# Initial Estimates
u    = ts.intersect(y, lag(y,-1), lag(y,-2))
varu = var(u); 
coru = cor(u)
phi  = coru[1,3]/coru[1,2]
q    = (1-phi^2)*varu[1,2]/phi
r    = varu[1,1] - q/(1-phi^2)

(init.par = c(phi, sqrt(q), sqrt(r)))
Linn  = function(para){
	phi = para[1]; 
	sigw = para[2]; 
	sigv = para[3]
	Sigma0 = (sigw^2)/(1-phi^2); 
	Sigma0[Sigma0<0]=0
	kf = astsa::Kfilter1(num, y, A = array(1, c(1,1,num)), mu0=0, Sigma0 = Sigma0,Phi = phi, cQ = sigw, cR  = sigv,
											 input = 0, Gam = 0, Ups = 0)
	return(kf$like)      }
# Estimation   (partial output shown)
(est = optim(init.par, Linn, gr=NULL, method='BFGS', hessian=TRUE,control=list(trace=1, REPORT=1)))
SE = sqrt(diag(solve(est$hessian)))
cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]),SE)


