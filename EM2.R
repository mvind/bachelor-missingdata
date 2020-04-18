# OLD 
EM2 <- function(num, y, A, mu0, Sigma0, Phi, cQ, cR, Ups, Gam, input, max.iter = 100, 
								tol = 0.001) {
	Phi = as.matrix(Phi)	
	pdim = nrow(Phi)	
	y = as.matrix(y)
	qdim = ncol(y)
	rdim = ncol(as.matrix(input))
	
	if (max(abs(Ups)) == 0) 
		Ups = matrix(0, pdim, rdim)
	if (max(abs(Gam)) == 0) 
		Gam = matrix(0, qdim, rdim)
	Ups = as.matrix(Ups)
	Gam = as.matrix(Gam)
	ut = matrix(input, num, rdim)
	
	like = matrix(0, nrow = max.iter, ncol = 1)
	
	miss = ifelse(abs(y) > 0, 0, 1)
	cvg = 1 + tol
	
	# Start EM Iterations 
	for(iter in 1:max.iter)	{
		ks = Ksmooth1(num, y, A, mu0, Sigma0, Phi, Ups = Ups, Gam = Gam, 
									cQ = cQ, cR = cR, input = input)	
		like[iter] = ks$like	
		
		if(iter > 1) {
			cvg = (like[iter - 1] - like[iter])/abs(like[iter - 1])
		}
		print(cvg)
		if(cvg < 0)	 stop("Likelihood not increasing")

		# one lag smoothed covariance 
		Pcs = array(NA, dim = c(pdim, pdim, num))
		eye=diag(1,pdim)
		B = matrix(A[,,num], nrow=qdim, ncol=pdim)
		Pcs[,,num]=(eye-ks$Kn%*%B)%*%Phi%*%ks$Pf[,,num-1]
		for (k in num:3) {
			Pcs[, , k - 1] = ks$Pf[, , k - 1] %*% t(ks$J[,,k - 2]) + 
				ks$J[,,k - 1] %*% (Pcs[, , k] - Phi %*%	ks$Pf[,,k - 1]) %*% t(ks$J[,,k - 2])
		}
		Pcs[, , 1] = ks$Pf[, , 1] %*% t(ks$J0) + ks$J[, , 1] %*% 
			(Pcs[, , 2] - Phi %*% ks$Pf[, , 1]) %*% t(ks$J0)
		
		Sxx = ks$xs[,, 1] %*% t(ks$xs[,, 1]) + ks$Ps[,, 1] # S00
		Sxb = ks$xs[,, 1] %*% t(ks$x0n)	+ Pcs[,, 1] # S10
		Sbb = ks$x0n %*% t(ks$x0n) + ks$P0n # S00
		 
		# We assume ut_0 = ut_1
		Syy = y[1, ] %*% t(y[1, ])
		Syx = y[1, ] %*% t(ks$xs[,, 1])
		Syu = y[1, ] %*% t(ut[1 ,])
		
		Sxu1 = ks$xs[,,1] %*% t(ut[1, ])
		Sbu = ks$x0n %*% t(ut[1, ])
		Suu1 = ut[1,] %*% t(ut[1, ])
		Sxu2 = Sxu1
		Suu2 = Suu1
		
		# Missing observation? Use previous R covariance 
		oldR = diag(miss[1, ], qdim)	%*% (t(cR) %*% cR)
		B = matrix(A[,, 1], nrow = qdim, ncol = pdim)
		u = y[1, ] - B %*% ks$xs[,, 1] - Gam %*% ut[1, ]
		R = u %*% t(u) + B %*% ks$Ps[,, 1] %*% t(B) + oldR
		# 	
		# Update parameters
		for(i in 2:num) {
			Sxx = Sxx + ks$xs[,,i] %*% t(ks$xs[,,i]) + ks$Ps[,,i]
			Sxb = Sxb + ks$xs[,, i] %*% t(ks$xs[,, i - 1]) + Pcs[,,i]
			Sbb = Sbb + ks$xs[,,i - 1] %*% t(ks$xs[,,i - 1]) + ks$Ps[,,i - 1]
			
			Syy = Syy + y[i, ] %*% t(y[i, ])
			Syx = Syx + y[i, ] %*% t(ks$xs[,, i])
			Syu = Syu + y[i, ] %*% t(ut[i, ])
			
			Sxu1 = Sxu1 + ks$xs[,,i] %*% t(ut[i - 1, ])
			Sbu = Sbu + ks$xs[,, i - 1] %*% t(ut[i - 1, ])
			Suu1 = Suu1 + ut[i - 1, ] %*% t(ut[i - 1, ])
			Sxu2 = Sxu2 + ks$xs[,, i] %*% t(ut[i, ])
			Suu2 = Suu2 + ut[i, ] %*% t(ut[i, ])
			
			B = matrix(A[,,i], nrow = qdim, ncol = pdim)
			u = y[i, ] - B %*% ks$xs[,, i] - Gam %*% ut[i, ]
			oldR = diag(miss[i, ], qdim) %*% (t(cR) %*% cR)
			R = R + u %*% t(u) + B %*% ks$Ps[,, i] %*% t(B) + oldR
		}
		# Update parameters // Maximation Step
		# if(input == 0) { # Correct matrix size for no input
		# 	Phi = Sxb %*% solve(Sbb)
		# 	Ups = matrix(0, pdim, rdim)
		# } else {
		# 	ab_inv = rbind(cbind(Sbb, Sbu), cbind(t(Sbu), Suu1))
		# 	ab_inv2 = cbind(Sxb, Sxu1)	
		# 	ab = ab_inv2 %*% solve(ab_inv) # How to extract Phi and Ups? 
		# 	Phi = ab[1:pdim, 1:pdim]
		# 	Ups = ab[1:pdim, (pdim + 1):ncol(ab)]
		# }
		
		Q = Sxx - Sxb %*% t(Phi) - Phi %*% t(Sxb) - Sxu1 %*% t(Ups) - Ups %*% t(Sxu1) +
			Phi %*% Sbu %*% t(Ups) + Ups %*% t(Sbu) %*% t(Phi) + Phi %*% Sbb %*% t(Phi) +
			Ups %*% Suu1 %*% t(Ups)
		Q = Q / num	
		cQ = chol(Q)	

		# if(input == 0) {
		# 	A = Syx %*% solve(Sxx)
		# 	A = array(A, dim = c(qdim, pdim, num))
		# 	Gam = matrix(0, qdim, rdim)
		# } else {
		# cd_inv = rbind(cbind(Sxx, Sxu2), cbind(t(Sxu2), Suu2))
		# 	cd_inv2 = cbind(Syx, Syu)
		# 	cd = cd_inv2 %*% solve(cd_inv)
		# 	#A = cd[1:qdim, 1:pdim]
		# 	#A = array(A, dim = c(qdim, pdim, num))
		# Gam = cd[1:qdim, (pdim + 1):ncol(cd)]
		# }
		# 
		R = Syy - Syx %*% t(A[,,num]) - A[,, num] %*% t(Syx) + Syu %*% t(Gam) - Gam %*% t(Syu) +
			A[,,num] %*% Sxu2 %*% t(Gam) + Gam %*% t(Sxu2) %*% t(A[,,num]) + A[,,num] %*% Sxx %*% t(A[,,num]) +
			Gam %*% Suu2 %*% t(Gam)

		R = R/num
		#R = diag(diag(R), qdim)
		diag(R) <- rep((pi^2) / 2, 3)
		cR = sqrt(abs(R))
		#cat("Phi: ", Phi, " Ups: ", Ups, "\n")
		mu0 = ks$x0n
		Sigma0 = ks$P0n
		# mu0 = mu0
		# Sigma0 = Sigma0
		#mu0 <- Ups/(1-Phi)
		#Sigma0 <- (cQ^2) / (1 - Phi^2)
		
	}
	# Return optimal parameters	
	return(list(A = A[,, 1], Gam = Gam, Ups = Ups, Phi = Phi, Q = cQ, R = cR, mu0 = mu0, Sigma0 = Sigma0,
			 like = like[1:iter], niter = iter, cvg = cvg))
}


# # Block matrix testing 
# a <- c(1, 1) %*% t(c(1, 1))
# b <- c(1, 1) %*% t(c(2, 2, 2))
# 
# c <- rbind(cbind(a, b), cbind(b, a))
# c[1:2, 1:2]
# 
# # EM0 Example from stoffer
# set.seed(999); num = 100
# x = arima.sim(n=num+1, list(ar = .8), sd=1) 
# y = ts(x[-1] + rnorm(num,0,1))
# # Initial Estimates (same as Example 6.6)
# u = ts.intersect(y, lag(y,-1), lag(y,-2))
# varu = var(u); coru = cor(u)
# phi = coru[1,3]/coru[1,2]
# q    = (1-phi^2)*varu[1,2]/phi
# r    = varu[1,1] - q/(1-phi^2)
# 
# em = astsa::EM1(num, y, A=array(1, c(1,1,num)), mu0=0, Sigma0=2.8, Phi=phi, cQ=sqrt(q), cR=sqrt(r),
# 				 max.iter=75, tol=.00001)
# debug(astsa::EM1)
# 
# lem2 = EM2(num, y, A=array(1, dim=c(1,1,num)), mu0=0, Sigma0 = 2.8, Phi = phi, cQ = sqrt(q), cR = sqrt(r),
# 					max.iter = 75, tol = .00001, input = 0, Gam = 3, Ups = 2)
# debug(EM2)
# 
# em$Phi - lem2$Phi
# em$mu0 - lem2$mu0
# em$Sigma0 - lem2$Sigma0
# em$R - lem2$R
# em$Q - lem2$Q
# 
# out <- Kfilter1(num, y, A = array(lem2$A, c(1,1,num)), mu0 = lem2$mu0, Sigma0 = lem2$Sigma0,
# 				Phi = lem2$Phi, Ups = lem2$Ups, Gam = lem2$Gam, cQ = sqrt(lem2$Q), cR = sqrt(lem2$R),
# 				input = 1)
# 
# 
# Kfilter1(num, y, A = array(1, dim = c(1, 1, num)), mu0 = 0, Sigma0 = 2.8, 
# 				 Phi = phi, Ups = 0, Gam = 0, cQ = sqrt(q), cR = sqrt(r), input = 0)
# undebug(Kfilter1)
# debug(Kfilter1)
