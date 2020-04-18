# This take middle point between the closet data points 
naive1 <- function(vec) {
	#i_mis <- which(is.na(vec) == TRUE)
	for (i in 1:length(vec)) {
		if (is.na(vec[i])) {
			# Find bot and top
			bot <-
				vec[i - 1] # Always assume that previous value has been filled
			
			top <- NA
			j <- 1
			while (is.na(top)) {
				if (!is.na(vec[i + j])) {
					top <- vec[i + j]
				} else {
					j <- j + 1
				}
			}
			# If j > 1 then we have missing data concurrently
			fill <- (top + bot) / 2
			vec[i:(i + j - 1)] <- fill
		}
	}
	vec
}

# This takes the last data point observed 
naive2 <- function(vec) {
	mis <- which(is.na(vec))
	fill <- vector("numeric", length(mis))
	
	for(m in 1:length(mis)) {
		fill[m] <- vec[mis[m] - 1]
		vec[mis[m]] <- fill[m]
	}

	#vec[mis] <- fill
	vec
} 

# Modified functions to handle NA observations
Kfilter1 <- function (num, y, A, mu0, Sigma0, Phi, Ups, Gam, cQ, cR, input) {
	Q = t(cQ) %*% cQ
	R = t(cR) %*% cR
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
	xp = array(NA, dim = c(pdim, 1, num))
	Pp = array(NA, dim = c(pdim, pdim, num))
	xf = array(NA, dim = c(pdim, 1, num))
	Pf = array(NA, dim = c(pdim, pdim, num))
	innov = array(NA, dim = c(qdim, 1, num))
	sig = array(NA, dim = c(qdim, qdim, num))
	x00 = as.matrix(mu0, nrow = pdim, ncol = 1)
	P00 = as.matrix(Sigma0, nrow = pdim, ncol = pdim)
	xp[, , 1] = Phi %*% x00 + Ups %*% ut[1, ]
	Pp[, , 1] = Phi %*% P00 %*% t(Phi) + Q
	B = matrix(A[, , 1], nrow = qdim, ncol = pdim)
	sigtemp = B %*% Pp[, , 1] %*% t(B) + R
	sig[, , 1] = (t(sigtemp) + sigtemp)/2
	siginv = solve(sig[, , 1])
	K = Pp[, , 1] %*% t(B) %*% siginv
	innov[, , 1] = y[1, ] - B %*% xp[, , 1] - Gam %*% ut[1, ]
	xf[, , 1] = xp[, , 1] + K %*% innov[, , 1]
	Pf[, , 1] = Pp[, , 1] - K %*% B %*% Pp[, , 1]
	sigmat = as.matrix(sig[, , 1], nrow = qdim, ncol = qdim)
	like = log(det(sigmat)) + t(innov[, , 1]) %*% siginv %*% 
		innov[, , 1]
	for (i in 2:num) {
		if (num < 2) 
			break
		xp[, , i] = Phi %*% xf[, , i - 1] + Ups %*% ut[i, ]
		Pp[, , i] = Phi %*% Pf[, , i - 1] %*% t(Phi) + Q
		
		if(is.na(y[i, ])) {
			innov[,,i] <- 0
			K <- 0
			sig[,, i] <- 0
		} else {
			B = matrix(A[, , i], nrow = qdim, ncol = pdim)
			siginv = B %*% Pp[, , i] %*% t(B) + R
			sig[, , i] = (t(siginv) + siginv)/2
			siginv = solve(sig[, , i])
			K = Pp[, , i] %*% t(B) %*% siginv
			
			innov[, , i] = y[i, ] - B %*% xp[, , i] - Gam %*% ut[i,]
		}
		
		xf[, , i] = xp[, , i] + K %*% innov[, , i]
		Pf[, , i] = Pp[, , i] - K %*% B %*% Pp[, , i]
		sigmat = matrix(sig[, , i], nrow = qdim, ncol = qdim)
		
		if(is.na(y[i, ])) {
			like = like
		} else {
			like = like + log(det(sigmat)) + t(innov[, , i]) %*% 
				siginv %*% innov[, , i]
			
		}
	}
	like = 0.5 * like
	list(xp = xp, Pp = Pp, xf = xf, Pf = Pf, like = like, innov = innov, 
			 sig = sig, Kn = K)
}
Ksmooth1 <- function (num, y, A, mu0, Sigma0, Phi, Ups, Gam, cQ, cR, input) {
	kf = Kfilter1(num, y, A, mu0, Sigma0, Phi, Ups, Gam, 
								cQ, cR, input)
	pdim = nrow(as.matrix(Phi))
	xs = array(NA, dim = c(pdim, 1, num))
	Ps = array(NA, dim = c(pdim, pdim, num))
	J = array(NA, dim = c(pdim, pdim, num))
	xs[, , num] = kf$xf[, , num]
	Ps[, , num] = kf$Pf[, , num]
	for (k in num:2) {
		J[, , k - 1] = (kf$Pf[, , k - 1] %*% t(Phi)) %*% solve(kf$Pp[,, k])
		xs[, , k - 1] = kf$xf[, , k - 1] + J[, , k - 1] %*% (xs[,, k] - kf$xp[,, k])
		Ps[, , k - 1] = kf$Pf[, , k - 1] + J[, , k - 1] %*% (Ps[,, k] - kf$Pp[,, k]) %*% t(J[,, k - 1])
	}
	x00 = mu0
	P00 = Sigma0
	J0 = as.matrix((P00 %*% t(Phi)) %*% solve(kf$Pp[, , 1]), 
								 nrow = pdim, ncol = pdim)
	x0n = as.matrix(x00 + J0 %*% (xs[, , 1] - kf$xp[, , 1]), 
									nrow = pdim, ncol = 1)
	P0n = P00 + J0 %*% (Ps[, , 1] - kf$Pp[, , 1]) %*% t(J0)
	list(xs = xs, Ps = Ps, x0n = x0n, P0n = P0n, J0 = J0, J = J, 
			 xp = kf$xp, Pp = kf$Pp, xf = kf$xf, Pf = kf$Pf, like = kf$like, 
			 Kn = kf$K)
}
EM1 <- function (num, y, A, mu0, Sigma0, Phi, cQ, cR, max.iter = 100, 
								 tol = 0.001) {
	Phi = as.matrix(Phi)
	pdim = nrow(Phi)
	y = as.matrix(y)
	qdim = ncol(y)
	cvg = 1 + tol
	like = matrix(0, max.iter, 1)
	miss = ifelse(abs(y) > 0, 0, 1)
	#cat("iteration", "   -loglikelihood", "\n")
	for (iter in 1:max.iter) {
		ks = astsa::Ksmooth1(num, y, A, mu0, Sigma0, Phi, Ups = 0, 
									Gam = 0, cQ, cR, input = 0)
		like[iter] = ks$like
		#cat("   ", iter, "        ", ks$like, "\n")
		if (iter > 1) 
			cvg = (like[iter - 1] - like[iter])/abs(like[iter - 
																									 	1])
		if (cvg < 0) 
			break
			#cat("Likelihood not increasing")
			return(	list(Phi = Phi, Q = cQ %*% t(cQ), R = cR %*% t(cR), mu0 = mu0, Sigma0 = Sigma0, 
									 like = like[1:iter], niter = iter, cvg = cvg))
			#stop("Likelihood Not Increasing")
		if (abs(cvg) < tol) 
			break
		
		# Expectation step
		
		Pcs = array(NA, dim = c(pdim, pdim, num))
		eye = diag(1, pdim)
		B = matrix(A[, , num], nrow = qdim, ncol = pdim)
		Pcs[, , num] = (eye - ks$Kn %*% B) %*% Phi %*% ks$Pf[,,num - 1]
		for (k in num:3) {
			Pcs[, , k - 1] = ks$Pf[, , k - 1] %*% t(ks$J[,,k - 2]) + 
				ks$J[,,k - 1] %*% (Pcs[, , k] - Phi %*%	ks$Pf[,,k - 1]) %*% t(ks$J[,,k - 2])
		}
		Pcs[, , 1] = ks$Pf[, , 1] %*% t(ks$J0) + ks$J[, , 1] %*% 
			(Pcs[, , 2] - Phi %*% ks$Pf[, , 1]) %*% t(ks$J0)
		
		S11 = ks$xs[, , 1] %*% t(ks$xs[, , 1]) + ks$Ps[, , 1]
		S10 = ks$xs[, , 1] %*% t(ks$x0n) + Pcs[, , 1]
		S00 = ks$x0n %*% t(ks$x0n) + ks$P0n
		
		B = matrix(A[, , 1], nrow = qdim, ncol = pdim)
		u = y[1, ] - B %*% ks$xs[, , 1] # (y_t - a_t x_t n)
		
		# Observed: oldR = 0
		oldR = diag(miss[1, ], qdim) %*% (t(cR) %*% cR) #?? whats this missign vector doing
		
		R = u %*% t(u) + B %*% ks$Ps[, , 1] %*% t(B) + oldR
		# Calculate
		for (i in 2:num) {
			S11 = S11 + ks$xs[,,i] %*% t(ks$xs[, , i]) + ks$Ps[,,i]
			S10 = S10 + ks$xs[,,i] %*% t(ks$xs[, , i - 1]) + Pcs[,,i]
			S00 = S00 + ks$xs[,,i - 1] %*% t(ks$xs[,,i-1]) + ks$Ps[,,i - 1]
			B = matrix(A[, , i], nrow = qdim, ncol = pdim)
			oldR = diag(miss[i, ], qdim) %*% (t(cR) %*% cR)
			# R = will be old_R   
			R = R + u %*% t(u) + B %*% ks$Ps[, , i] %*% t(B) + 
				oldR
		}
		# Update i.e Maximation step
		Phi = S10 %*% solve(S00)
		Q = (S11 - Phi %*% t(S10))/num
		Q = (Q + t(Q))/2
		cQ = chol(Q)
		R = R/num
		R = diag(diag(R), qdim)
		cR = sqrt(R)
		mu0 = ks$x0n
		Sigma0 = ks$P0n
	}
	list(Phi = Phi, Q = Q, R = R, mu0 = mu0, Sigma0 = Sigma0, 
			 like = like[1:iter], niter = iter, cvg = cvg)
}

