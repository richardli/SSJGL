SSJGL <- function(Y,penalty="fused",lambda0,lambda1,lambda2,
				 v1 = 1, 
				 v0s = seq(0.0001, 0.01, len = 10),
				 doubly=FALSE,
				 rho=1, a=1, b =1,
				 maxitr.em=500, tol.em=1e-4,
				 maxitr.jgl=500,tol.jgl=1e-5,
				 warm=NULL, warm.connected=NULL, 
				 truncate=1e-5, 
				 normalize=FALSE, 
				 c=0.1,
				 impute=TRUE){

	# decreasing v0, i.e., increasing penalty, warm start the 0 elements
	if(length(v0s) > 1){
		if(v0s[1] > v0s[2]){
			use.warm.connected <- TRUE
		}else{
			use.warm.connected <- FALSE
		}
	}else{
		use.warm.connected <- FALSE		
	}

	# get dimensions
	p <- dim(Y[[1]])[2]
	K <- length(Y)
	n <- rep(0, K)
	for(k in 1:K) n[k] = dim(Y[[k]])[1]

	# normalize Y?
	meanj <- matrix(NA, K, p)
	if(normalize){
		for(k in 1:K){
		for(j in 1:p){
			meanj[k, j] <- mean(Y[[k]][,j], na.rm = TRUE)
			Y[[k]][,j] = Y[[k]][,j] - meanj[k, j]
		}}
	}
	if(impute){
		klist <- ilist <- jlist <- NULL
		for(k in 1:K){
		for(i in 1:dim(Y[[k]])[1]){
			tmp <- which(is.na(Y[[k]][i, ]))
			if(length(tmp) > 0){
				klist <- c(klist, rep(k, length(tmp)))
				ilist <- c(ilist, rep(i, length(tmp)))
				jlist <- c(jlist, tmp)
			}
		}}
		if(length(klist) == 0) impute <- FALSE
		missed <- cbind(klist, ilist, jlist)
		imputed <- rep(NA, length(klist))
	}else{
		missed <- NULL
		imputed <- NULL
	}
	
	warm.connected <- NULL
	trace_theta <- trace_d <- NULL
	trace_prob_si <- trace_d_si <- NULL
	trace_pi1 <- trace_pi2 <- NULL
	trace_fit <- trace_prob <- NULL
	trace_itr <- trace_diff <- rep(NA, length(v0s))
	theta <- theta_last <- NULL
	d1 <- d2 <- matrix(1, p, p)
	prob1 <- matrix(a/(a+b), p, p)
	prob2 <- matrix(a/(a+b), p, p)
	diag(d1) <- diag(prob1) <- diag(d2) <- diag(prob2) <- 0
	for(k in 1:K){
		theta[[k]] <- theta_last[[k]] <- solve(cov(Y[[k]], use='complete.obs') + diag(c, p))  
	}
	time <- rep(0, length(v0s))
	for(i in 1:length(v0s)){
		start_time <- Sys.time()

		# initialize parameters
		pi_delta <- a/(a+b)
		pi_xi <- a/(a+b)
		pi_delta_last <- NULL
		itr <- 1
		diff <- 1
		v0 <- v0s[i]
		# re-initiate if all 0 in previous case
		tmp <- theta
		klist <- NULL
		for(k in 1:K) diag(tmp[[k]]) <- 0
		for(k in 1:K){
			theta[[k]] <- solve(cov(Y[[k]], use='complete.obs') + diag(c, p)) 
		}
		# estep0 <- gete(p, theta, lambda1, lambda2, v0, v1, pi_delta, pi_xi, penalty, doubly)
		# prob1 <- estep0$prob1

		for(k in 1:K){
			if(sum(unlist(tmp[[k]])) != 0 && sum(diag(theta_last[[k]])) > 1){
				# bring back slab elements over median model?
				theta_last[[k]][prob1 > 0.5] <- 1
				theta[[k]][theta_last[[k]] == 0] <- 0
			}else{
				klist <- c(klist, k)
			}
		}
		theta_last <- NULL
		if(length(klist) > 0) message(paste("Re-initiated to full precision matrices for group", paste(klist, collapse=",")))

		for(itr in 1:maxitr.em){
			if(diff < tol.em) break
			# missing impute step
			if(impute){
				tmp <- getmissing(Y, theta, missed)
				YY <- tmp$YY
				addvar <- tmp$addvar
			}else{
				YY <- Y
				addvar <- NULL
			}
			# E-step
			estep <- gete(p, theta, lambda1, lambda2, v0, v1, pi_delta, pi_xi, penalty, doubly)
			d1 <- estep$d1
			d2 <- estep$d2
			prob1 <- estep$prob1
			prob2 <- estep$prob2
			diag(prob1) <- 0
			if(doubly) diag(prob2) <- 0
			pi_delta <- (a + sum(prob1) - 1) / (a + b + p * (p-1)- 2)
			pi_xi <- (a + sum(prob2) - 1) / (a + b + p * (p-1)- 2)

			# print(summary(as.numeric(prob1)))

			# M-step
			lambda1_current <- lambda1 * d1
			if(doubly){
				lambda2_current <- lambda2 * d2
			}else{
				lambda2_current <- lambda2 * d1
				d2 <- prob2 <- NULL
			}
			mstep <- JGL.adaptive(YY, addvar = addvar, penalty=penalty,lambda0=lambda0, lambda1=lambda1_current,lambda2=lambda2_current,rho=rho, maxiter=maxitr.jgl,tol=tol.jgl,warm=NULL, warm.connected=NULL, return.whole.theta=TRUE, truncate=truncate, normalize=FALSE)
			theta <- mstep$theta

			# compare difference
			if(!is.null(pi_delta_last)){
				diff <- 0
				# diff <- max(abs(pi_delta_last - pi_delta))
				for(k in 1:length(theta_last)) diff <- max(diff, (theta_last[[k]] - theta[[k]])^2)
				if(doubly){
					cat(paste0("Itr ", itr, "  Difference: ", round(diff,6), "  p.slab1: ", round(pi_delta, 4), "  p.slab2: ", round(pi_xi, 10), "\n"))
				}else{
					cat(paste0("Itr ", itr, "  Difference: ", round(diff,6), "  p.slab: ", round(pi_delta, 4), "\n"))
				}
			}
			pi_delta_last <- pi_delta
			theta_last <- theta
		}
		trace_fit[[i]] <- mstep
		trace_prob[[i]] <- prob1
		trace_d[[i]] <- d1
		trace_prob_si[[i]] <- prob2
		trace_d_si[[i]] <- d2
		trace_theta[[i]] <- theta
		trace_pi1[[i]] <- pi_delta
		trace_pi2[[i]] <- pi_xi
		trace_itr[i] <- itr
		trace_diff <- diff
		if(use.warm.connected) warm.connected <- mstep$connected
		time[i] <- as.numeric(Sys.time() - start_time, units="secs")
		cat(paste0("Ladder= ", i, " v0 = ", round(v0,5), " done. Time: ", round(time[i]), "\n"))
	}

	if(impute){
		imputed <- rep(NA, dim(missed)[1])
		for(i in 1:dim(missed)[1]){
			imputed[i] <- YY[[missed[i, 1]]][missed[i, 2], missed[i, 3]]
		}
		imputed <- imputed + meanj[missed[, 3]]
	}

	out <- list(thetalist = trace_theta, pilist = trace_pi1, pilist = trace_pi2, fitlist = trace_fit, itrlist = trace_itr, problist1 = trace_prob, penlist1 = trace_d, problist2 = trace_prob_si, penlist2 = trace_d_si, timelist = time, 
		imputed = imputed, missed = missed)
	return(out)
}
