gete <- function(p, theta, lambda1, lambda2, v0, v1, pi_delta, pi_xi, penalty, doubly){
	if(doubly) return(gete.doubly(p, theta, lambda1, lambda2, v0, v1, pi_delta, pi_xi, penalty))
	# this function assumes sum_{g1 < g2}
	G <- length(theta)
	pen <- matrix(0, p, p)
	if(penalty == "fused"){
		for(i in 1:(G-1)){
			for(j in (i+1):G){
				pen <- pen + abs(theta[[i]] - theta[[j]])			
			}
		}
		const <- G *(G-1) / 2
	}else if(penalty == "group"){
		for(i in 1:length(theta)){
			pen <- pen + theta[[i]]^2
		}
		pen <- sqrt(pen)
		const <- 1
	}else{
		stop("Penalty needs to be fused or group.")
	}

	abssum <- matrix(0, p, p)
	for(i in 1:G){
		abssum <- abssum + abs(theta[[i]])
	}

	logp0 <- -lambda1/v0 * abssum + log(lambda1/v0) * G - lambda2/v0 * pen + log(lambda2/v0) * const
	logp1 <- -lambda1/v1 * abssum + log(lambda1/v1) * G - lambda2/v1 * pen + log(lambda2/v1) * const

	prob <- exp(logp1) * pi_delta / (exp(logp0) * (1-pi_delta) + exp(logp1) * pi_delta)
	d <- (1-prob)/v0 + prob/v1
	diag(prob) <- 0
	diag(d) <- 0

	d <- (d + t(d)) / 2	
	prob <- (prob + t(prob)) / 2
	out <- list(d1 = d, prob1 = prob, d2=0, prob2 = 0)
	return(out)
}

gete.doubly <- function(p, theta, lambda1, lambda2, v0, v1, pi_delta, pi_xi, penalty){
	# this function assumes sum_{g1 < g2}
	G <- length(theta)
	pen <- matrix(0, p, p)
	if(penalty == "fused"){
		for(i in 1:(G-1)){
			for(j in (i+1):G){
				pen <- pen + abs(theta[[i]] - theta[[j]])			
			}
		}
		const <- G *(G-1) / 2
	}else if(penalty == "group"){
		for(i in 1:length(theta)){
			pen <- pen + theta[[i]]^2
		}
		pen <- sqrt(pen)
		const <- 1
	}else{
		stop("Penalty needs to be fused or group.")
	}

	abssum <- matrix(0, p, p)
	for(i in 1:G){
		abssum <- abssum + abs(theta[[i]])
	}

	prob <- array(NA, dim=c(p,p,3))

	# 00
	prob[,,1] <- -lambda1/v0 * abssum + log(lambda1/v0) * G - lambda2/v0 * pen + log(lambda2/v0) * const + log(1-pi_delta) + log(1-pi_xi) 
	# 10
	prob[,,2] <- -lambda1/v1 * abssum + log(lambda1/v1) * G - lambda2/v0 * pen + log(lambda2/v0) * const + log(pi_delta) + log(1-pi_xi) 
	# 11
	prob[,,3] <- -lambda1/v1 * abssum + log(lambda1/v1) * G - lambda2/v1 * pen + log(lambda2/v1) * const + log(pi_delta) + log(pi_xi) 

	# z <- apply(prob, c(1,2), function(x) sum(exp(x)))
	# for(i in 1:3) prob[,,i] <- exp(prob[,,i]) / z

	# log-sum-exp trick
	mx <- apply(prob, c(1, 2), max)
	prob2 <- prob
	prob2[,,1] <- prob[,,1] - mx
	prob2[,,2] <- prob[,,2] - mx
	prob2[,,3] <- prob[,,3] - mx
	logz <- log(apply(prob2, c(1, 2), function(x){sum(exp(x))})) + mx
	for(i in 1:3) prob[,,i] <- exp(prob[,,i] - logz)


	prob1 <- prob[,,2] + prob[,,3] 
	prob2 <- prob[,,3]
	d1 <- (1-prob1)/v0 + prob1/v1
	d2 <- (1-prob2)/v0 + prob2/v1
	diag(prob1) <- diag(prob2) <- 0
	diag(d1) <- diag(d2) <- 0

	d1 <- (d1 + t(d1)) / 2	
	prob1 <- (prob1 + t(prob1)) / 2
	d2 <- (d2 + t(d2)) / 2	
	prob2 <- (prob2 + t(prob2)) / 2
	out <- list(d1 = d1, prob1 = prob1, d2=d2, prob2 = prob2)

	return(out)
}

# missed has three cols: k, i, j
getmissing <- function(Y, theta, missed){
	YY <- Y
	addvar <- NULL
	for(i in 1:length(Y)){
		addvar[[i]] <- matrix(0, dim(Y[[i]])[2], dim(Y[[i]])[2])
	}
	for(k in unique(missed[,1])){
		tmp <- missed[missed[,1] == k, ]
		cov <- solve(theta[[k]])

		for(i in unique(tmp[, 2])){
			jj <- tmp[tmp[,2] == i, 3]
			oldy <- Y[[k]][i, -jj]
			tmpcov = solve(cov[-jj, -jj])
			newy = as.numeric(t(oldy) %*% t(cov[jj, -jj] %*% tmpcov))
			YY[[k]][i, jj] <- newy
			addvar[[k]][jj, jj] <- addvar[[k]][jj, jj] + cov[jj, jj] - cov[jj, -jj] %*% tmpcov %*% cov[-jj,jj]
		}
	}
	return(list(YY=YY, addvar = addvar))
}