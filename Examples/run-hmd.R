#########################################################
## Download data from HMD and fit Lee-Carter model
#########################################################
library(ilc)
library(demography)
# download two countries
# Need to register with HMD
username <- ""
password <- ""
US <- extract.years(extract.ages(hmd.mx("USA", username, password, "US"), 0:100), 1960:2010)
seeds <- 1

for(rep in seeds){
	set.seed(rep)
	# remove some elements to be missing
	US2 <- US
	for(i in 1:2){
		damaged <- sample(1:51, 25)
		US2$rate[[i]][sample(1:101, 50), damaged] <- NA
	}
	
	Y0 <- NULL
	Y0[[1]] <- log(t(US$rate[["female"]]))
	Y0[[2]] <- log(t(US$rate[["male"]]))

	Ymiss <- NULL
	Ymiss[[1]] <- log(t(US2$rate[["female"]]))
	Ymiss[[2]] <- log(t(US2$rate[["male"]]))

	mod <- NULL
	mod[[1]] <- lca.rh(US2, mod='lc', interpolate=T, verbose=F, error="gaussian", restype="logrates", series="female")
	mod[[2]] <- lca.rh(US2, mod='lc', interpolate=T, verbose=F, error="gaussian", restype="logrates", series="male")
	
	Ym <- Yr <- NULL
	for(i in 1:2){
		Ym[[i]] <- log(t(mod[[i]]$y$y)) - t(mod[[i]]$residuals$y)
		Yr[[i]] <- as.matrix(t(mod[[i]]$residuals$y))
		for(j in 1:dim(Yr[[i]])[2]){
			Yr[[i]][which(is.na(Ymiss[[i]][,j])), j] <- NA
		} 
	}
	dat <- list(Y0 = Y0, Ymiss=Ymiss, Ym = Ym, Yr = Yr)
}

#########################################################
## Fit model
#########################################################
source("../R/JGL.R")
source("../R/admm.iters.R")
source("../R/gete.R")
source("../R/SSJGL.R")
source("../R/eval.R")
library(Matrix)
library(JGL)
library(corrplot)

post <- "US"
names <- c("USA Female", "USA Male")
for(seed in seeds){
	
	#####################################################
	Y <- NULL
	Yc <- NULL
	for(i in whichgroups){
		Y[[i]] <- dat$Yr[[i]]
		damaged <- which(apply(dat$Yr[[i]], 1, function(x){sum(is.na(x))}) > 0)
		Yc[[i]] <- dat$Yr[[i]][-damaged, ]
	}
	names(Y) <- names
	colnames <- c(0:100) 
	
	#####################################################
	lam1 <- .01
	lam2 <- .01
	v1 <- 1
	lam.eff <- lam1 + c(1:20) * 5
	v0s <- lam1/lam.eff
	fit1 = SSJGL(Y=Y,penalty="fused",lambda0=1, lambda1=lam1,lambda2=lam2, v1 = v1, v0s = v0s, tol.em=1E-4, a=1, b=dim(Y[[1]])[2], doubly=TRUE, normalize=TRUE, c = 0.01)
	nE <- 0
	for(k in 1:length(Y)){
		tmp <- fit1$thetalist[[length(v0s)]][[k]]
		diag(tmp) <- 0
		nE <- nE + sum(tmp!=0)/2
	}
	print(nE)
	########################################################
	if(length(fit1$imputed) > 0){
		truevalues <- rep(NA, length(fit1$imputed))
		meanimpute <- rep(NA, length(fit1$imputed))
		intimpute <- rep(NA, length(fit1$imputed))
		ssimputed <- rep(NA, length(fit1$imputed))
		for(i in 1:dim(fit1$missed)[1]){
			truevalues[i] <- dat$Y0[[fit1$missed[i,1]]][fit1$missed[i,2], fit1$missed[i,3]]
			intimpute[i] <- dat$Ym[[fit1$missed[i,1]]][fit1$missed[i,2], fit1$missed[i,3]]
			meanimpute[i] <- mean(Y[[fit1$missed[i,1]]][fit1$missed[i,2],], na.rm=TRUE)
			ssimputed[i] <- fit1$imputed[i] + dat$Ym[[fit1$missed[i,1]]][fit1$missed[i,2], fit1$missed[i,3]]
		}
	}
	########################################################
	lambda1 <- exp(seq(log(0.001), log(0.01), len = 20)) 
	lambda2<- exp(seq(log(0.001), log(0.01), len = 20)) 
	AIC <- E <- matrix(0, length(lambda1), length(lambda2))
	min <- 1e10
	for(i in 1:length(lambda1)){
		for(j in 1:length(lambda2)){
			fit.cv <- JGL(Y=Yc,penalty="fused",lambda1=lambda1[i],lambda2=lambda2[j], return.whole.theta=TRUE, weights = "sample.size")
			for(k in 1:length(Yc)){
				Sk <- 1/dim(Yc[[k]])[1] * (t(Yc[[k]]) %*% Yc[[k]])
				thetak <- fit.cv$theta[[k]]
				tmp <- thetak
				diag(tmp) <- 0
				Ek <- sum(tmp != 0)/2
				AIC[i, j] <- AIC[i, j] + dim(Yc[[k]])[1] * sum(diag(Sk %*% thetak)) - dim(Yc[[k]])[1] * log(det(thetak)) + 2 * Ek
				E[i, j] <- E[i, j] + Ek
			}
			if(AIC[i, j] < min){
				fit <- fit.cv
				minx <- lambda1[i]
				miny <- lambda2[j]
				min <- AIC[i, j]
			}
			cat(".")
		}
	}
	########################################################
	x <- which.min(abs(E - nE)) %% 20
	y <- (which.min(abs(E - nE)) - x) / 20 + 1
	fit2 <- JGL(Y=Yc,penalty="fused",lambda1=lambda1[x],lambda2=lambda2[y], return.whole.theta=TRUE, weights = "sample.size")
	########################################################
	if(length(fit1$imputed) > 0){
		jglimpute <- rep(NA, length(fit1$imputed))
		aicimpute <- rep(NA, length(fit1$imputed))
		rawimpute <- rep(NA, length(fit1$imputed))
		newY <- newYaic <- newY0 <-  Y
		for(i in 1:length(Y)){
			covAIC <- solve(fit$theta[[i]])
			cov <- solve(fit2$theta[[i]])
			cov0 <- diag(diag(cov(Yc[[i]])))
			cov0 <- nearPD(cov(Yc[[i]]))$mat #+ 0.01 * diag(1, dim(Yc[[i]])[2])
			for(j in 1:dim(Y[[i]])[1]){
				if(sum(is.na(Y[[i]][j, ])) > 0){
					imp <- which(is.na(Y[[i]][j, ]))
					newY[[i]][j, imp] <- as.numeric(t(newY[[i]][j, -imp]) %*% t(cov[imp, -imp] %*% solve(cov[-imp, -imp])))	+ dat$Ym[[i]][j, imp]			
					newYaic[[i]][j, imp] <- as.numeric(t(newYaic[[i]][j, -imp]) %*% t(covAIC[imp, -imp] %*% solve(covAIC[-imp, -imp])))	+ dat$Ym[[i]][j, imp]			
					newY0[[i]][j, imp] <- as.numeric(t(newY[[i]][j, -imp]) %*% t(cov0[imp, -imp] %*% solve(cov0[-imp, -imp]))) + dat$Ym[[i]][j, imp]		
				}
			}
		}
		for(i in 1:dim(fit1$missed)[1]){
			jglimpute[i] <- newY[[fit1$missed[i,1]]][fit1$missed[i,2], fit1$missed[i,3]]
			aicimpute[i] <- newYaic[[fit1$missed[i,1]]][fit1$missed[i,2], fit1$missed[i,3]]
			rawimpute[i] <- newY0[[fit1$missed[i,1]]][fit1$missed[i,2], fit1$missed[i,3]]
		}
	}
	########################################################################
	imputed <- cbind(truevalues, ssimputed, jglimpute, aicimpute, intimpute, rawimpute)
	print(round(c(
			mean(abs(imputed[,1] - imputed[,2])^2),
			mean(abs(imputed[,1] - imputed[,3])^2),
			mean(abs(imputed[,1] - imputed[,4])^2),
			mean(abs(imputed[,1] - imputed[,5])^2),
			mean(abs(imputed[,1] - imputed[,6])^2)
		 ), 6) * 100)
}

#########################################################
# Visualization
#########################################################
methods <- c( "FGL", "DSS-FGL")
names <- c("Female", "Male")
pdf("figures/hmd-cvprec.pdf", width=16, height=16)
par(mfrow = c(2, 2))
for(i in 1:length(Y)){
	corrplot(-cov2cor(fit2$theta[[i]]), method="color", diag=F, is.corr=F)
	title(paste0(methods[1], ": ", names[i]), cex.main = 2)
}
for(i in 1:length(Y)){
	corrplot(-cov2cor(fit1$thetalist[[20]][[i]]), is.corr=F, method='color', diag=F)
	title(paste0(methods[2], ": ", names[i]), cex.main = 2)
}
dev.off()

pdf("figures/hmd-cvcor.pdf", width=16, height=16)
par(mfrow = c(2, 2))
for(i in 1:length(Y)){
	corrplot(cov2cor(solve(fit2$theta[[i]])), method="color", diag=F, is.corr=F)
	title(paste0(methods[1], ": ", names[i]), cex.main = 2)
}
for(i in 1:length(Y)){
	corrplot(cov2cor(solve(fit1$thetalist[[20]][[i]])), is.corr=F, method='color', diag=F)
	title(paste0(methods[2], ": ", names[i]), cex.main = 2)
}
dev.off()
