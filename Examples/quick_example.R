###################################################
##  Load dependencies 
###################################################
remove(list = ls())
source("../R/JGL.R")
source("../R/admm.iters.R")
source("../R/gete.R")
source("../R/SSJGL.R")
source("../R/eval.R")
library(BDgraph)
library(tmvtnorm)
library(JGL)
library(corrplot)
library(Matrix)
penalty <- "fused"
# penalty <- "group"

###################################################
##  Simulate data
###################################################
set.seed(1)
P <- 10
PP <- 100
N <- rep(150, 2)
rho <- c(0.7, 0.9)
remove <- 4
Y <- g <- thetatrue <- NULL
s1 <- bdgraph.sim.rho(n = 1, p = P, vis = FALSE, graph = "AR1")$G
diag(s1) <- 0
cat(paste("total edges in Class 1:", sum(s1 /2), "\n"))
s2 <- bdgraph.sim.rho(n = 1, p = P/2, vis = FALSE, graph = "AR1")$G
diag(s2) <- 0
s2 <- bdiag(s2, diag(0, P/2))
cat(paste("total edges Class 2:", sum(s2/2), "\n"))
s0 <- matrix(0, PP - P, PP - P)
g[[1]] <- as.matrix(bdiag(s1, s0))
g[[2]] <- as.matrix(bdiag(s2, s0))
diag(s1) <- diag(s2) <- 1
theta1 <- s1 * bdgraph.sim.rho(n = 1, p = P, vis = FALSE, graph = "AR1", rho=rho[1])$K
theta2 <- s2 * bdgraph.sim.rho(n = 1, p = P, vis = FALSE, graph = "AR1", rho=rho[2])$K
thetatrue[[1]] <- (as.matrix(bdiag(theta1, diag(1, PP-P))))
thetatrue[[2]] <- (as.matrix(bdiag(theta2, diag(1, PP-P))))
mains <- c("Class 1", "Class 2")

Y[[1]] <- rtmvnorm(n = N[1], mean = rep(0, PP), sigma = solve(thetatrue[[1]]))
Y[[2]] <- rtmvnorm(n = N[2], mean  = rep(0, PP), sigma = solve(thetatrue[[2]]))

pdf(paste0("figures/example1-", penalty, "-sim.pdf"), width = 12, height=6.2)
par(mfrow=c(1,2))
for(i in 1:2){
	tmp <- thetatrue[[i]][1:(P*2), 1:(P*2)]
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp, is.corr=F, method="color", diag=T, type="full", cl.lim = range(thetatrue))
	title(paste0("Class", i, ": Simulated precision matrix"), cex.main=2, line=2.5)
}
dev.off()


###################################################
##  Joint graphical lasso
###################################################
lambda.eff <- seq(0.03, 0.3, len = 20)
fit0 <- NULL
fit0$thetalist <- NULL
fit0$problist <- NULL
for(i in 1:length(lambda.eff)){
	fit00 <- JGL(Y=Y,penalty=penalty,lambda1=lambda.eff[i],lambda2=0.1, return.whole.theta=TRUE, weights = "sample.size")
	fit0$thetalist[[i]] <- fit00$theta
	fit0$problist[[i]] <- matrix(1, PP, PP)
	cat(".")
	print(fit00)
}
pdf(paste0("figures/example1-", penalty, "-lasso-prec.pdf"), width = 12, height=6)
path.plot(lambda.eff, fit0, G=g, thres = 0, normalize = F, par = c(1, 2), xlab = expression(lambda[1]), ylab = expression(omega[jk]), main = mains, reverse=T, color1="red1", cex=0.9, cex.main = 1, cex.lab=1.5)
dev.off()

pdf(paste0("figures/example1-", penalty, "-lasso-parcor.pdf"), width = 12, height=6)
path.plot(lambda.eff, fit0, G=g, thres = 0, normalize = T, xlab = expression(lambda[1]), ylab = "Partial correlations",  main = mains, par = c(1, 2), reverse=T, cex=0.9, cex.main = 1, cex.lab=1.5)
dev.off()


###################################################
##  SS Joint graphical lasso
###################################################
lambda1 <- 1 
lambda2 <- 1 
v1 <- 1
lambda.eff2 <- lambda1 + 1 + c(0:20)*10
v0s <- lambda1/lambda.eff2
fit = SSJGL(Y=Y,penalty=penalty,lambda0=1, lambda1=lambda1,lambda2=lambda2, v1 = v1, v0s = v0s, tol.em=0.0001, a=1, b=1, doubly=FALSE, c = 0.01)
pdf(paste0("figures/example1-", penalty, "-ss-prec.pdf"), width = 12, height=6)
path.plot(lambda1/v0s, fit, G=g, thres = 0, normalize = F, xlab = expression(lambda[1]/v[0]), ylab = expression(omega[jk]),  main = mains, par = c(1, 2), color1="red1",cex=0.9, cex.main = 2, cex.lab=1.5)
dev.off()
pdf(paste0("figures/example1-", penalty, "-ss-parcor.pdf"), width = 12, height=6)
path.plot(lambda1/v0s, fit, G=g, thres = 0, normalize = T, xlab = expression(lambda[1]/v[0]), ylab = "Partial correlations",  main = mains, par = c(1, 2), cex=0.9, cex.main = 2, cex.lab=1.5)
dev.off()



###################################################
##  DSS Joint graphical lasso
###################################################
lambda1 <- 1 
lambda2 <- 1 
v1 <- 1
lambda.eff2 <- lambda1 + 1 + c(0:20) * 10
v0s <- (lambda1/lambda.eff2)
fit2 = SSJGL(Y=Y,penalty=penalty,lambda0=1, lambda1=lambda1,lambda2=lambda2, v1 = v1, v0s = v0s, tol.em=0.0001, a=1, b=PP, doubly=TRUE, c = 0.01)

pdf(paste0("figures/example1-", penalty, "-ssd-prec.pdf"), width = 12, height=6)
path.plot(lambda1/v0s, fit2, G=g, thres = 0, normalize = F, xlab = expression(lambda[1]/v[0]), ylab = expression(omega[jk]),  main = mains, par = c(1, 2), color1="red1", cex=0.9, cex.main = 2, cex.lab=1.5)
dev.off()
pdf(paste0("figures/example1-", penalty, "-ssd-parcor.pdf"), width = 12, height=6)
path.plot(lambda1/v0s, fit2, G=g, thres = 0, normalize = T, xlab = expression(lambda[1]/v[0]), ylab = "Partial correlations",  main = mains, par = c(1, 2), cex=0.9, cex.main = 2, cex.lab=1.5)
dev.off()


###################################################
##  Calculate bias
###################################################
SS1 <- SS2 <- rep(0, 2)
for(i in 1:2){
	tmp1 <- fit$thetalist[[20]][[i]]
	tmp2 <- fit2$thetalist[[20]][[i]]
	tmp0 <- thetatrue[[i]]
	SS1[i] <- norm(as.matrix(tmp1) - as.matrix(tmp0), type = "f")#sum((tmp1 - tmp0)^2)
	SS2[i] <- norm(as.matrix(tmp2) - as.matrix(tmp0), type = "f")#sum((tmp2 - tmp0)^2)
}

label <- expression(hat(Omega)[SS-FGL])
if(penalty == "group") label <- expression(hat(Omega)[SS-GGL])
pdf(paste0("figures/example1-", penalty, "-ss-compare.pdf"), width = 12, height=6.3)
par(mfrow=c(1,2))
for(i in 1:2){
	tmp <- (fit$thetalist[[20]][[i]])
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:(P*2), 1:(P*2)], is.corr=F, method="color", diag=T, type="upper")
	tmp <- (thetatrue[[i]])
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:(P*2), 1:(P*2)], is.corr=F, method="color", diag=T, type="lower", add=T)
	text((P*2)-4,(P*2)-4,label, cex=3)
	text(5,5, expression(Omega[true]), cex=3)
	title(paste0(mains[i], " : F-norm=", round(SS1[i],1)), cex.main=2, line=2.5)
}
dev.off()


label <- expression(hat(Omega)[DSS-FGL])
if(penalty == "group") label <- expression(hat(Omega)[DSS-GGL])
pdf(paste0("figures/example1-", penalty, "-ssd-compare.pdf"), width = 12, height=6.3)
par(mfrow=c(1,2))
for(i in 1:2){
	tmp <- (fit2$thetalist[[20]][[i]])
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:(P*2), 1:(P*2)], is.corr=F, method="color", diag=T, type="upper")
	tmp <- (thetatrue[[i]])
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:(P*2), 1:(P*2)], is.corr=F, method="color", diag=T, type="lower", add=T)
	text((P*2)-4,(P*2)-4, label, cex=3)
	text(5,5, expression(Omega[true]), cex=3)
	title(paste0(mains[i], " : F-norm=", round(SS2[i],1)), cex.main=2, line=2.5)
}
dev.off()

###################################################
##  Model selection for JGL
###################################################
min <- 1e10
AIC <- E <- rep(0, length(fit0$thetalist))
for(i in 1:length(fit0$thetalist)){
		for(k in 1:2){
			Sk <- 1/N[k] * (t(Y[[k]]) %*% Y[[k]])
			tmp <- fit0$thetalist[[i]][[k]]
			diag(tmp) <- 0
			Ek <- sum(tmp != 0)/2
			thetak <- fit0$thetalist[[i]][[k]]
			AIC[i] <- AIC[i] + N[k] * sum(diag(Sk %*% thetak)) - N[k] * log(det(thetak)) + 2 * Ek
			E[i] <- E[i] + Ek
		}
}
which <- which.min(AIC)
SS0 <- rep(0, 2)
which2 <- which.min(abs(E - (sum(s1 + s2)- P*2)/2))
for(i in 1:2){
	tmp1 <- fit0$thetalist[[which2]][[i]]
	tmp0 <- thetatrue[[i]]
	# diag(tmp1)  <- diag(tmp0) <- 0
	SS0[i] <- norm(as.matrix(tmp1) - as.matrix(tmp0), type = "f")#sum((tmp1 - tmp0)^2)
}

label <- expression(hat(Omega)[FGL])
if(penalty == "group") label <- expression(hat(Omega)[GGL])
pdf(paste0("figures/example1-", penalty, "-AIC-compare.pdf"), width = 12, height=6.3)
par(mfrow=c(1,2))
for(i in 1:2){
	tmp <- (fit0$thetalist[[which2]][[i]])
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:(P*2), 1:(P*2)], is.corr=F, method="color", diag=T, type="upper")
	tmp <- (thetatrue[[i]])
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:(P*2), 1:(P*2)], is.corr=F, method="color", diag=T, type="lower", add=T)
	text((P*2)-4,(P*2)-4, label, cex=3)
	text(5,5, expression(Omega[true]), cex=3)
	title(paste0(mains[i], " : F-norm=", round(SS0[i],1)), cex.main=2, line=2.5)
}
dev.off()

###################################################
##  Combine plots
###################################################
ylim1 <- c(-4, 0.5)
jpeg(paste0("figures/example1-", penalty, "-combine2.jpeg"), width = 1200*.9, height=900*.9)
label <- expression(hat(Omega)[FGL])
par(mfrow = c(3, 4), mar = c(4,5,4,2))
path.plot(lambda.eff, fit0,  G=g, thres = 0, normalize = F, par = NULL, xlab = expression(lambda[1]), ylab = expression(omega[jk]), main = mains, reverse=T, color1="red1", cex=0.9, cex.main = 2.2, cex.lab=1.5, vline=c(which, which2), position = "none", ylim =ylim1)
par(mar=c(3,2,2,2)+2)
for(i in 1:2){
	tmp <- -cov2cor(fit0$thetalist[[which2]][[i]])
	diag(tmp) <- 0
	tmp[tmp == 0] <- NA
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:P, 1:P], is.corr=T, method="color", diag=T, type="upper", mar=c(3,2,2,2)+1, na.label = "square", na.label.col = "gray80")
	tmp <- -cov2cor(thetatrue[[i]])
	diag(tmp) <- 0
	tmp[tmp == 0] <- NA
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:P, 1:P], is.corr=T, method="color", diag=F, type="lower", add=T, na.label = "square", na.label.col = "gray80")
	text(P-2.5,P-2.5, label, cex=2)
	text(2.5,2.5, expression(Omega[true]), cex=2)
	title(paste0(mains[i], " : F-norm=", round(SS0[i],1)), cex.main=2.4, line=1)
}
label <- expression(hat(Omega)[SS-FGL])
par(mar = c(4,5,4,2))
path.plot(lambda1/v0s, fit, G=g, thres = 0, normalize = F, xlab = expression(lambda[1]/v[0]), ylab = expression(omega[jk]),  main = mains, par = NULL, color1="red1", cex=0.9, cex.main = 2.2, cex.lab=1.5, position = "none", ylim =ylim1)
par(mar=c(3,2,2,2)+2)
for(i in 1:2){
	tmp <- -cov2cor(fit$thetalist[[20]][[i]])
	diag(tmp) <- 0
	tmp[tmp == 0] <- NA
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:P, 1:P], is.corr=T, method="color", diag=T, type="upper", , mar=c(3,2,2,2)+1, na.label = "square", na.label.col = "gray80")
	tmp <- -cov2cor(thetatrue[[i]])
	diag(tmp) <- 0
	tmp[tmp == 0] <- NA
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:P, 1:P], is.corr=T, method="color", diag=F, type="lower", add=T, na.label = "square", na.label.col = "gray80")
	text(P-1.5,P-2.5, label, cex=2)
	text(2.5,2.5, expression(Omega[true]), cex=2)
	title(paste0(mains[i], " : F-norm=", round(SS1[i],1)), cex.main=2.4, line=1)
}
label <- expression(hat(Omega)[DSS-FGL])
par(mar = c(4,5,4,2))
path.plot(lambda1/v0s, fit2, G=g, thres = 0, normalize = F, xlab = expression(lambda[1]/v[0]), ylab = expression(omega[jk]),  main = mains, par = NULL, color1="red1", cex=0.9, cex.main = 2.2, cex.lab=1.5, position = "none", ylim =ylim1)
par(mar=c(3,2,2,2)+2)
for(i in 1:2){
	tmp <- -cov2cor(fit2$thetalist[[20]][[i]])
	diag(tmp) <- 0
	tmp[tmp == 0] <- NA
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:P, 1:P], is.corr=T, method="color", diag=T, type="upper", , mar=c(3,2,2,2)+1, na.label = "square", na.label.col = "gray80")
	tmp <- -cov2cor(thetatrue[[i]])
	diag(tmp) <- 0
	tmp[tmp == 0] <- NA
	colnames(tmp) <- rownames(tmp) <- rep("", dim(tmp)[1])
	corrplot(tmp[1:P, 1:P], is.corr=T, method="color", diag=F, type="lower", add=T, na.label = "square", na.label.col = "gray80")
	text(P-1.5,P-2.5, label, cex=2)
	text(2.5,2.5, expression(Omega[true]), cex=2)
	title(paste0(mains[i], " : F-norm=", round(SS2[i],1)), cex.main=2.4, line=1)
}
dev.off()

###################################################
##  Print out results
###################################################
# AIC
getmetric(est=fit0$thetalist[[which]][[1]], truth=thetatrue[[1]], graph=g[[1]])
getmetric(est=fit0$thetalist[[which]][[2]], truth=thetatrue[[2]], graph=g[[2]])

# true sparsity
getmetric(est=fit0$thetalist[[which2]][[1]], truth=thetatrue[[1]], graph=g[[1]])
getmetric(est=fit0$thetalist[[which2]][[2]], truth=thetatrue[[2]], graph=g[[2]])

# dss-fgl
getmetric(est=fit2$thetalist[[length(v0s)]][[1]], truth=thetatrue[[1]], graph=g[[1]])
getmetric(est=fit2$thetalist[[length(v0s)]][[2]], truth=thetatrue[[2]], graph=g[[2]])