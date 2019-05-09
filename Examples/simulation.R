#################################################################
## This script perform one iteration of simulation analysis
## The setup is similar as in Danaher 2012, 
##     except Omega is generated with G-Wishart distribution
##
#################################################################
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
set.seed(1)
G <- 3

## Specify N, P
p <- 100
N <- rep(150, G)
# smaller block dimension
p0 <- p/10
prep <- p/p0
 
## Generate network
g <- NULL
Graph <- matrix(0, p, p)
for(i in 1:prep){
	sub0 <- matrix(.C("scale_free", G = as.integer(matrix(0, p0, p0)), as.integer(p0), 
            PACKAGE = "BDgraph")$G, p0, p0)
	Graph[((i-1)*p0+1):(i*p0), ((i-1)*p0+1):(i*p0)] <- sub0
}
s1 <- Graph
s2 <- s1
i <- 1
s2[((i-1)*p0+1):(i*p0), ((i-1)*p0+1):(i*p0)] <- 0
s3 <- s2
i <- 2
s3[((i-1)*p0+1):(i*p0), ((i-1)*p0+1):(i*p0)] <- 0
graph.plot <- igraph::graph.adjacency(s1, mode = "undirected", diag = FALSE)
V(graph.plot)$color <- c(rep("red", p0), rep("blue", p0), rep("white", p-p0*2))


# Visualization of the plot
# set.seed(123)
# igraph::plot.igraph(graph.plot, main = "Graph structure", 
# 			layout = layout.fruchterman.reingold,
#             vertex.size = 2, vertex.label = NA, vertex.color=V(graph.plot)$color)

g[[1]] <- s1
g[[2]] <- s2
g[[3]] <- s3

# Simulate from G-Wishart
Theta <- Sigma <- NULL
result = .C("rgwish_c", as.integer(s1), as.double(chol(diag(p))), 
    K = as.double(matrix(0, p, p)), as.integer(3), as.integer(p), as.double(1E-8),
    PACKAGE = "BDgraph")
Theta[[1]] <- matrix(result$K, p, p)
diag <- diag(Theta[[1]])
Theta[[1]] <- Theta[[1]] * s1 + diag(diag)
Sigma[[1]] <- solve(Theta[[1]])
Sigma[[1]] <- cov2cor(Sigma[[1]])
Theta[[1]] <- solve(Sigma[[1]])
Theta[[1]] <- Theta[[1]] * s1 + diag(diag(Theta[[1]]))

# look at summary statistics
if(FALSE){
	absvals <- cov2cor(Theta[[1]])
	diag(absvals) <- 0
	summary(abs(absvals[absvals!=0]))
	sum(Theta[[1]]^2) - sum(diag(Theta[[1]])^2)
}

# Make Sigma correlation matrix
Theta[[2]] <- Theta[[1]] * s2 + diag(diag(Theta[[1]]))
Theta[[3]] <- Theta[[1]] * s3 + diag(diag(Theta[[1]]))
Sigma[[2]] <- solve(Theta[[2]])
Sigma[[3]] <- solve(Theta[[3]])

# Simulate Data
Y <- NULL
Y[[1]] <- rmvnorm(n = N[1], mean = rep(0, p), sigma = Sigma[[1]])
Y[[2]] <- rmvnorm(n = N[2], mean = rep(0, p), sigma = Sigma[[2]])
Y[[3]] <- rmvnorm(n = N[3], mean = rep(0, p), sigma = Sigma[[3]])

# Get plot titles 
mains <- c(paste0("Class 1: n=", N[1], " , edges=", sum(s1)/2), 
		   paste0("Class 2: n=", N[2], " , edges=", sum(s2)/2), 
		   paste0("Class 3: n=", N[3], " , edges=", sum(s3)/2))

penalty <- "fused"
lambda2.fixed <- 0.02

lambda1<- exp(seq(log(0.02), log(1), len = 10)) 
lambda2 <- seq(lambda2.fixed, lambda2.fixed, len=20)
fit0 <- NULL
fit0$fit <- fit0$thetalist <- fit0$problist <- NULL
fit0$nedges <- fit0$SSE <- fit0$fp <- fit0$tp <- fit0$fn <- fit0$tn <- fit0$dKL <- fit0$L1 <- matrix(NA, length(lambda1), G)
fit0$tp.gdiff <- fit0$fp.gdiff <- fit0$fn.gdiff <- rep(NA, length(lambda1))
fit0$tp.diff <- fit0$fp.diff <- fit0$fn.diff <- rep(NA, length(lambda1))
fit0$time <- rep(NA, length(lambda1))
for(i in 1:length(lambda1)){
	start_time <- Sys.time()
	fit <- JGL(Y=Y,penalty=penalty,lambda1=lambda1[i],lambda2=lambda2[i], return.whole.theta=TRUE, weights = "equal")
	print(fit)
	fit0$time[i] <- as.numeric(Sys.time() - start_time, units="secs")
	fit0$fit[[i]] <- fit
	fit0$thetalist[[i]] <-fit$theta
	fit0$problist[[i]] <- matrix(1, p, p)
	for(k in 1:G){
		met <- getmetric(est=fit0$thetalist[[i]][[k]], truth=Theta[[k]], graph=g[[k]])
		fit0$nedges[i,k] <- met$nedges
		fit0$SSE[i,k] <- met$SSE
		fit0$fp[i,k] <- met$fp
		fit0$tp[i,k] <- met$tp
		fit0$fn[i,k] <- met$fn
		fit0$tn[i,k] <- met$tn
		fit0$dKL[i,k] <- met$dKL
		fit0$L1[i,k] <- met$L1
	}
	met2 <- getdiffmetric(est=fit0$thetalist[[i]], truth = Theta, graph=g, tol=1e-2)
	fit0$tp.gdiff[i] <- met2$tp.gdiff
	fit0$fp.gdiff[i] <- met2$fp.gdiff
	fit0$fn.gdiff[i] <- met2$fn.gdiff
	fit0$tp.diff[i] <- met2$tp.diff
	fit0$fp.diff[i] <- met2$fp.diff
	fit0$fn.diff[i] <- met2$fn.diff
	cat(".")
}

## Visualation of partial correlation
# Highlight non-edges
# path.plot(lambda1, fit0, G=g, thres = 0, normalize = T, xlab = expression(lambda[1]), ylab = "Partial correlations",  main = mains, par = c(1, 3), ylim = c(-1, 1), color0="red", color1="gray50", which.top = 0)
# # Highlight edges
# path.plot(lambda1, fit0, G=g, thres = 0, normalize = T, xlab = expression(lambda[1]), ylab = "Partial correlations",  main = mains, par = c(1, 3), ylim = c(-1, 1))

## Fig SSJGL
lam1 <- 1
lam2 <- 1
v1 <- 1
lam.eff <- lam1 + c(1:10) * 5
v0s <- lam1/lam.eff
fit1 = SSJGL(Y=Y,penalty=penalty,lambda0=1, lambda1=lam1,lambda2=lam2, v1 = v1, v0s = v0s, tol.em=1E-4, a=1, b=p, doubly=TRUE, normalize=TRUE)

## Calculate metrics for SSJGL
fit1$nedges <- fit1$SSE <- fit1$fp <- fit1$tp <- fit1$fn <- fit1$tn <- fit1$dKL <- fit1$L1 <- matrix(NA, length(lam.eff), G)
fit1$tp.gdiff <- fit1$fp.gdiff <- fit1$fn.gdiff <- rep(NA, length(lam.eff))
fit1$tp.diff <- fit1$fp.diff <- fit1$fn.diff <- rep(NA, length(lam.eff))
for(i in 1:length(v0s)){
	for(k in 1:G){
		met <- getmetric(est=fit1$thetalist[[i]][[k]], truth=Theta[[k]], graph=g[[k]])
		fit1$nedges[i,k] <- met$nedges
		fit1$SSE[i,k] <- met$SSE
		fit1$fp[i,k] <- met$fp
		fit1$tp[i,k] <- met$tp
		fit1$fn[i,k] <- met$fn
		fit1$tn[i,k] <- met$tn
		fit1$dKL[i,k] <- met$dKL
		fit1$L1[i,k] <- met$L1
	}
	met2 <- getdiffmetric(est=fit1$thetalist[[i]], truth = Theta, graph=g, tol=1e-2)
	fit1$tp.gdiff[i] <- met2$tp.gdiff
	fit1$fp.gdiff[i] <- met2$fp.gdiff
	fit1$fn.gdiff[i] <- met2$fn.gdiff
	fit1$tp.diff[i] <- met2$tp.diff
	fit1$fp.diff[i] <- met2$fp.diff
	fit1$fn.diff[i] <- met2$fn.diff
	cat(".")
}
## Visualization
# path.plot(lam.eff, fit1, G=g, thres = 0, normalize = T, xlab = expression(lambda[1]/v[0]^2), ylab = "Partial correlations",  main = mains, par = c(1, 3), ylim = c(-1, 1), color0="red", color1="gray50", which.top = 0)
# path.plot(lam.eff, fit1, G=g, thres = 0, normalize = T, xlab = expression(lambda[1]/v[0]^2), ylab = "Partial correlations",  main = mains, par = c(1, 3), ylim = c(-1, 1))

# Make the plots for one fitted objects
# The x-axis range is chosen similarly as in Danaher 2012
pdf("figures/simulation1.pdf", width=20, height = 5)
par(mfrow=c(1,4))
# figure 1
xx <- apply(fit0$fp, 1, sum)
yy <- apply(fit0$tp, 1, sum)
xx1 <- sum(fit1$fp[dim(fit1$fp)[1], ])
yy1 <- sum(fit1$tp[dim(fit1$tp)[1], ])
plot(xx[order(xx)], yy[order(xx)], type='l', xlab="FP edges", ylab = "TP edges", xlim = c(0,p*(p-1)*.03), ylim = range(c(yy, yy1)))
points(xx1, yy1, col = "red")

# figure 2
xx <- apply(fit0$nedges, 1, sum)
yy <- apply(fit0$SSE, 1, sum)
xx1 <- sum(fit1$nedges[dim(fit1$nedges)[1], ])
yy1 <- sum(fit1$SSE[dim(fit1$SSE)[1], ])
plot(xx[order(xx)], yy[order(xx)], type='l', xlab="Total edges selected", ylab = "Sum of squared errors", xlim = c(0,p*(p-1)*.03), ylim = range(c(yy, yy1)))
points(xx1, yy1, col = "red")

# figure 3
xx <- fit0$fp.diff 
yy <- fit0$tp.diff 
xx1 <- fit1$fp.diff[length(fit1$fp.diff)]
yy1 <- fit1$tp.diff[length(fit1$tp.diff)] 
plot(xx[order(xx, decreasing = T)], yy[order(xx, decreasing = T)], type='l', xlab="FP differential edges (tol = 0.01)", ylab = "TP differential edges (tol = 0.01)", xlim = c(0, p*(p-1)*.03), ylim = range(c(yy, yy1)))
points(xx1, yy1, col = "red")

# figure 4
xx <- apply(fit0$L1, 1, sum)
yy <- apply(fit0$dKL, 1, sum)
xx1 <- sum(fit1$L1[dim(fit1$L1)[1], ])
yy1 <- sum(fit1$dKL[dim(fit1$dKL)[1], ])
plot(xx[order(xx)], yy[order(xx)], type='l', xlab="L1 norm", ylab = "KL-divergence", xlim = range(c(xx, xx1)), ylim = range(c(yy, yy1)))
points(xx1, yy1, col = "red")
dev.off()





