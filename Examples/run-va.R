source("../R/JGL.R")
source("../R/admm.iters.R")
source("../R/gete.R")
source("../R/SSJGL.R")
source("../R/eval.R")
library(openVA)
library(corrplot)
library(JGL)
set.seed(123)
PHMRC_all <- read.csv(getPHMRC_url("adult"))
causes <-  as.character(PHMRC_all[(match(1:34, PHMRC_all$va34)), "gs_text34"])
symps <- c(paste0("a2_0", c(1,3,8)), paste0("a2_", c(15,22,24,26,28,33,37,41,48,54,58,62,65,68,70,73,76,79,83,86)), c("a3_08", "a4_03", "a4_04", "g1_07a"))
data <- PHMRC_all[, c("gs_text34", "site", symps)]
sympstext <- c("ill", "fever", "rash", "ulcer", "yellow discoloration", "ankle swelling", "puffiness face", "puffiness body", "cough", "difficulty breathing", "fast breathing", "liquid stools", "vomit", "difficulty swallowing", "belly pain", "protruding belly", "mass belly", "headaches", "stiff neck", "unconsciousness", "confusion", "convulsion", "paralysis", "period overdue", "tobacco", "cigarettes", "age")
cbind(symps, sympstext)
head(data)
for(i in 3:(length(symps)+2)) data[,i] <- as.numeric(data[,i])
colnames(data)[3:(length(symps)+2)] <- sympstext

par(mfrow = c(2,3))
causes2 <- c("Stroke", "Pneumonia", "AIDS")
group <- 1
Y <- NULL
counter <- 1
for(i in 1:length(causes2)){
	tmp<- log(as.matrix(data[data[,group] == causes2[i], -c(1,2)]) + 1)
	print(paste0("Class ", i, " n = ", dim(tmp)[1]))
	if(dim(tmp)[1] > 0){
		Y[[counter]] <- tmp
		counter <- counter + 1
	}
}

penalty <- "fused"

# SSJGL
g <- NULL
for(i in 1:length(Y)) g[[i]] <- matrix(1, dim(Y[[1]])[2], dim(Y[[1]])[2])
lam1 <- 1
lam2 <- 0.5
v1 <- 1
lam.eff <- lam1 + c(0:20)
v0s <- (lam1/lam.eff)
fit1 = SSJGL(Y=Y,penalty="fused",lambda0=1, lambda1=lam1,lambda2=lam2, v1 = v1, v0s = v0s, tol.em=1E-4, a=1, b=dim(Y[[1]])[2], doubly=TRUE, normalize=TRUE)
fit1g = SSJGL(Y=Y,penalty="group",lambda0=1, lambda1=lam1,lambda2=lam2*2, v1 = v1, v0s = v0s, tol.em=1E-4, a=1, b=dim(Y[[1]])[2], doubly=TRUE, normalize=TRUE)

###########################################################
# DSS-FGL
###########################################################
pdf("figures/phmrc-DSSFGL.pdf", width=18, height=6)
set.seed(1)
par(mfrow = c(1,3))
common <- matrix(0, dim(Y[[1]])[2], dim(Y[[1]])[2])
for(i in 1:length(Y)){
	sign <- fit1$thetalist[[length(v0s)]][[i]]
	sign[sign != 0] <- 1
	common <- common + sign
}
common[common != length(Y)] <- 0
common[common == length(Y)] <- 1
common.g <- graph.adjacency(common, mode = "upper", weighted = TRUE)
for(i in 1:length(Y)){
	adj <- -cov2cor(fit1$thetalist[[length(v0s)]][[i]]) 
	diag(adj) = 0
	gadj = graph.adjacency(adj, mode = "upper", weighted = TRUE)
	E(gadj)$weight = get.edge.attribute(gadj, "weight") * 15
	E(gadj)$color = 6
	ee <- get.edgelist(gadj)
	e0 <- get.edgelist(common.g)
	is.common <- which(apply(ee, 1, paste0, collapse="-") %in% apply(e0, 1, paste0, collapse="-"))
	E(gadj)$color[is.common] = 2
	# by hand!
	order <- c(7:11, 1:6, 12:27)
	if(i == 1)	layout <- layout_in_circle(gadj, order = V(gadj)[order])
	vlocs <- rep(1, length(V(gadj)))
	vlocs[20] = 0.8
	vlocs[21] = 1.1
	vlocs[22] = 0.9
	 radian.rescale <- function(x, start=0, direction=1) {
	   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
	   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
	 }
	locs <- radian.rescale(x=1:length(V(gadj)), direction=-1, start=0)[order]

	plot(gadj, vertex.frame.color = "white", layout = layout, 
	    vertex.size = 4, 
	    vertex.color= "darkred",
	    vertex.label.cex = 1.6, 
	    vertex.label.color= "black",
	    vertex.label.dist=vlocs,
	    vertex.label.degree=locs,
	    edge.width=abs(E(gadj)$weight), 
	    edge.color = E(gadj)$color, 
	    edge.lty = E(gadj)$lty)
	title(causes2[i],cex.main=2)

}
dev.off()
###########################################################
# DSS-GGL
###########################################################
pdf("figures/phmrc-DSSGGL.pdf", width=18, height=6)
set.seed(1)
par(mfrow = c(1,3))
common <- matrix(0, dim(Y[[1]])[2], dim(Y[[1]])[2])
for(i in 1:length(Y)){
	sign <- fit1g$thetalist[[length(v0s)]][[i]]
	sign[sign != 0] <- 1
	common <- common + sign
}
common[common != length(Y)] <- 0
common[common == length(Y)] <- 1
common.g <- graph.adjacency(common, mode = "upper", weighted = TRUE)
for(i in 1:length(Y)){
	adj <- -cov2cor(fit1g$thetalist[[length(v0s)]][[i]]) 
	gadj = graph.adjacency(adj, mode = "upper", weighted = TRUE)
	E(gadj)$weight = get.edge.attribute(gadj, "weight") * 15
	E(gadj)$color = 6
	ee <- get.edgelist(gadj)
	e0 <- get.edgelist(common.g)
	is.common <- which(apply(ee, 1, paste0, collapse="-") %in% apply(e0, 1, paste0, collapse="-"))
	E(gadj)$color[is.common] = 2
	# by hand!
	order <- c(7:11, 1:6, 12:27)
	if(i == 1)	layout <- layout_in_circle(gadj, order = V(gadj)[order])
	vlocs <- rep(1, length(V(gadj)))
	vlocs[20] = 0.8
	vlocs[21] = 1.1
	vlocs[22] = 0.9
	 radian.rescale <- function(x, start=0, direction=1) {
	   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
	   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
	 }
	locs <- radian.rescale(x=1:length(V(gadj)), direction=-1, start=0)[order]

	plot(gadj, vertex.frame.color = "white", layout = layout, 
	    vertex.size = 4, 
	    vertex.color= "darkred",
	    vertex.label.cex = 1.6, 
	    vertex.label.color= "black",
	    vertex.label.dist=vlocs,
	    vertex.label.degree=locs,
	    edge.width=abs(E(gadj)$weight), 
	    edge.color = E(gadj)$color, 
	    edge.lty = E(gadj)$lty)
	title(causes2[i],cex.main=2)

}
dev.off()

nE <- 0
for(k in 1:length(Y)){
	tmp <- fit1$thetalist[[length(v0s)]][[k]]
	diag(tmp) <- 0
	nE <- nE + sum(tmp!=0)/2
}
print(nE)

nEg <- 0
for(k in 1:length(Y)){
	tmp <- fit1g$thetalist[[length(v0s)]][[k]]
	diag(tmp) <- 0
	nEg <- nEg + sum(tmp!=0)/2
}
print(nEg)


for(penalty in c("fused", "group")){
	post <- "FGL"
	if(penalty == "group") post <- "GGL"

	lambda1 <- exp(seq(log(0.01), log(0.2), len = 20)) 
	lambda2<- exp(seq(log(0.01), log(0.2), len = 20)) 
	AIC <- E <- matrix(0, length(lambda1), length(lambda2))
	min <- 1e10
	for(i in 1:length(lambda1)){
		for(j in 1:length(lambda2)){
			fit.cv <- JGL(Y=Y,penalty=penalty,lambda1=lambda1[i],lambda2=lambda2[j], return.whole.theta=TRUE, weights = "sample.size")
			for(k in 1:length(Y)){
				Sk <- 1/dim(Y[[k]])[1] * (t(Y[[k]]) %*% Y[[k]])
				thetak <- fit.cv$theta[[k]]
				tmp <- thetak
				diag(tmp) <- 0
				Ek <- sum(tmp != 0)/2
				AIC[i, j] <- AIC[i, j] + dim(Y[[k]])[1] * sum(diag(Sk %*% thetak)) - dim(Y[[k]])[1] * log(det(thetak)) + 2 * Ek
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
	print(minx)
	print(miny)
	####################################################### # FGL/GGL with AIC
	########################################################
	pdf(paste0("figures/phmrc-", post, "-AIC.pdf"), width=18, height=6)
	set.seed(1)
	par(mfrow = c(1,3))
	common <- matrix(0, dim(Y[[1]])[2], dim(Y[[1]])[2])
	for(i in 1:length(Y)){
		sign <- fit$theta[[i]]
		sign[sign != 0] <- 1
		common <- common + sign
	}
	common[common != length(Y)] <- 0
	common[common == length(Y)] <- 1
	common.g <- graph.adjacency(common, mode = "upper", weighted = TRUE)
	for(i in 1:length(Y)){
		adj <- -cov2cor(fit$theta[[i]]) #adjs[[i]]
		# adj <- adj * (1-common)
		diag(adj) = 0
		gadj = graph.adjacency(adj, mode = "upper", weighted = TRUE)
		E(gadj)$weight = get.edge.attribute(gadj, "weight") * 15
		# E(gadj)$lty = sign(get.edge.attribute(gadj, "weight")) 
		E(gadj)$color = 6
		ee <- get.edgelist(gadj)
		e0 <- get.edgelist(common.g)
		is.common <- which(apply(ee, 1, paste0, collapse="-") %in% apply(e0, 1, paste0, collapse="-"))
		E(gadj)$color[is.common] = 2
		# by hand!
		order <- c(7:11, 1:6, 12:27)
		if(i == 1)	layout <- layout_in_circle(gadj, order = V(gadj)[order])
		vlocs <- rep(1, length(V(gadj)))
		vlocs[20] = 0.8
		vlocs[21] = 1.1
		vlocs[22] = 0.9
		 radian.rescale <- function(x, start=0, direction=1) {
		   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
		   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
		 }
		locs <- radian.rescale(x=1:length(V(gadj)), direction=-1, start=0)[order]

		plot(gadj, vertex.frame.color = "white", layout = layout, 
		    vertex.size = 4, 
		    vertex.color= "darkred",
		    vertex.label.cex = 1.6, 
		    vertex.label.color= "black",
		    vertex.label.dist=vlocs,
		    vertex.label.degree=locs,
		    edge.width=abs(E(gadj)$weight), 
		    edge.color = E(gadj)$color, 
		    edge.lty = E(gadj)$lty)
		title(causes2[i],cex.main=2)

	}
	dev.off()

	# Same number of edges
	if(penalty == "fused"){
		EE <- nE
	}else{
		EE <- nEg
	}
	x <- which.min(abs(E - EE)) %% 20
	y <- (which.min(abs(E - EE)) - x) / 20 + 1
	fit2 <- JGL(Y=Y,penalty=penalty,lambda1=lambda1[x],lambda2=lambda2[y], return.whole.theta=TRUE, weights = "sample.size")
	
	#######################################################
	# FGL/GGL same number of edges
	#######################################################
	pdf(paste0("figures/phmrc-", post, "-nE.pdf"), width=18, height=6)
	set.seed(1)
	par(mfrow = c(1,3))
	common <- matrix(0, dim(Y[[1]])[2], dim(Y[[1]])[2])
	for(i in 1:length(Y)){
		sign <- fit2$theta[[i]]
		sign[sign != 0] <- 1
		common <- common + sign
	}
	common[common != length(Y)] <- 0
	common[common == length(Y)] <- 1
	common.g <- graph.adjacency(common, mode = "upper", weighted = TRUE)
	for(i in 1:length(Y)){
		adj <- -cov2cor(fit2$theta[[i]])  
		diag(adj) = 0
		gadj = graph.adjacency(adj, mode = "upper", weighted = TRUE)
		E(gadj)$weight = get.edge.attribute(gadj, "weight") * 15
		E(gadj)$color = 6
		ee <- get.edgelist(gadj)
		e0 <- get.edgelist(common.g)
		is.common <- which(apply(ee, 1, paste0, collapse="-") %in% apply(e0, 1, paste0, collapse="-"))
		E(gadj)$color[is.common] = 2
		# by hand!
		order <- c(7:11, 1:6, 12:27)
		if(i == 1)	layout <- layout_in_circle(gadj, order = V(gadj)[order])
		vlocs <- rep(1, length(V(gadj)))
		vlocs[20] = 0.8
		vlocs[21] = 1.1
		vlocs[22] = 0.9
		 radian.rescale <- function(x, start=0, direction=1) {
		   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
		   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
		 }
		locs <- radian.rescale(x=1:length(V(gadj)), direction=-1, start=0)[order]

		plot(gadj, vertex.frame.color = "white", layout = layout, 
		    vertex.size = 4, 
		    vertex.color= "darkred",
		    vertex.label.cex = 1.6, 
		    vertex.label.color= "black",
		    vertex.label.dist=vlocs,
		    vertex.label.degree=locs,
		    edge.width=abs(E(gadj)$weight), 
		    edge.color = E(gadj)$color, 
		    edge.lty = E(gadj)$lty)
		title(causes2[i],cex.main=2)

	}
	dev.off()
}
