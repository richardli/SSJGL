path.plot <- function(v0s, obj, G, thres = 0.5, normalize = FALSE, xlab="", ylab = "", main = "", par = c(1, 2), ylim = NULL, reverse=FALSE, position="bottomright", only = NULL, color0 = "gray70", color1 = "darkblue", which.top = c(0,1)[2], cex=.5, cex.xlab=1, cex.ylab=1, vline = NULL, ...){
	NV <- length(obj$thetalist)
	K <- length(obj$thetalist[[1]])
	P <- dim(obj$thetalist[[1]][[1]])[1]
	if(is.null(obj$problist1)) obj$problist1 <- obj$problist	
	if(is.null(only)) only <- c(0, 1)
	if(!is.null(par)) par(mfrow = par)
	
	if(is.null(ylim)){
	    for(k in 1:K){
			prob <- omega <- array(NA, dim = c(P, P, NV))
			for(i in 1:NV){ 
				prob[, , i] <- obj$problist1[[i]]
				if(normalize){
					omega[, , i] <- -cov2cor(obj$thetalist[[i]][[k]]) 
				}else{
					omega[, , i] <- (obj$thetalist[[i]][[k]])
				}
			}
			if(!is.null(prob)){
				for(i in 1:dim(omega)[3]){ 
					omega[, , i] <- omega[, , i] * (prob[, , i] >= thres)
				}
			}
			for(i in 1:dim(omega)[3]){ 
				diag(omega[, , i]) <- NA 
			}
			ylim <- range(c(ylim, omega), na.rm=TRUE)
		}
	}

	for(k in 1:K){
		prob <- omega <- array(NA, dim = c(P, P, NV))
		for(i in 1:NV){ 
			prob[, , i] <- obj$problist1[[i]]
			if(normalize){
				omega[, , i] <- -cov2cor(obj$thetalist[[i]][[k]]) 
			}else{
				omega[, , i] <- (obj$thetalist[[i]][[k]])
			}
		}
		if(!is.null(prob)){
			for(i in 1:dim(omega)[3]){ 
				omega[, , i] <- omega[, , i] * (prob[, , i] >= thres)
			}
		}
		for(i in 1:dim(omega)[3]){ 
			diag(omega[, , i]) <- NA 
		}
		xlim <- range(v0s)
		if(reverse) xlim = rev(range(v0s))
		if(is.null(ylim)) ylim <- range(omega, na.rm=T)	

		plot(v0s, rep(0, length(v0s)), ylim = ylim, col = "white", xlab = xlab, main = main[k], ylab = ylab, xlim = xlim, ...)
		if(!is.null(vline)) abline(v = v0s[vline], col = "darkblue", lty = 2)

		M <- dim(G[[k]])[1]
	
		if(1 %in% only && which.top == 0){
			for(i in 2:M){
				for(j in 1:(i-1)){
					if(G[[k]][i, j] != 0){
						lines(v0s, omega[i, j, ], type = "b", col = color1, cex=cex, pch = 20)
					}
				}
			}
		}
		if(0 %in% only){	
			for(i in 2:M){
				for(j in 1:(i-1)){
					if(G[[k]][i, j] == 0 && sum(abs(omega[i,j,])) > 0){
					# if(G[[k]][i, j] == 0 ){
						lines(v0s, omega[i, j, ], type = "b", col = color0, lty = 4, cex=cex, pch = 20)
					}
				}
			}
		}
		if(1 %in% only && which.top == 1){
			for(i in 2:M){
				for(j in 1:(i-1)){
					if(G[[k]][i, j] != 0){
						lines(v0s, omega[i, j, ], type = "b", col = color1, cex=cex, pch = 20)
					}
				}
			}
		}
		if(0 %in% only && 1 %in% only){
			if(position != "none") legend(position, lty = c(1, 4), pch= c(20, 20), col = c(color1, color0), c("Edge", "Non-edge"))
		}
	}
}



curveplot <- function(fitlist, xobj, yobj, type='l', xlab="", ylab = "", main="", grid=100, color="black", add = FALSE, xlim = NULL, ylim = NULL, lty=1, which.density=NULL, lastcol=FALSE, ggplot=TRUE, labels=NULL, dotcolor="blue", dotalpha=0.1, dotsize = 0.1, which.ss = NULL, extrap = TRUE, ...){
	require("Hmisc")
	N <- length(fitlist)
	xm <- matrix(NA, N, grid)
	ym <- matrix(NA, N, grid)
	update_xlim <- is.null(xlim)
	update_ylim <- is.null(ylim)
	dotdat <- NULL
	dotcounter <- 1
	for(rep in 1:N){
		fit <- fitlist[[rep]]
		L <- length(fit)
		if(L == 0) next
		dim <- c(1, 1)
		if(class(fit[[1]][[xobj]]) == "matrix"){
			dim[1] <- 2
			M <- dim(fit[[1]][[xobj]])[1]
		}else{
			dim[1] <- 1
			M <- length(fit[[1]][[xobj]])
		}
		if(class(fit[[1]][[yobj]]) == "matrix"){
			dim[2] <- 2
		}else{
			dim[2] <- 1
		}
		xx <- yy <- matrix(NA, L, M)
		if(dim[1] == 1){
			for(i in 1:L){
				xx[i, ] <-  fit[[i]][[xobj]]
			}		
		}else{
			for(i in 1:L){
				xx[i, ] <-  apply(fit[[i]][[xobj]], 1, sum)
			}
		}
		if(dim[2] == 1){
			for(i in 1:L){
				yy[i, ] <-  fit[[i]][[yobj]]
			}		
		}else{
			for(i in 1:L){
				yy[i, ] <-  apply(fit[[i]][[yobj]], 1, sum)
			}
		}
		if(!is.null(which.density)){
			if(rep %in% which.density){
				if(lastcol){
				dat0 <- data.frame(x=as.numeric(xx[, dim(xx)[2]]), y = as.numeric(yy[, dim(yy)[2]]))
				}else{
				dat0 <- data.frame(x=as.numeric(xx), y = as.numeric(yy))	
				}
				if(!ggplot) kk <- with(dat0,MASS:::kde2d(x,y))
				dat0$color <- dotcolor[dotcounter]
				dotcounter <- dotcounter + 1
				dotdat <- rbind(dotdat, dat0)
			}
		}

		if(rep %in% which.ss){
			xx <- xx[, -1]
			yy <- yy[, -1]
		}
		min <- min(xx)
		max <- max(xx)
		yy2 <- matrix(NA, L, grid)
		xout <- seq(min, max, len=grid)
		for(i in 1:L){
			if(length(unique(xx[i,])) == 1) next
			if(extrap){
				tmp <- approxExtrap(xx[i, ], yy[i, ], xout = xout)
				tmp$y[xout > max(xx[i, ])] <- NA

			}else{
				tmp <- approx(xx[i, ], yy[i, ], xout = xout)
			}
			yy2[i, ] <- tmp$y
			# out <- union(which(xout > max(xx[i,])), which(xout < min(xx[i,])))
			# yy2[i, out] <- NA
		}
		if(sum(!is.na(yy2)) == 0) next

		xm[rep, ] <- xout
		ym[rep, ] <- apply(yy2, 2, mean, na.rm=TRUE)
		# toosmall <- apply(yy2, 2, function(x){sum(!is.na(x)) / length(x)})
		# ym[rep, toosmall > 0.9] <- NA

		if(update_xlim) xlim <- range(c(xlim, xm[rep, ]), na.rm=T)
		if(update_ylim) {
			shown <- intersect(which(xm[rep, ] <= xlim[2]), which(xm[rep, ] >= xlim[1]))
			if(extrap){
				boundary <- approxExtrap(xm[rep, ], ym[rep, ], xout =xlim)$y
			}else{
				boundary <- approx(xm[rep, ], ym[rep, ], xout =xlim)$y
			}
			ylim <- range(c(ylim, ym[rep, shown], boundary), na.rm=T)
		}
	}
	
	if(!ggplot){
		if(!is.null(which.density)){
			filled.contour.nol(kk, xlim = xlim, ylim = ylim, xlab=xlab, ylab=ylab, main = main, col=colorRampPalette(brewer.pal(9, "Blues"))(20))
		}
		lty <- as.numeric(lty)
		color <- as.character(color)

		if(!add && is.null(which.density)){
			plot(xm[1,], ym[1,], xlim = xlim, ylim = ylim, type="l", col = "white", xlab=xlab, ylab=ylab, main = main, lty=lty, ...)
		}
		for(i in 1:N){
			if(!is.null(which.density) && i == which.density){
				next
			}
			lines(xm[i,], ym[i,], type="l", col = color[i], lty=lty[i], ...)
		}
	}else{
		require(reshape)
		if(!is.null(which.density)){
			ym[which.density, ] <- NA
			xm[which.density, ] <- NA
		}
		dat <- t(ym)
		colnames(dat) <- labels
		dat <- melt(dat)
		dat[,1] <- as.numeric(t(xm))
		dat$lty <- rep(as.character(lty), each=dim(ym)[2])
		dat$col <- rep(as.character(color), each=dim(ym)[2])
		colnames(dat) <- c("x", "class", "y", "lty", "color")
		g <- ggplot()

		g <- g + geom_line(data=dat, aes(x=x, y=y, group=class, color=color, linetype=lty)) + xlab(xlab) + ylab(ylab) + xlim(xlim) 
		if(!is.null(which.density)){
			# g <- g + stat_density2d(data=dat0, aes(x=x, y = y, fill=..level..), geom="polygon")
			g <- g + geom_point(data=dotdat, aes(x=x, y = y, color=color), alpha=dotalpha, size=dotsize)
		}
		return(g)
	}
}





# grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

#   plots <- list(...)
#   position <- match.arg(position)
#   g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
#   legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#   lheight <- sum(legend$height)
#   lwidth <- sum(legend$width)
#   gl <- lapply(plots, function(x) x + theme(legend.position="none"))
#   gl <- c(gl, ncol = ncol, nrow = nrow)

#   combined <- switch(position,
#                      "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
#                                             legend,
#                                             ncol = 1,
#                                             heights = unit.c(unit(1, "npc") - lheight, lheight)),
#                      "right" = arrangeGrob(do.call(arrangeGrob, gl),
#                                            legend,
#                                            ncol = 2,
#                                            widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

#   grid.newpage()
#   grid.draw(combined)

#   # return gtable invisibly
#   invisible(combined)

# }

densityplot <- function(fit, xobj, yobj, type='l', xlab="", ylab = "", main="", grid=100, color="black", add = FALSE, xlim = NULL, ylim = NULL, lty=1, palette="Blues", lastcol=FALSE, ...){
	update_xlim <- is.null(xlim)
	update_ylim <- is.null(ylim)
	L <- length(fit)
	if(L == 0) next
	dim <- c(1, 1)
	if(class(fit[[1]][[xobj]]) == "matrix"){
		dim[1] <- 2
		M <- dim(fit[[1]][[xobj]])[1]
	}else{
		dim[1] <- 1
		M <- length(fit[[1]][[xobj]])
	}
	if(class(fit[[1]][[yobj]]) == "matrix"){
		dim[2] <- 2
	}else{
		dim[2] <- 1
	}
	xx <- yy <- matrix(NA, L, M)
	if(dim[1] == 1){
		for(i in 1:L){
			xx[i, ] <-  fit[[i]][[xobj]]
		}		
	}else{
		for(i in 1:L){
			xx[i, ] <-  apply(fit[[i]][[xobj]], 1, sum)
		}
	}
	if(dim[2] == 1){
		for(i in 1:L){
			yy[i, ] <-  fit[[i]][[yobj]]
		}		
	}else{
		for(i in 1:L){
			yy[i, ] <-  apply(fit[[i]][[yobj]], 1, sum)
		}
	}
	if(lastcol){
		xx = xx[, dim(xx)[2]]
		yy = yy[, dim(yy)[2]]
	}
	dat <- data.frame(x=as.numeric(xx), y = as.numeric(yy))
 	# k <- with(dat,MASS:::kde2d(x,y))

 	if(update_xlim) xlim <- range(c(xlim, range(xx)), na.rm=T)
	if(update_ylim) ylim <- range(c(ylim, range(yy)), na.rm=T)

	# filled.contour.nol(k, xlim = xlim, ylim = ylim, xlab=xlab, ylab=ylab, main = main, col=colorRampPalette(brewer.pal(9, palette))(20))
	g <- ggplot(dat, aes(x=x, y=y))+ stat_density_2d(geom = "raster", aes(fill = ..density..), contour=FALSE)
	return(g)
}




densityplot_multi <- function(fitlist, xobj, yobj, type='l', xlab="", ylab = "", main="", grid=100, color="black", add = FALSE, xlim = NULL, ylim = NULL, lty=1, palette="Blues", lastcol=FALSE, labels = NULL, ...){
	update_xlim <- is.null(xlim)
	update_ylim <- is.null(ylim)
	if(is.null(labels)) labels <- as.character(1:length(fitlist))
	dat <- NULL
	if(length(lastcol) == 1) lastcol = rep(lastcol, length(fitlist))
	for(rep in 1:length(fitlist)){
		fit <- fitlist[[rep]]
		L <- length(fit)
		if(L == 0) next
		dim <- c(1, 1)
		if(class(fit[[1]][[xobj]]) == "matrix"){
			dim[1] <- 2
			M <- dim(fit[[1]][[xobj]])[1]
		}else{
			dim[1] <- 1
			M <- length(fit[[1]][[xobj]])
		}
		if(class(fit[[1]][[yobj]]) == "matrix"){
			dim[2] <- 2
		}else{
			dim[2] <- 1
		}
		xx <- yy <- matrix(NA, L, M)
		if(dim[1] == 1){
			for(i in 1:L){
				xx[i, ] <-  fit[[i]][[xobj]]
			}		
		}else{
			for(i in 1:L){
				xx[i, ] <-  apply(fit[[i]][[xobj]], 1, sum)
			}
		}
		if(dim[2] == 1){
			for(i in 1:L){
				yy[i, ] <-  fit[[i]][[yobj]]
			}		
		}else{
			for(i in 1:L){
				yy[i, ] <-  apply(fit[[i]][[yobj]], 1, sum)
			}
		}
		if(lastcol[rep]){
			xx = xx[, dim(xx)[2]]
			yy = yy[, dim(yy)[2]]
		}
		dat <- rbind(dat, data.frame(x=as.numeric(xx), y = as.numeric(yy), group=labels[rep]))
	 	if(update_xlim) xlim <- range(c(xlim, range(xx)), na.rm=T)
		if(update_ylim) ylim <- range(c(ylim, range(yy)), na.rm=T)
	}
 # 	g <- ggplot() + xlim(xlim) + ylim(ylim)
	# for(sub in labels){
	# 	g <- g + stat_density2d(geom="tile", aes(x=x, y=y, fill = group, alpha=..density..), contour=FALSE, data=subset(dat, group==sub))
	# }
	g <- ggplot(dat, aes(x=x, y=y)) + xlim(xlim) + ylim(ylim) + stat_density2d(geom="tile", aes(fill = group), contour=FALSE)
   #                + 
  	# scale_fill_manual(values=c("a"="#FF0000", "b"="#00FF00"))
	return(g)
}

filled.contour.nol <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
    col = colorRampPalette(length(levels) - 1), plot.title, plot.axes, 
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, main = "", xlab=xlab, ylab = ylab, ...) 
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    # if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    #     stop("increasing 'x' and 'y' values expected")
    # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    # on.exit(par(par.orig))
    # w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    # layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    # par(las = las)
    # mar <- mar.orig
    # mar[4L] <- mar[2L]
    # mar[2L] <- 1
    # par(mar = mar)
    # plot.new()
    # plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
    #     yaxs = "i")
    # rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    # if (missing(key.axes)) {
    #     if (axes) 
    #         axis(4)
    # }
    # else key.axes
    # box()
    # if (!missing(key.title)) 
    #     key.title
    # mar <- mar.orig
    # mar[4L] <- 1
    # par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    # plot(NA, xlim=xlim, ylim=ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    .filled.contour(x, y, z, levels, col)
    # if (missing(plot.axes)) {
        # if (axes) {
            title(main = main, xlab = xlab, ylab = ylab)
            Axis(x, side = 1)
            Axis(y, side = 2)
        # }
    # }
    # else plot.axes
    # if (frame.plot) 
        # box()
    # if (missing(plot.title)) 
        # title(...)
    # else plot.title
    # invisible()
}


getmetric <- function(est, truth, graph){
	out <- NULL
	gtmp <- graph
	tmp <- est
	tmp0 <- truth
	diag(tmp) <- diag(tmp0) <- diag(gtmp) <- NA
	out$SSE<- sum((tmp-tmp0)^2, na.rm=T)
	out$L1 <- sum(abs(tmp), na.rm=T)

	tmp[tmp!=0] <- 1
	tmp0[tmp0!=0] <- 1
	diag(tmp) <- diag(tmp0) <- diag(gtmp) <- NA
	out$nedges <- sum(tmp != 0, na.rm=T)/2
	out$tp <- sum(tmp * gtmp > 0, na.rm=T)/2
	out$fp <- sum(tmp * (1-gtmp) > 0, na.rm=T)/2
	out$fn <- sum((1-tmp) * gtmp > 0, na.rm=T)/2
	out$tn <- sum((1-tmp) * (1-gtmp) > 0, na.rm=T)/2
	out$dKL <- 0.5 * (-log(det(est %*% solve(truth))) + sum(diag(est %*% solve(truth))))
	return(out)
}


getmetric2 <- function(est, truth, graph){
	out <- NULL
	gtmp <- graph
	tmp <- est
	tmp0 <- truth

	out$mNorm <- norm(as.matrix(tmp) - as.matrix(tmp0), type = "m")
	out$sNorm <- base::norm(as.matrix(tmp) - as.matrix(tmp0), type = "2")
	out$fNorm <- norm(as.matrix(tmp) - as.matrix(tmp0), type = "f")

	tmp <- cov2cor(tmp)
	tmp0 <- cov2cor(tmp0)
	 
	out$mNorm2 <- norm(as.matrix(tmp) - as.matrix(tmp0), type = "m")
	out$sNorm2 <- base::norm(as.matrix(tmp) - as.matrix(tmp0), type = "2")
	out$fNorm2 <- norm(as.matrix(tmp) - as.matrix(tmp0), type = "f")

	return(out)
}

getdiffmetric <- function(est, truth, graph, tol){
	out <- NULL
	out$tp.diff=out$fp.diff=out$fn.diff=0
	out$tp.gdiff=out$fp.gdiff=out$fn.gdiff=0
	G <- length(est)
	for(i in 1:(G-1)){
		for(j in (i+1) : G){
				tmp.a <- est[[i]]
				tmp.a[tmp.a!=0] <- 1
				gtmp.a <- graph[[i]]
				tmp.b <- est[[j]]
				tmp.b[tmp.b!=0] <- 1
				gtmp.b <- graph[[j]]
				diag(tmp.a) <- diag(tmp.b) <- diag(gtmp.a) <- diag(gtmp.b) <- NA

				out$tp.gdiff <- out$tp.gdiff+sum((tmp.a != tmp.b) * (gtmp.a != gtmp.b), na.rm=T)/2
				out$fp.gdiff <- out$fp.gdiff+sum((tmp.a != tmp.b) * (gtmp.a == gtmp.b), na.rm=T)/2
				out$fn.gdiff <- out$fn.gdiff+sum((tmp.a == tmp.b) * (gtmp.a != gtmp.b), na.rm=T)/2
		}
	}
	for(i in 1:(G-1)){
		for(j in (i+1) : G){
				tmp.a <- est[[i]]
				gtmp.a <- truth[[i]]
				tmp.b <- est[[j]]
				gtmp.b <- truth[[j]]
				diag(tmp.a) <- diag(tmp.b) <- diag(gtmp.a) <- diag(gtmp.b) <- NA
				diff <- abs(tmp.a - tmp.b) > tol 
				tdiff <- abs(gtmp.a - gtmp.b) > tol

				out$tp.diff <- out$tp.diff+sum(diff * tdiff, na.rm=T)/2
				out$fp.diff <- out$fp.diff+sum(diff * (1-tdiff), na.rm=T)/2
				out$fn.diff <- out$fn.diff+sum((1-diff) * tdiff, na.rm=T)/2
		}
	}
	return(out)
}




bdgraph.sim.rho <- function (p = 10, graph = "random", n = 0, type = "Gaussian", 
    prob = 0.2, size = NULL, mean = 0, class = NULL, cut = 4, 
    b = 3, D = diag(p), K = NULL, sigma = NULL, vis = FALSE, rho = 0.7) 
{
    if (is.matrix(K)) 
        graph <- "fixed"
    if (type == "normal") 
        type = "Gaussian"
    if (type == "non-normal") 
        type = "non-Gaussian"
    if (is.matrix(graph)) {
        G <- graph
        graph <- "fixed"
    }
    if (sum(graph != c("fixed", "AR1", "AR2", "star")) == 0) {
        G <- BDgraph::graph.sim(p = p, graph = graph, prob = prob, 
            size = size, class = class, vis = vis)
    }
    if (graph == "AR1") {
        sigma = matrix(0, p, p)
        for (i in 1:(p - 1)) for (j in (i + 1):p) sigma[i, j] = (rho)^abs(i - 
            j)
        sigma = sigma + t(sigma) + diag(p)
        K = solve(sigma)
    }
    if (graph == "AR2") {
        K = toeplitz(c(1, 0.5, 0.25, rep(0, p - 3)))
    }
    if (graph == "star") {
        K <- diag(p)
        K[1, (2:p)] <- 0.1
        K[(2:p), 1] <- 0.1
    }
    if (n != 0) {
        if (!is.null(sigma)) 
            K <- solve(sigma)
        if (is.matrix(K)) {
            G <- 1 * (abs(K) > 0.02)
            diag(G) <- 0
            if (is.null(sigma)) 
                sigma <- solve(K)
        }
        else {
            Ti = chol(solve(D))
            diag(G) = 0
            K = matrix(0, p, p)
            result = .C("rgwish_c", as.integer(G), as.double(Ti), 
                K = as.double(K), as.integer(b), as.integer(p), 
                PACKAGE = "BDgraph")
            K = matrix(result$K, p, p)
            sigma = solve(K)
        }
        d <- BDgraph::rmvnorm(n = n, mean = mean, sigma = sigma)
        if (type == "mixed") {
            ps = floor(p/5)
            col_number <- c(1:ps)
            prob <- pnorm(d[, col_number])
            d[, col_number] <- qpois(p = prob, lambda = 10)
            col_number <- c((ps + 1):(2 * ps))
            prob <- pnorm(d[, col_number])
            d[, col_number] <- qpois(p = prob, lambda = 2)
            col_number <- c((2 * ps + 1):(3 * ps))
            prob <- pnorm(d[, col_number])
            d[, col_number] <- qexp(p = prob, rate = 10)
            col_number <- c((3 * ps + 1):(4 * ps))
            prob <- pnorm(d[, col_number])
            d[, col_number] <- qbinom(p = prob, size = 1, prob = 0.5)
        }
        if (type == "non-Gaussian") {
            prob <- pnorm(d)
            d <- qexp(p = prob, rate = 10)
        }
        if (type == "discrete") {
            runif_m <- matrix(runif(cut * p), nrow = p, ncol = cut)
            marginals <- apply(runif_m, 1, function(x) {
                qnorm(cumsum(x/sum(x))[-length(x)])
            })
            if (cut == 2) 
                marginals = matrix(marginals, nrow = 1, ncol = p)
            for (j in 1:p) {
                breaks <- c(min(d[, j]) - 1, marginals[, j], 
                  max(d[, j]) + 1)
                d[, j] <- as.integer(cut(d[, j], breaks = breaks, 
                  right = FALSE))
            }
            d = d - 1
        }
        if (type == "binary") {
            if (p > 16) 
                stop("For type 'binary', number of nodes (p) must be less than 16")
            clique_factors = generate_clique_factors(ug = G)
            d = sample_ug(n = n, ug = G, clique_factors = clique_factors)
            d = d - 1
        }
    }
    if (n != 0) {
        simulation <- list(G = G, data = d, sigma = sigma, K = K, 
            graph = graph, type = type)
    }
    else {
        simulation <- list(G = G, graph = graph)
    }
    class(simulation) <- "sim"
    return(simulation)
}
