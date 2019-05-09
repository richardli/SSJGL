# This function is adapted from the JGL package (version 2.3, 2013-04-16)
# It takes lambda1 and lambda2 as matrices
# It also takes a lambda0 variable (scalar) for penalization of the diagonals
JGL.adaptive <-
function(Y, addvar = NULL, penalty="fused",lambda0,lambda1,lambda2,rho=1,weights="sample.size", maxiter=500,tol=1e-5,warm=NULL, warm.connected=NULL, return.whole.theta=FALSE, truncate=0, normalize=FALSE)
{
	## initialize:
	p = dim(Y[[1]])[2]
	K = length(Y)
	n = rep(0,K)
	for(k in 1:K) {n[k] = dim(Y[[k]])[1]}

	# assign feature names if none exist:
	if(length(dimnames(Y[[1]])[[2]])==0)
	{
		for(k in 1:K)
		{
			dimnames(Y[[k]])[[2]]=paste("V",1:p,sep="")
		}
	}
	
	# mean-normalize Y:
	if(normalize){
		for(k in 1:K){
		for(j in 1:p){
			Y[[k]][,j] = Y[[k]][,j]-mean(Y[[k]][,j])
		}}
	}	

	# set weights:
	if(length(weights)==1){if(weights == "equal"){
		weights = rep(1,K)
	}}
	if(length(weights)==1){if(weights == "sample.size"){
		weights = n #/sum(n)
	}}


	connected = rep(TRUE,p) 
	
	### now get criterion over connected S:
	## define S
	S = vector("list",length=K)
	for(k in 1:K)
	{
		ntemp = dim(Y[[k]])[1]
		S[[k]] = cov(Y[[k]][,connected])*(ntemp-1)/ntemp
		if(!is.null(addvar)) S[[k]] <- S[[k]] + addvar[[k]]/ntemp
		# if warm.connected specified, make those variance to be 0.
		if(!is.null(warm.connected)){
			dd <- diag(S[[k]])
			S[[k]][, !warm.connected] <- 0
			S[[k]][!warm.connected, ] <- 0
			diag(S[[k]]) <- dd
		}
	}

	# if a penalty matrix is entered, only take its appropriate rows:
	lam1 = lambda1
	lam2 = lambda2
	if(length(lam1)>1) {lam1 = lam1[connected,connected]}
	if(length(lam2)>1) {lam2 = lam2[connected,connected]}

	# ## examine criteria:  (value 0 where S allows theta=0 to satisfy KKT; value 1 where theta must be connected)
	if(penalty=="fused"){
  	if(K==2){  #use bi-conditional screening rule to identify block structure exactly
  	crit1 = list()
  	for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 + lam2 }  
  	S.sum = matrix(0,sum(connected),sum(connected))
  	for(k in 1:K) {S.sum = S.sum + weights[k]*S[[k]]}
  	S.sum = abs(S.sum)
  	crit2 = S.sum > 2*lam1
  	}
  	
    if(K>2){  #use sufficient screening rule to identify larger-grained block structure
  	crit1 = list()
  	for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 }  
  	crit2 = matrix(0,sum(connected),sum(connected))
  	}
  
  	# are both criteria met?
  	critboth = crit2
  	for(k in 1:K) {critboth = critboth + crit1[[k]]}
  	critboth = (critboth!=0)				
  	diag(critboth) = 1
	}
	
	if(penalty=="group")
	{
   	## examine criteria:  (value 0 where S allows theta=0 to satisfy KKT; value 1 where theta must be connected)
  	tempsum = matrix(0,sum(connected),sum(connected))
  	for(k in 1:K) {tempsum = tempsum + (pmax(weights[k]*abs(S[[k]]) - lam1,0))^2 }    
  	critboth = tempsum > lam2^2
  	diag(critboth) = 1
 	}
  
	## now identify block structure using igraph:
	g1 <- graph.adjacency(critboth)	
	cout = clusters(g1)
	blocklist = list()
	# identify unconnected elements, and get blocks:
	unconnected = c()

	# adapt cout$membership to start with index 1:
	if(min(cout$membership)==0){cout$membership=cout$membership+1}
	for(i in 1:(cout$no))
	{
		if(sum(cout$membership==i)==1) { unconnected <- c(unconnected,which(cout$membership==i)) }
		if(sum(cout$membership==i)>1) { blocklist[[length(blocklist)+1]] <- which(cout$membership==i) }
	}

	# final set of connected nodes
	connected[unconnected] = FALSE

	# connected indices of connected nodes:  0 for unconnected nodes, and 1:length(connected) for the rest.  
	# maps features 1:p to their order in the connected features
	connected.index = rep(0,p)
	connected.index[connected] = 1:sum(connected)
	# regular indices of connected nodes: map connected nodes onto 1:p indexing:

	# redefine unconnected as !connected (up until now it's been extra nodes caught as unconnected)
	unconnected=!connected

	# unconnected <- rep(FALSE, p)
	# connected <- rep(TRUE, p)

	## define theta on all connected:   (so theta is really theta.connected).
	theta = list()
	for(k in 1:K) 
	{
		theta[[k]] = matrix(0,sum(connected),sum(connected))
		if(sum(connected)>0)
		{
			dimnames(theta[[k]])[[1]]=dimnames(theta[[k]])[[2]]=dimnames(Y[[k]])[[2]][connected]	
		}
	}

	## get solution on unconnected nodes
	# data:
	Yu = list()
	for(k in 1:K){Yu[[k]] = Y[[k]][,unconnected]}
	# penalty vectors:
  	# note: for admm.iters.unconnected, we use the penalize.diagonal argument before calling the function.  for admm.iters, we use it IN the function.  
  	
  	#### This works fine for scalar, vector and matrix lambda1!
	#### if(length(lambda1)==1) { lam1.unconnected = lambda1 }
	#### if(length(lambda1)>1) { lam1.unconnected = diag(lambda1)[unconnected] }
	#### if(length(lambda2)==1) { lam2.unconnected = lambda2 }
	#### if(length(lambda2)>1) { lam2.unconnected = diag(lambda2)[unconnected] }
 	# if penalize.diagonal==FALSE, then set the appropriate penalty vectors to zero:
	#### if(!penalize.diagonal){ 
	#### 	lam1.unconnected = lam1.unconnected * 0
	#### 	if(penalty=="group") {lam2.unconnected = lam2.unconnected * 0}
	#### }
	# get the unconnected portion of theta:
	lam1.unconnected <- lambda0
	lam2.unconnected <- 0

	if(sum(unconnected)>0) 
	{
		theta.unconnected = JGL:::admm.iters.unconnected(Yu,lambda1=lam1.unconnected,lambda2=lam2.unconnected,penalty=penalty,rho=rho,weights=weights,maxiter=maxiter,tol=tol)$Z
		for(k in 1:K) { names(theta.unconnected[[k]])=dimnames(Y[[k]])[[2]][!connected] }
	}
	if(sum(unconnected)==0) {theta.unconnected = NULL}

	## now run JGL on each block of the connected nodes to fill in theta:
	if(length(blocklist)>0){
	for(i in 1:length(blocklist)){
		# the variables in the block
		bl <- blocklist[[i]] 
		Ybl = list()
		# get the data on only those variables
		for(k in 1:K) 
		{
			Ybl[[k]] = Y[[k]][,bl]
		}  
		# penalty matrices:
		if(length(lambda1)==1) { lam1.bl = lambda1 }
		if(length(lambda1)>1) { lam1.bl = lambda1[bl,bl] }
		if(length(lambda2)==1) { lam2.bl = lambda2 }
		if(length(lambda2)>1) { lam2.bl = lambda2[bl,bl] }
	  	# initialize lambdas:
 	 	lam1.bl = JGL:::penalty.as.matrix(lam1.bl,dim(Ybl[[1]])[2],penalize.diagonal=TRUE)
 	 	diag(lam1.bl) <- lambda0

		if(penalty=="fused") {lam2.bl = JGL:::penalty.as.matrix(lam2.bl,dim(Ybl[[1]])[2],penalize.diagonal=FALSE)}
		if(penalty=="group") {lam2.bl = JGL:::penalty.as.matrix(lam2.bl,dim(Ybl[[1]])[2],penalize.diagonal=FALSE)}
    
		# implement warm start if desired
		if(length(warm)==0) {warm.bl = NULL}
		if(length(warm)>0)
		{
			warm.bl = list()
			for(k in 1:K) { warm.bl[[k]] = warm[[k]][bl,bl] }
		}
		# run JGL on the block:
		Thetabl = admm.iters.adaptive(Ybl,lam1.bl,lam2.bl,penalty=penalty,rho=rho,weights=weights,maxiter=maxiter,tol=tol,warm=warm.bl)
		# update Theta with Thetabl's results:
		for(k in 1:K) {theta[[k]][connected.index[bl],connected.index[bl]] = Thetabl$Z[[k]]}   
	}}

	# round very small theta entries down to zero:
	if(dim(theta[[1]])[1]>0)
	{
	for(k in 1:K)
	{
		rounddown = abs(theta[[k]])<truncate; diag(rounddown)=FALSE
		theta[[k]]=theta[[k]]*(1-rounddown)
	}}

	# return output: theta on connected nodes, diagonal theta on unconnected nodes, and the identities of the connected nodes
	if(!return.whole.theta) 
	{
		out = list(theta=theta,theta.unconnected=theta.unconnected,connected=connected)
	}
	if(return.whole.theta) 
	{
		whole.theta = list()
		for(k in 1:K) 
		{
			whole.theta[[k]] = matrix(0,p,p)
			if(sum(unconnected) > 0) diag(whole.theta[[k]])[unconnected] = theta.unconnected[[k]]
			whole.theta[[k]][connected,connected] = theta[[k]]
			dimnames(whole.theta[[k]])[[1]] = dimnames(whole.theta[[k]])[[2]] = dimnames(Y[[k]])[[2]]
		}
		out = list(theta=whole.theta,connected=connected)
	}
	class(out)="jgl"
	return(out)
}

