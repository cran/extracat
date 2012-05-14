
dcor = function(x){
	n <- nrow(x)
	m <- ncol(x)
	
	cases <- expand.grid(1:n,1:m)
	DMY <- as.matrix(dist(cases[,1]))
	DMX <- as.matrix(dist(cases[,2]))
	
	pv <- as.vector(x)/sum(x)
	
	S1 <- sum(outer(pv,pv) * DMY*DMX)
	S2 <- sum(outer(pv,pv) * DMY)*sum(outer(pv,pv) * DMX)
		v1 <- DMY %*% pv
		v2 <- DMX %*% pv
	S3 <- t(v1*v2) %*% pv
	
	S1X <- sum(outer(pv,pv) * DMX*DMX)
	S2X <- sum(outer(pv,pv) * DMX)^2
	S3X <- t(v2^2) %*% pv
	
	S1Y <- sum(outer(pv,pv) * DMY*DMY)
	S2Y <- sum(outer(pv,pv) * DMY)^2
	S3Y <- t(v1^2) %*% pv
	
	dcor <- (S1+S2-2*S3)/sqrt(  (S1X+S2X-2*S3X)* (S1Y+S2Y-2*S3Y) )
	#scl <- sqrt(  (S1X+S2X-2*S3X)* (S1Y+S2Y-2*S3Y) )
	#return(c(dcor,S1,S2,S3,S1X,S1Y,scl))
	return(dcor)
}

dcor2.table = function(env,xi,yi){
	n <- length(yi)
	m <- length(xi)
	
	cases <- expand.grid(1:n,1:m)
	DMY <- as.matrix(dist(cases[,1]))
	DMX <- as.matrix(dist(cases[,2]))
	
	pv <- as.vector(env$mat[yi,xi])/sum(env$mat)
	
	S1 <- sum(outer(pv,pv) * DMY*DMX)
	S2 <- sum(outer(pv,pv) * DMY)*sum(outer(pv,pv) * DMX)
		v1 <- DMY %*% pv
		v2 <- DMX %*% pv
	S3 <- t(v1*v2) %*% pv
	
	S1X <- sum(outer(pv,pv) * DMX*DMX)
	S2X <- sum(outer(pv,pv) * DMX)^2
	S3X <- t(v2^2) %*% pv
	
	S1Y <- sum(outer(pv,pv) * DMY*DMY)
	S2Y <- sum(outer(pv,pv) * DMY)^2
	S3Y <- t(v1^2) %*% pv
	
	dcor <- (S1+S2-2*S3)/sqrt(  (S1X+S2X-2*S3X)* (S1Y+S2Y-2*S3Y) )
	#scl <- sqrt(  (S1X+S2X-2*S3X)* (S1Y+S2Y-2*S3Y) )
	#return(c(dcor,S1,S2,S3,S1X,S1Y,scl))
	return(dcor)
}

distcor = function(data, dims, perm.cat, ... ){
	nd <- length(dims)
	stopifnot(nd == 2)
	optimal <- FALSE
	globalbest <- -1
	currcrit <- -1
	bestcrit <- -1
	ci <- 1:dims[2]
	ri <- 1:dims[1]
	tri <- ri
	tci <- ci
	e1 <- new.env()
	e1$mat <- xtabs(Freq~.,data=data)
	#print(e1$mat)
	while(!optimal){
		
		#columns
		for(i in 1:dims[2]){
			for(j in 1:dims[2]){
				if(i < j){
					tci[i:(j-1)] <- tci[(i+1):j]
					tci[j] <- ci[i]
					currcrit <- dcor2.table(e1, tci, ri)
				}
				if(i > j){
					tci[(j+1):i] <- tci[j:(i-1)]
					tci[j] <- ci[i]
					currcrit <- dcor2.table(e1, tci, ri)
				}
				if(currcrit > bestcrit){
					bestcrit <- currcrit
					ci <- tci	
				}else{
					tci <- ci	
				}	
			}	
		}
		#rows
		for(i in 1:dims[1]){
			for(j in 1:dims[1]){
				if(i < j){
					tri[i:(j-1)] <- tri[(i+1):j]
					tri[j] <- ri[i]
					currcrit <- dcor2.table(e1, ci, tri)
				}
				if(i > j){
					tri[(j+1):i] <- tri[j:(i-1)]
					tri[j] <- ri[i]
					currcrit <- dcor2.table(e1, ci, tri)
				}
				if(currcrit > bestcrit){
					bestcrit <- currcrit
					ri <- tri	
				}else{
					tri <- ri	
				}	
			}	
		}
		if(bestcrit <= globalbest){
			optimal	<- TRUE
		}else{
			globalbest <- bestcrit	
		}
		#print(e1$mat[ri,ci])
		#print(globalbest)
	}
    if(kendalls(e1$mat[ri,ci])< 0){
        ri <- rev(ri)
    }
	#print(c(ri-1,ci-1,globalbest))	
	return(c(ri-1,ci-1,globalbest))
}

dcov2 = function(x){
	n <- nrow(x)
	stopifnot(ncol(x)==2)
	
	MY <- as.matrix(dist(x[,1]))
	MX <- as.matrix(dist(x[,2]))
	
	xc <- apply(MX,2,mean)
	xr <- apply(MX,1,mean)
	xm <- mean(MX)
	
	yc <- apply(MY,2,mean)
	yr <- apply(MY,1,mean)
	ym <- mean(MY)
	
	MX <- (MX - xr) - rep(xc, each=n) + xm
	MY <- (MY - yr) - rep(yc, each=n) + ym
	 dcov <- sum(MX*MY)/n^2
	return(dcov)	
}

# a test
#
#M <- bloma(3,3,1)
#dim(M)
#
#dcs <- replicate(1000,{ 
#	M2<- M[sample(1:nrow(M)),sample(1:ncol(M))]
#	dcor.table(M2)
#
#})
#dcor.table(optile(M,iter=40))
#
#
#
#tt <- as.table(M)
#
#dset <- as.data.frame(tt)
#dset <- subtable(dset,1:2)
#dset <- untableSet(dset)
#dset <- sapply(dset, as.integer)
#cor(dset[,1], dset[,2])
#
#
#tt <- as.table(optile(M,iter=40))
#
#dset <- as.data.frame(tt)
#dset <- subtable(dset,1:2)
#dset <- untableSet(dset)
#dset <- sapply(dset, as.integer)
#cor(dset[,1], dset[,2])
#
#
#cs1 <- replicate(1000,{ 
#
#dset <- as.data.frame(tt[sample(1:nrow(M)),sample(1:ncol(M))])
#dset <- subtable(dset,1:2)
#dset <- untableSet(dset)
#dset <- sapply(dset, as.integer)
#cor(dset[,1], dset[,2], method="spearman")
#
#})
#
