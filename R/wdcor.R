
dcorOld = function(x){
	n <- nrow(x)
	m <- ncol(x)
	
	cases <- expand.grid(1:n,1:m)
	DMY <- as.matrix(dist(cases[,1]))
	DMX <- as.matrix(dist(cases[,2]))
	
	pv <- as.vector(x)/sum(x)
	
	S1 <- pv %*% (DMX*DMY) %*% pv
	S2 <- (pv %*% DMX %*% pv) * (pv %*% DMY %*% pv)
		v1 <- DMY %*% pv
		v2 <- DMX %*% pv
	S3 <- t(v1*v2) %*% pv
	print(S1)
	print(S2)
	print(S3)
	S1X <- pv %*% (DMX*DMX) %*% pv
	S2X <- (pv %*% DMX %*% pv)^2
	S3X <- t(v2^2) %*% pv
	
	S1Y <- pv %*% (DMY*DMY) %*% pv
	S2Y <- (pv %*% DMY %*% pv)^2
	S3Y <- t(v1^2) %*% pv
	
	print(S1X)
	print(S2X)
	print(S3X)
	print(S1Y)
	print(S2Y)
	print(S3Y)
	
	dcor <- (S1+S2-2*S3)/sqrt(  (S1X+S2X-2*S3X) * (S1Y+S2Y-2*S3Y) )
	#scl <- sqrt(  (S1X+S2X-2*S3X)* (S1Y+S2Y-2*S3Y) )
	#return(c(dcor,S1,S2,S3,S1X,S1Y,scl))
	return(sqrt(dcor))
}
wdcor = function(x,...){
UseMethod("wdcor")
	
}

wdcor.default = function(x,y,w = NULL, ep = 1,...){
	
	if(is.null(w)){
		w <- rep(1,length(x))
	}
	stopifnot(length(x) == length(w))
	stopifnot(length(x) == length(y))
	stopifnot(all(w >= 0))
	
	storage.mode(w) <- "double"
	storage.mode(x) <- "double"
	storage.mode(y) <- "double"
	storage.mode(ep) <- "double"

ret <- .Call("dcorR",x,y,w/sum(w),ep)
	
	return(ret)
	
}
wdcor.table = function(x, ep = 1, ...){
	stopifnot( length(dim(x)) == 2) 
	dx <- as.data.frame(x)
	dx <- dx[dx$Freq > 0,]
	dx <- sapply(dx,as.numeric)
	
		NextMethod("wdcor", x = dx[,1], y = dx[,2], w = dx[,3], ep = ep)
}


approx.dcor = function(x,y, n = 50, correct = FALSE, ep = 1){
	stopifnot(length(x) == length(y))
	x <- cut(x,n)
	y <- cut(y,n)
	z <- table(x,y)
	
	ret <- wdcor(z, ep = ep)
	if(correct) ret <- ret*(1 + 2*sqrt(8)/n^2)
#if(correct) ret <- ret*(1 + 1/n/sqrt(12)/2)
#if(correct) ret <- ret+1.74/n^(15/8)
#if(correct) ret <- ret*(1+3/n^(1.84))
	return(ret)
}







opt.wdcor = function(env,xi,yi){
	return(wdcor(as.table(env$mat[yi,xi])))
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
					currcrit <- opt.wdcor(e1, tci, ri)
				}
				if(i > j){
					tci[(j+1):i] <- tci[j:(i-1)]
					tci[j] <- ci[i]
					currcrit <- opt.wdcor(e1, tci, ri)
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
					currcrit <- opt.wdcor(e1, ci, tri)
				}
				if(i > j){
					tri[(j+1):i] <- tri[j:(i-1)]
					tri[j] <- ri[i]
					currcrit <- opt.wdcor(e1, ci, tri)
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

#is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
#	abs(x - round(x)) < tol
#}
#	if(all(is.wholenumber(x))){
#		storage.mode(x) <- "integer"
#	}
