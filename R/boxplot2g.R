
boxplot2g = function(x,y=NULL, groups = NULL, smooth = loess, smooth.args = list(span = 0.1),  colv = NULL, alpha = 1, n = 360,...){

prbs <- c(0.25,0.5,0.75)

if(is.null(y)){
	stopifnot(ncol(x)==2)	
	data <- as.data.frame(x)
}else{
	data <- as.data.frame(cbind(x,y))	
}

if(is.null(groups)){
	data$groups <- as.factor(0)
}else{
	data$groups <- as.factor(groups)
}

labs <- names(data)
names(data) <- c("x","y","groups")
DM <- data.matrix(data)
#require(ggplot2)


# initiate the smoother
if(is.logical(smooth)){
	do.smooth <- smooth	
}else{
	do.smooth <- TRUE	
}
if(do.smooth){
	poss.args <- names(formals(smooth))
	spec.args <- names(smooth.args)


	ind <- match(spec.args, poss.args)
	for(i in seq_along(ind)){
		formals(smooth)[ind[i]] <- smooth.args[[i]]	
	}	
	if("span" %in% poss.args){
		formals(smooth)$span <- formals(smooth)$span/3
	}
}else{
	smooth <- NULL
}

	phi = seq(360/n, 360, 360/n)/180*pi
	e1 <- new.env()
	e1$vectors <- cbind(sin(phi),cos(phi))

ntv <- nlevels(data$groups)
if(is.null(colv)){
#print(ntv)
	if(ntv == 1){
		colv = 1	
	}else{
		colv <- rainbow(ntv)	
	}
}

e1$colv <- colv
e1$lvls <- levels(data$groups)
#colv <- colv[match(groups,levels(as.factor(data$groups)))]
#e1$gp <- qplot(data$x, data$y, colour = data$groups)	
e1$gp <- ggplot(data=data,aes(x=x,y=y,colour=groups))+geom_point(alpha=alpha)	
#print(formals(smooth))
if(ntv == 1){
	groupbox2d(x=data,env=e1,prbs=prbs,smooth=smooth,do.smooth)
}else{
	by(data,groups, groupbox2d, env= e1, prbs = prbs, smooth = smooth)
}
#e1$gp <- e1$gp  + opts(legend.position = "none") 
return(e1$gp)
	
}


groupbox2d = function( x, env, prbs, past, smooth){
	
	grp <- x[1,3]	
	
	colid <- match(grp, env$lvls)
	
	if(any(colid)){
		colv <- env$colv[]
		}else{
		colv <- env$col[1]	
		}
	xs <- x[,1:2]
	mm <- apply(xs,2,mean)
	xs <-  data.matrix(xs) - rep(mm,each=nrow(xs))
	
	S <- cov(xs)
	
	if (requireNamespace("MASS", quietly = TRUE)) {
		Sinv <- MASS::ginv(S)
		SSinv <- svd(Sinv)
		SSinv <- SSinv$u %*% diag(sqrt(SSinv$d))
		SS <- MASS::ginv(SSinv)
	}else{
		Sinv <- solve(S)
		SSinv <- svd(Sinv)
		SSinv <- SSinv$u %*% diag(sqrt(SSinv$d))
		SS <- solve(SSinv)	
	}

	xs <- xs %*% SSinv
	
	
	prj <- xs %*% t(env$vectors)
	
	qut <- t(apply(prj,2, function(z){
		quarts <- quantile(z, probs = prbs)
		iqr <- quarts[3]-quarts[1]
		w1 <- min(z[which(z >= quarts[1] - 1.5*iqr)])
		#w2 <- max(z[which(z <= quarts[3] + 1.5*iqr)])
		#return(c(w1,quarts,w2))
		return(c(w1,quarts))
	}))
	#print(formals(smooth))
	if( !is.null(smooth) ){
		n <- nrow(qut)
		qut <- apply(qut,2,function(z){
			x <- 1:(3*n)
			z <- rep(z,3)
			ys <- predict(smooth(z~x))
			return(ys[(n+1):(2*n)])
		})
		#print(dim(qut))
	}
	
	
	ccBox <- env$vectors*qut[,2]
	md <- data.frame((env$vectors*qut[,3])%*%SS)
	md <- sapply(md,mean)+mm		
	
	md[3] <- grp
	
	ccWsk <- env$vectors*qut[,1]
	
	ccc <- data.frame(rbind(ccBox,ccWsk) %*% SS + rep(mm,each=2*nrow(ccBox)))
	
	ccc$grp <- as.factor(rep(c("box","wsk"),each=nrow(ccBox)))
	ccc$groups <- factor(grp)
	md <- data.frame(md[1],md[2],grp)
	names(md) <- names(ccc)[-3]
	X1 <- NULL
	X2 <- NULL
	groups <- NULL
	#env$gp <- env$gp + geom_point(x=md[1],y=md[2],colour=md[3])
	env$gp <- env$gp + geom_point(data=md,aes(x=X1,y=X2, colour = groups),size=5) +  scale_colour_manual(values = colv) 
	
	env$gp <- env$gp + geom_path(data=ccc, aes(x=X1,y=X2,group=grp, colour = groups), alpha = 1/8)
	env$gp <- env$gp + geom_polygon(data=ccc,aes(x=X1,y=X2,group=grp, colour = groups, fill = groups), alpha = 1/8)
	env$gp <- env$gp + geom_point(data=md,aes(x=X1,y=X2),size=3,alpha=1,colour="white")
	env$gp <- env$gp + geom_point(data=md,aes(x=X1,y=X2),size=1,alpha=1)
return( invisible(TRUE) )

}
