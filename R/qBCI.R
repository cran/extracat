

qBCI <- function(x,...){
	UseMethod("qBCI")
}

qBCI.default <- function(x,y,p = NULL, k = 5,iter=20, ...){
	N <- length(x)
	
	if(is.null(p)){
		p <- sqrt(k/N)
		p <- 1/floor(1/p)
	}
	if(!is.factor(x)){
		qx <- quantile(x,seq(0,1,p))
		if( length(qxx <- unique(qx)) < length(qx)){
			simpleWarning("non-unique quantiles detected")	
		}
		if(length(qxx) == 2){
			x <- as.factor(x)	
		}else{
			x <- cut(x, unique(qxx) )
		}
	}
	if(!is.factor(y)){
		qy <- quantile(y,seq(0,1,p))
		if( length(qyy <- unique(qy)) < length(qy)){
			simpleWarning("non-unique quantiles detected")	
		}
		if(length(qyy) == 2){
			y <- as.factor(y)	
		}else{
			y <- cut(y, unique(qyy) )
		}
	}
	tt <- table(x,y)
	BCI(optile(tt,iter=iter))
}

qBCI.data.matrix <- function(x,p = NULL, k = 5, sort = TRUE, iter=20, ...){
	
	N <- nrow(x)
	nd <- ncol(x)
	
	
	if(is.null(p)){
		p <- (k/N)^(1/nd)
		p <- 1/floor(1/p)
	}
	
	x <- sapply(x, function(z) {
		cut(z, quantile(z,seq(0,1,p)))
		})
	x <- subtable(x,1:nd)
	x <- xtabs(Freq~.,data=x)
	if(sort){
		return(BCI(optile(x,iter=iter)))
	}else{
		return(BCI(x))
	}
	
}

qBCI.data.frame <- function(x,p = NULL, k = 5, sort = TRUE, iter=20, ...){
	x <- data.matrix(x)
	NextMethod("qBCI",x = x, p = p, k = k, sort = sort, iter = iter)
}


BCImat <- function(x, k = 5, iter = 20, p = NULL){

	nd <- ncol(x)
	ids <- combn(1:nd,2)
	
	values <- apply(ids,2,function(id){
		z <- na.omit(x[,id])
		qBCI(z[,1],z[,2], k = k, p = p, iter = iter)
	})
	M <- matrix(0,nd,nd)
	M[lower.tri(M)] <- values
	M <- M + t(M)
	colnames(M) <- rownames(M) <- names(x)
		return(M)
}

#BCImat2 <- function(x, k = 5, iter = 20, p = NULL){
	
#	nd <- ncol(x)
#	ids <- combn(1:nd,2)
	
#	types <- sapply(x, typeof)
#	if(any(types != "factor")){
#		vi <- which(types != "factor")
#		x[,vi] <- sapply(x[,vi], function(z) qbin(z, k = k, p = p))
#	}
	
	
#	values <- apply(ids,2,function(id){
#					tt <- table(x[,id[1]],x[,id[2]])
#					BCI(optile(tt, iter = iter))
#					})
#	M <- matrix(0,nd,nd)
#	M[lower.tri(M)] <- values
#	M <- M + t(M)
#	colnames(M) <- rownames(M) <- names(x)
#	return(M)
#}

wdcor.data.frame = function(x, approx = TRUE, ...){
	nmz <- names(x)
	x <- data.matrix(x)
	nd <- ncol(x)
	ids <- combn(1:nd,2)
	if(approx){
	values <- apply(ids,2,function(id){
		z <- na.omit(x[,id])
		approx.dcor(z[,1],z[,2], n=100)
	})
	}else{
		values <- apply(ids,2,function(id){
		z <- na.omit(x[,id])
		wdcor(z[,1],z[,2], n=100)
	})
	}
	M <- matrix(0,nd,nd)
	M[lower.tri(M)] <- values
	diag(M) <- 1
	M <- M + t(M)
	colnames(M) <- rownames(M) <- nmz
		return(M)
}

#mcrit <- function(x, fun = "BCI", fn = NULL, ff = NULL, all.numeric = FALSE, all.factor= FALSE, weights = "Freq", symmetric = TRUE , args = list()){
#	
#	if( weights %in% names(x) ){
#		w <- which(names(x) == weights)
#		x <- x[,-w]
#		w <- x[,w]	
#	}
#	nd <- ncol(x)
#	level <- 2
#	
#	
#	
##if(symmetric){
###	ids <- do.call(c,lapply(1:level,function(z) combn(1:nd,z,simplify=FALSE)))
##}else{
#		ids <- expand.grid(replicate(level,list(1:nd)))
##ids <- lapply(ids, unique)
##}
#	if(symmetric){
#		ids <- ids[which(ids[,1] >= ids[,2]), ]	
#	}
#	
#	isnum <- which(sapply(x, is.numeric))
#	
#	if(all.numeric){
#		x <- data.matrix(x)	
#		isnum <- 1:ncol(x)
#	}
#	if(all.factor){
#			x[,isnum] <- sapply(x[,isnum], qbin)
#	}
#	
#	types <- sapply(x,typeof)
#	
#	if( is.null(ff) ) ff <- fun
#	if( is.null(fn) ) fn <- fun
#	
#print(ids)	
#	values <- apply(ids,1, function(z){
#		fun1 <- fn
#		if( all (types[z] %in% c("numeric","integer") ) ){
#			fun1 <- fun	
#		}
#		if( all (types[z] %in% c("factor") ) ){
#			fun1 <- ff	
#		}
#		
#		cc <- call(fun1,x[,z])
#		for( i in seq_along(args) ){
#			 cc[[i+2]] <- args[[i]]
#		}
#		
#		eval(cc)
#	})
#	print(values)
#	if(symmetric){
#		M <- matrix(0,ncol=nd,nrow=nd)
#		M[lower.tri(M, diag = TRUE)] <- values	
#		M <- M + t(M)
#		diag(M) <- diag(M)/2
#	}else{
#		M <- values
#		dim(M) <- rep(nd, level)
#	}
#	return(M)
#	
#	
#}

qbin <- function(x,p = NULL, k = 5, d = 2){
	N <- length(x)
	if(is.null(p)){
		p <- (k/N)^(1/d)
		p <- 1/floor(1/p)
	}
	if(!is.factor(x)){
		qx <- quantile(x,seq(0,1,p), na.rm = TRUE)
		if( length(qxx <- unique(qx)) < length(qx)){
			simpleWarning("non-unique quantiles detected")	
		}
		if(length(qxx) == 2){
			x <- as.factor(x)	
		}else{
			x <- cut(x, unique(qxx) )
		}
	}
}
