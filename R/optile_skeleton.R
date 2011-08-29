
optile = function(x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
freqvar = NULL, return.data = TRUE, return.type = "data.frame", vs = 0, ...){
	UseMethod("optile")	
}

optile.default = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "data.frame", vs = 0, ...){

	
stepwise <- ifelse(is.null(method), FALSE, method %in% c("stepwise","sw","step"))
joint <- ifelse(is.null(method), FALSE, method %in% c("joint","pairwise","pw"))
# mv <- method %in% c("mv","multivar","multivariate")

fi <- which(names(x) %in% c("Freq",freqvar))
nd <- ncol(x)-any(fi)

# workarounds:
if( fun == "rmca" & nd == 2 ){
	fun <- "casort"	
}
if( fun == "casort" ){
	joint <- TRUE	
	foreign <- NULL
	iter <- 1
}
if( fun %in% c("casort","rmca")){
	presort <- FALSE
	iter <- 1	
}
if( fun == "preclass" ){
	presort <- FALSE	
}
if( fun %in% c("hamming","hamm") ){
	fun <- "mvclass"
	vs <- as.integer(1)	
}
if( (fun %in% c("class","mvclass") ) ){
	foreign <- ".Call"
	fun <- "mvclass"
	if(joint){
		foreign <- NULL
		fun <- "jointclass"
	}
	if(stepwise){	
		foreign <- NULL
		fun <- "stepclass"
	}
}
if( nd == 2 & fun %in% c("class","mvclass") & vs < 1 ){
	fun <- "quicktile"
	foreign <- ".Call"
	neg <- -1
}else{
	neg <- 1	
}

	
args = c(args,as.integer(vs))


stopifnot( (stepwise & joint) == FALSE )


if(!any(fi)){
	data <- subtable(x, 1:nd, keep.zero = FALSE, allfactor = TRUE)	
	fi <- nd+1
}else{
	data <- subtable(x, 1:nd[-fi], keep.zero = FALSE, freqvar = names(data)[fi], allfactor = TRUE)	
	fi <- nd+1
}
rm(x)
gc()

dims <- sapply(data[,-fi], nlevels )

if( is.logical(perm.cat) ){
	if( all(perm.cat == TRUE) ){
		perm.cat <- rep(1,nd)
	}
	perm.cat <- as.integer( perm.cat )
}else{
	stopifnot(is.numeric(perm.cat) & length(perm.cat) == nd )
	pv <- rep(0,nd)
	pv[perm.cat] <- 1
	perm.cat <- as.integer( pv )
}

	


	
# indicator matrix (code from package nnet):
	imat <- function(s)
	{
		s <- as.factor(s)
		n <- length(s)
		s <- as.factor(s)
		x <- matrix(0, n, length(levels( s )) )
		x[(1:n) + n*(unclass(s)-1)] <- 1
		dimnames(x) <- list(names(s), levels(s))
		return(x)
	}	
	

if( stepwise ){
	
	# a presort step is included in the steptile function
	#if( fun == "class" ){
	#	fun <- "quicktile"
	#	foreign <- ".Call"	
	#}
	
	# preparing for steptile( args )	
	args <- c(list(data, dims, perm.cat, fun, foreign, return.data, return.type, presort, iter), args)
	

	ret <- steptile( args )
	
	return( ret )	
}

# presort (class) for joint and mv

# TODO: add custom presort functions
	orders <- vector(mode="list",length=nd)
	st= Sys.time()
	
	if( presort & !stepwise ){
		Z = xtabs(Freq~.,data = data)
		storage.mode(Z) = "integer"
		dims <- as.integer( dim(Z) )
		cumdims <- c(0,cumsum(dims))
		
		res0 = .Call("preclass",Z,dims,perm.cat)
		
		for( i in 1:nd ){
			orders[[i]] <- res0[ (1+cumdims[i]):cumdims[i+1]]+1
			data[,i] <- factor(data[,i], levels = levels(data[,i])[ orders[[i]] ])	
		}
		# Z = xtabs(Freq~.,data = data)
		#print(Z)
		cat("presort time = ",Sys.time()-st)
		st = Sys.time()
	}else{
		orders <- lapply(dims, function(z) 1:z )
	}


if ( joint ){
	
	# indicator matrix and Burt matrix
	Z <- as.data.frame(do.call(cbind,sapply(data[,1:nd],imat, simplify = FALSE)))
	Z <- t(Z * data[,fi]) %*% as.matrix(Z)
	dims <- as.integer( sapply(data[,1:nd], nlevels) )
	storage.mode(Z) = "integer"	

	
	# preparing mv optimization call
	args <- c(list(as.matrix(Z),dims,perm.cat),args)
		
}else{
	# multidimensional table
	Z <- xtabs(Freq~.,data = data)
	dims <- as.integer( dim(Z) )
	storage.mode(Z) = "integer"	
		
	# preparing mv optimization call
	args <- c(list(Z,dims,perm.cat),args)
}
	
	if(is.null(foreign)){
		optcall <- call( fun )
	}else{
		optcall <- call( foreign,fun )
	}
	for(i in seq_along(args)){
		optcall[[1+i+!is.null(foreign)]] = args[[i]]	
	}
	
	cumdims <- c(0,cumsum(dims))
	
	#data <- as.data.frame(data)
	#data <- data.table(data)
	 
	try(
		res <- eval(optcall)
	)
	
	
	orders2 <- lapply(dims, function(s) 1:s)
	preorders <- list()
	
	if( iter > 1 ){
	for( h in 2:iter ){	
		data0 <- data
		rnd <- lapply(dims, function(s) sample(1:s))
		for(i in 1:nd){
			data0[,i] <- factor(data0[,i], levels <- levels(data0[,i])[rnd[[i]]] )	
		}
		if( presort & !stepwise ){
			Z = xtabs(Freq~.,data = data0)
			storage.mode(Z) = "integer"
			dims <- as.integer( dim(Z) )
			cumdims <- c(0,cumsum(dims))
		
			res0 = .Call("preclass",Z,dims,perm.cat)
		
			for( i in 1:nd ){
				preorders[[i]] <- res0[ (1+cumdims[i]):cumdims[i+1]]+1
				data0[,i] <- factor(data0[,i], levels = levels(data0[,i])[ preorders[[i]] ])	
			}
		}else{
			preorders <- lapply(dims, function(z) 1:z )
		}
		if ( joint ){
		# indicator matrix and Burt matrix
			Z <- as.data.frame(do.call(cbind,sapply(data0[,1:nd],imat, simplify = FALSE)))
			Z <- t(Z * data0[,fi]) %*% as.matrix(Z)
			dims <- as.integer( sapply(data0[,1:nd], nlevels) )
			storage.mode(Z) = "integer"	
		# preparing mv optimization call
			args[[1]] <- as.matrix(Z)
		}else{
		# multidimensional table
			Z <- xtabs(Freq~.,data = data0)
			dims <- as.integer( dim(Z) )
			storage.mode(Z) = "integer"	
		
		# preparing mv optimization call
			args[[1]] <- Z
		}
		
		optcall[[1+1+!is.null(foreign)]] <- args[[1]]	
		
		try(
			res0 <- eval(optcall)
		)
		
		if(neg*res0[length(res0)] > neg*res[length(res)]){
			# res is the result for data[rnd][preorders]
			res <- res0
			orders2 <- mapply(function(a,b) a[b], a = rnd, b = preorders, SIMPLIFY = FALSE)
			
		}
	}
		
	}
		for( i in 1:nd ){
			orders2[[i]] <-  orders2[[i]][res[ (1+cumdims[i]):cumdims[i+1]]+1] 
			data[,i] <- factor(data[,i], levels = levels(data[,i])[ orders2[[i]] ])
			orders[[i]] = 	orders[[i]][ orders2[[i]] ]
		}
	
	#cat("finish time = ",Sys.time()-st)
	crit <- res[length(res)]
	#cat("crit = ",crit)
			
	if( return.data ){
		if(nd > 2 & return.type == "matrix"){
			return.type <- "array"
		}
		return(
			switch(return.type, 
				data.frame 	=  as.data.frame(data)			,
				matrix 		=  as.matrix( xtabs(Freq~.,data=data)) 	,
				table 		=  xtabs(Freq~.,data=data) 			 	,
				ftable 		=  ftable( xtabs(Freq~.,data=data)) 	 	,
				array 		=  array(data, dim = dims) 		 	 	)
		)
	}else{
		return( list(crit, orders) )
		
	}
	
} 



optile.matrix = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "matrix", vs = 0, ...){
	
	dx <- dim(x)
	
	x <- as.data.frame(as.table(x))
	
	NextMethod("optile",object = x, fun = fun, presort = presort, foreign = foreign, args = args, perm.cat = perm.cat, method = method, iter = iter, past = past, scale = scale, 
    freqvar = freqvar, return.data = return.data, return.type = return.type, vs = vs )
} 

optile.array = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "array", vs = 0, ...){
	
	dx <- dim(x)
	
	x <- as.data.frame(x)
	
	NextMethod("optile",object = x, fun = fun, presort = presort, foreign = foreign, args = args, perm.cat = perm.cat, method = method, iter = iter, past = past, scale = scale, 
    freqvar = freqvar, return.data = return.data, return.type = return.type, vs = vs )
} 

optile.table = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "table", vs = 0, ...){
	print("spectabular")
	dx <- dim(x)
	
	x <- as.data.frame(x)
	
	NextMethod("optile",object = x, fun = fun, presort = presort, foreign = foreign, args = args, perm.cat = perm.cat, method = method, iter = iter, past = past, scale = scale, 
    freqvar = freqvar, return.data = return.data, return.type = return.type, vs = vs )
} 

optile.ftable = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "ftable", vs = 0, ...){
	
	dx <- dim(x)
	
	x <- as.data.frame(x)
	
	NextMethod("optile",object = x, fun = fun, presort = presort, foreign = foreign, args = args, perm.cat = perm.cat, method = method, iter = iter, past = past, scale = scale, 
    freqvar = freqvar, return.data = return.data, return.type = return.type, vs = vs )
}    

# ----------------------------------------------------------------------------------------------- #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> stepwise algorithm <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# ----------------------------------------------------------------------------------------------- #

steptile = function( args ){
	data <- args[[1]] 
	dims <- as.integer(args[[2]])
	cumdims <- c(0,cumsum(dims))
	perm.cat <- as.integer(args[[3]])
	return.data <- args[[6]]
	return.type <- args[[7]]
	presort <- args[[8]] 
	iter <- args[[9]]
	fun <- args[[4]]
	foreign <- args[[5]]
	args <- args[-c(4:9)]
		
	nd <- length(dims)	
		#print(" CHECK steptile parameters ")
		#print(args)
	if( fun %in% c("quicktile","stepclass") ){
		fun <- "stepclass"
		foreign <- NULL
		neg <- -1	
	}else{
		neg <- 1	
	}
		
	if(is.null(foreign)){
		optcall <- call( fun )
	}else{
		optcall <- call( foreign,fun )
	}
	# after presort data is sorted by orders
	# orders2 is for other orders from iterations and is based on data[orders]
	orders <- lapply(dims, function(s) 1:s )
	orders2 = lapply(dims, function(s) 1:s)
	
	# initial pair
		tab <- xtabs(data$Freq~data[,1]+data[,2])
		storage.mode(tab) <- "integer"
		args[[1]] <- tab
			
		args[[2]] <- as.integer(dims[1:2])
		args[[3]] <- as.integer(perm.cat[1:2])
		
		if( iter > 1 ){
			data0 <- data	
		}
		if( presort ){
			res0 = .Call("preclass", args[[1]], args[[2]], args[[3]])
		
			for( i in 1:2 ){
				orders[[i]] <- res0[ (1+cumdims[i]):cumdims[i+1] ]+1
				data[,i] <- factor(data[,i], levels = levels(data[,i])[ orders[[i]] ])	
			}
			tab <- xtabs(data$Freq~data[,1]+data[,2])
			storage.mode(tab) <- "integer"
			args[[1]] <- tab

		}
		
		optcall2 <- optcall
		
		for(i in seq_along(args)){
			optcall2[[i+1+!is.null(foreign)]] <- args[[i]]	
		}
		try( 
			res <- eval(optcall2)
		)
		
		if(iter > 1){
			# repeat the procedure for iter random initial orders
			# based on data = data[orders]
			preorders <- lapply(dims, function(z) 1:z )
			for( h in 2:iter ){
				data0 <- data
				#random category orders
				rnd1 <- sample(1:dims[1])
				rnd2 <- sample(1:dims[2])
				
				data0[,1] <- factor(data0[,1], levels <- levels(data0[,1])[rnd1] )	
				data0[,2] <- factor(data0[,2], levels <- levels(data0[,2])[rnd2] )	
				
				tab0 <- xtabs(data0$Freq~data0[,1]+data0[,2])
				storage.mode(tab0) <- "integer"
				args[[1]] <- tab0
			
				if( presort ){
					res0 = .Call("preclass", args[[1]], args[[2]], args[[3]])
		
					for( i in 1:2 ){
						preorders[[i]] <- res0[ (1+cumdims[i]):cumdims[i+1] ]+1
						data0[,i] <- factor(data0[,i], levels = levels(data0[,i])[ preorders[[i]] ])	
					}
					tab0 <- xtabs(data0$Freq~data0[,1]+data0[,2])
					storage.mode(tab0) <- "integer"
					args[[1]] <- tab0
				}
			
				optcall2 <- optcall
				for(i in seq_along(args)){
					optcall2[[i+1+!is.null(foreign)]] <- args[[i]]	
				}
				try( 
					res0 <- eval(optcall2)
				)
				if(neg*res0[length(res0)] > neg*res[length(res0)]){
					# res is the result for data[rnd][preorders] with data = data[orders]
					res <- res0
					orders2[[1]] <- rnd1[ preorders[[1]] ]
					orders2[[2]] <- rnd2[ preorders[[2]] ]
				}
			}
			
		}
			
			for( i in 1:2 ){
				orders2[[i]] <-  orders2[[i]][res[ (1+cumdims[i]):cumdims[i+1]]+1]  
				data[,i] <- factor(data[,i], levels = levels(data[,i])[ orders2[[i]] ])
				orders[[i]] = 	orders[[i]][ orders2[[i]] ]
			}
				
	if( nd > 2 ){	
	# add other variables stepwise
	perm.cat[1:2] <- 0
	f = "data$Freq~data[,1]+data[,2]"
	
	for( s in seq_along(dims[-c(1,2)]) ){
		k <- s+2
		
		f = paste(f,"+ data[,",k, "]",sep="")
		#cat("formula = ",f)
		tab <- xtabs(as.formula(f))
		storage.mode(tab) <- "integer"
		args[[1]] <- tab
		args[[2]] <- as.integer(dims[1:k])
		args[[3]] <- as.integer(perm.cat[1:k])
		
		if( presort ){
			res0 = .Call("preclass", args[[1]], args[[2]], args[[3]])
		
			orders[[k]] <- res0[ (1+cumdims[k]):cumdims[k+1] ]+1
			data[,k] <- factor(data[,k], levels = levels(data[,k])[ orders[[k]] ])	
			
			args[[1]] <- xtabs(as.formula(f))
			storage.mode(args[[1]]) <- "integer"
		}
		
		
		optcall2 <- optcall
		for(i in seq_along(args)){
			optcall2[[1+i+!is.null(foreign)]] <- args[[i]]	
		}
		#print(optcall2)
		try( 
			res <- eval(optcall2)
		)
		
		if( iter > 1 ){
			# repeat the procedure for iter random initial orders
			# data is in orders[[]]
			preorders <- list()
			f0 = "data0$Freq~data0[,1]+data0[,2]"
	
			for( h in 2:iter ){
				data0 <- data
				rnd <- sample(1:dims[[k]])
				data0[,k] <- factor(data0[,k], levels <- levels(data0[,k])[rnd] )	
				
				f = paste(f0,"+ data0[,",k, "]",sep="")
				tab <- xtabs(as.formula(f))
				storage.mode(tab) <- "integer"
				args[[1]] <- tab
						
				if( presort ){
					res0 = .Call("preclass", args[[1]], args[[2]], args[[3]])
		
					preorders[[k]] <- res0[ (1+cumdims[k]):cumdims[k+1] ]+1
					data0[,k] <- factor(data0[,k], levels = levels(data0[,k])[ preorders[[k]] ])	
			
					args[[1]] <- xtabs(as.formula(f))
					storage.mode(args[[1]]) <- "integer"
				}
		
				optcall2 <- optcall
				for(i in seq_along(args)){
					optcall2[[1+i+!is.null(foreign)]] <- args[[i]]	
				}
				try( 
					res0 <- eval(optcall2)
				)
				if(neg*res0[length(res0)] > neg*res[length(res0)]){
					# res is the result for data[rnd][preorders]
					res <- res0
					orders2[[k]] <- rnd[preorders]
				}
			}
		}
			orders2[[k]] <-  orders2[[k]][res[ (1+cumdims[k]):cumdims[k+1]]+1]
			data[,k] <- factor(data[,k], levels = levels(data[,k])[ orders2[[k]] ])
			orders[[k]] <- 	orders[[k]][ orders2[[k]] ]
		
	}
	}
	#print("finish orders")
	#print(orders)
	if( return.data ){
		if(nd > 2 & return.type == "matrix"){
			return.type <- "array"
		}
		return(
			switch(return.type, 
				data.frame 	=  as.data.frame(data)			,
				matrix 		=  as.matrix( xtabs(Freq~.,data=data)) 	,
				table 		=  xtabs(Freq~.,data=data) 			 	,
				ftable 		=  ftable( xtabs(Freq~.,data=data)) 	 	,
				array 		=  array(data, dim = dims) 		 	 	)
		)
	}else{
		return( list( orders) )
		
	}
}




# ----------------------------------------------------------------------------------------------- #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  ca algorithms  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# ----------------------------------------------------------------------------------------------- #


casort <- function(data, dims, perm.cat, ...){
	
	if(! "ca" %in% installed.packages() ){
		install.packages("ca")	
	}
	require(ca)
	
	nd <- length(dims)
		
	ca1 <- ca(data)
	cd2 <- cumsum(dims)
	cd1 <- c(1,cd2[-nd]+1)
	orders <- list()
	
	for( k in 1:nd ){
		orders[[k]] <- 	order( ca1$colcoord[cd1[k]:cd2[k],1])-1
	}
	crit <- ca1$sv[1]
	return(c( unlist(orders),crit ))
}



# ----------------------------------------------------------------------------------------------- #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> iterative joint <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# ----------------------------------------------------------------------------------------------- #


jointclass <- function( Z, dims, perm.cat, ... ){
	
	Z0 <- Z
	cumdims <- c(0,cumsum(dims))
	
			
	nd <- length(dims)
	opt <- FALSE
	cumdims <- c(0,cumsum(dims))
	
	# the current matrix acc to the predefined orders
	Zind <- 1:ncol(Z)
	
	
	
	currcrit <- Inf
	crit <- Inf
	orders <- lapply(dims, function(s) 1:s)
	
	
	while( !opt ){
		opt <- TRUE
		for( k in 1:nd ){
		if( perm.cat[k] > 0 ){
			iv <- (1+cumdims[k]):cumdims[k+1]
			ZX <- Z[ iv, -iv ]
			
			#print(ZX)
			storage.mode(ZX) <- "integer"
			crit0 <- .Call("simplecrit",ZX, as.integer( c(dim(ZX),1) ), as.integer(1) )
			try( 
				res <- .Call("quicktile", ZX, as.integer(dim(ZX)), as.integer(c(1,0)))
			)
			neworder <- res[ 1:dims[k] ] + 1 
			Zind[ iv ] <- Zind[ iv[neworder] ]
			Z <- Z0[Zind, Zind]
			#print("Z")
			#print(Z)
			#orders2[[k]] <-  orders2[[k]][ res[ 1:dims[k] ] + 1 ]
			
			# sort data and original order anew
			orders[[k]] <- 	orders[[k]][ neworder ]
			# crt1 <- res[length(res)]
			crt1 <- .Call("simplecrit",ZX[neworder,], as.integer( c(dim(ZX),1) ), as.integer(1) )
			if( crt1 < crit0 ){
				#print( res[length(res)] / crit0)
				#print( crt1 / crit0)
				#cat("in dim ",k)
				opt <- FALSE	
			}
			#else no dimension has changed => optimum reached
		}}
	}

	return(c( unlist(orders)-1, 0) )
	
}


# ----------------------------------------------------------------------------------------------- #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> iterative steps <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# ----------------------------------------------------------------------------------------------- #


stepclass <- function( tab, dims, perm.cat, ... ){
	nd <- length(dims)
	ft <- ftable(tab)
	class(ft) <- "matrix"
	tt2 <- t(ft)
	m <- ncol(tt2)
	n <- nrow(tt2)
	
	dims2 <- dim(tt2)
	zid <- which(colSums(tt2)==0)
		if(any(zid)){
			tt2 <- tt2[,-zid]
			dims2[2] <- ncol(tt2)
		}
		storage.mode(tt2) = "integer"
		# perm.cat[nd-1] is 0 in all but the initial step
		pc <- as.integer(c(perm.cat[nd],perm.cat[nd-1]))
		try( 
			res <- .Call("quicktile", tt2, as.integer(dims2), pc)
		)
		if(nd > 2){
			res <- c( unlist(lapply(dims[-nd],function(s) 1:s) )-1, res[1:n],res[length(res)])
		}else{
			res <- res[c((n+1):(n+m),1:n,n+m+1)]
		}
	return( res )
}


# ----------------------------------------------------------------------------------------------- #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  ca algorithms  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# ----------------------------------------------------------------------------------------------- #


csvd <- function(data, dims, perm.cat, ...){

	nd <- length(dims)
	di <- c(1:nd)
	
	coords <- list()
	cumdims <- c(0,cumsum(dims))
	
	#presort with ca algorithm
	
	# indicator matrix (code from package nnet):
	imat <- function(s)
	{
		s <- as.factor(s)
		n <- length(s)
		s <- as.factor(s)
		x <- matrix(0, n, length(levels( s )) )
		x[(1:n) + n*(unclass(s)-1)] <- 1
		dimnames(x) <- list(names(s), levels(s))
		return(x)
	}
	# back from table to
	tt <- data
	
	#print("tt after")
	#print(tt)
	
	N <- sum(tt)	
	data <- as.data.frame(data)
	
	Z <- as.data.frame(do.call(cbind,sapply(data[,1:nd],imat, simplify = FALSE)))
	Z <- t(Z * data[,nd+1]) %*% as.matrix(Z)
	storage.mode(Z) = "integer"	
	
	res0 <- casort(Z,dims,perm.cat)

	
	orders = list()
	for( i in 1:nd ){
		orders[[i]] <-  res0[ (1+cumdims[i]):cumdims[i+1]]+1 
		data[,i] <- factor(data[,i], levels = levels(data[,i])[ orders[[i]] ])
	}
	
	# table in new order
	tt <- xtabs(Freq~.,data=data)
	
	
	#data <- as.data.frame(data)
	#for( i in 1:nd ){
	#	data[,i] <- as.factor(data[,i])	
	#}
	
	
	opt <- FALSE
	while(!opt){
		opt <- TRUE
	for( i in 1:nd ){
		m <- dim(tt)[i]
		if(nd < 3){
			profs <- apply(tt,i, function(z){
				z <- cumsum(z)
				z/sum(z)	
			})
		}else{
			profs <- apply(tt,i, function(z){
					for(s in 1:(nd-1)){
						z <- apply(z,-s,cumsum)
					}
			z/max(z)#sum	
			})
		}
		#ex.profs = rowMeans(profs)
		#M = (profs -	ex.profs)/sqrt(ex.profs)
		M <- profs
		M[is.nan(M)] <- 0

		
		S <- svd(M)
		coords[[i]] <- S$v[,1]
		#print(dim(S$u))
		#print(dim(S$v))
		#ord = order(S$v[,1])
		#if(ord[1] > ord[length(ord)]){
		#	coords[[i]] = rev(coords[[i]])	
		#}
		#cat("rows ",order(S$u[,1]))
		#cat("cols ",order(S$v[,1]))
		
		neworder <- order(coords[[i]])
		data[,i] = factor(data[,i],levels = levels(data[,i])[neworder])
		
		if( !all(orders[[i]] == orders[[i]][neworder]) ){
			opt <- FALSE	
		}
		orders[[i]] <- orders[[i]][neworder]
		
		tt <- xtabs(Freq~.,data=data)
		#print(tt)
		
	}	
	}
	
	# handle reverse orders
	tt2 <- xtabs(data$Freq~data[,1]+data[,2])
	if( (tt2[1,dims[2]] + tt2[dims[1],1]) > (tt2[1,1] + tt2[dims[1],dims[2]]) ){
		orders[[1]] <- rev(orders[[1]])	
	}
	
	return(c(unlist(orders)-1,0))
}


rmca <- function( data, dims, perm.cat, ... ){

	data <- as.data.frame(data)
	A <- xtabs(Freq~.,data = data)
	nd <- length(dims)

if(nd > 2){
	N <- sum(A)
	coords <- list()
ind <- 1:nd

#if(.Platform$OS.type == "unix" && "multicore" %in% .packages(all.available = TRUE) ){
#	require("multicore")
#	coords <- mclapply(ind, function(i){
		
#	mz <- apply(A, -i, sum)
#	mz <- mz/N

#	M <- matrix( apply(A, i, function(s){
#		(s/sum(s) - mz)
#		}), ncol= dims[i])

#	S <- svd(M[which(rowSums(M) > 0),])
#	return(S$v[,1])
#})	

#}else{
	coords <- lapply(ind, function(i){
		
	mz <- apply(A, -i, sum)
	mz <- mz/N

	M <- matrix( apply(A, i, function(s){
		(s/sum(s) - mz)
		}), ncol= dims[i])

	S <- svd(M[which(rowSums(M) > 0),])
	return(S$v[,1])
	})
#}

		
		
#for( i in 1:nd ){	
#	mz <- apply(A, -i, sum)
#	mz <- mz/N
#
#	M <- matrix( apply(A, i, function(s){
#		(s/sum(s) - mz)
#		}), ncol= dims[i])
#
#	S <- svd(M[which(rowSums(M) > 0),])
#	coords[[i]] <- S$v[,1]
#}

	for( i in 1:nd ){
		data[,i] = factor(data[,i],levels = levels(data[,i])[order(coords[[i]])])
	}	
	return(c(unlist(lapply(coords,order))-1,0))
}else{
	# TODO: better solution here
	return(extracat:::casort(data,dims,perm.cat))	
}
}
