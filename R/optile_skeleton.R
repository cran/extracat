
optile = function(x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
freqvar = NULL, return.data = TRUE, return.type = "data.frame", vs = 0, tree = NULL, ...){
	UseMethod("optile")	
}
#optile = function(x,...){
#UseMethod("optile")	
#}

optile.default = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "data.frame", vs = 0, tree = NULL, ...){



fi <- which(names(x) %in% c("Freq",freqvar))
nd <- ncol(x)-any(fi)


if(!is.null(tree)){
	stopifnot(length(tree) == nd | all(tree == "hc"))
	if(all(tree == "hc")){
		tree <- rep("hc", nd)	
	}
	if(!(fun %in% c("treeclass","class"))){
		cat("trees defined -> using fun = treeclass.")	
	}
	fun <- "treeclass"
	foreign <- NULL
	method <- NULL
	args <- c( 0, args ) # WARUM?
	args[[1]] <- tree
	names(args)[[1]] <- "treelist"
	presort <- FALSE
}
	
stepwise <- ifelse(is.null(method), FALSE, method %in% c("stepwise","sw","step"))
joint <- ifelse(is.null(method), FALSE, method %in% c("joint","pairwise","pw"))
# mv <- method %in% c("mv","multivar","multivariate")



# workarounds:
if(nd == 2){
	if( fun == "rmca" ){
		fun <- "casort"	
	}
	
}

if( fun %in% c("ca","casort") ){
	fun <- "casort"
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
	if(nd == 2){
		fun = "quickhamm2d"
		foreign=".Call"
	}else{
		simpleWarning("hamming for more than 2 dimensions is under construction. code is still erroneous. Please use the unweighted classification criterion.")	
		fun <- "mvclass"
		#vs <- as.integer(1)
	}
	
		
}
if( (fun %in% c("class","mvclass") ) ){
	foreign <- ".Call"
	fun <- "quickmv"
	if(joint){
		foreign <- NULL
		fun <- "jointclass"
	}
	if(stepwise){	
		foreign <- NULL
		fun <- "stepclass"
	}
}
if( fun %in% c("quickmv","quickhamm") ){
	args <- c( args, as.integer(0), as.integer(0), as.integer(0)) 
}
neg <- 1
if( nd == 2 & fun %in% c("class","quickmv") & vs < 1 ){
	fun <- "quicktile"
	foreign <- ".Call"
	neg <- -1
}
if( fun == "quickhamm2d" ){
	neg <- -1	
}

# vs is the 4th standard parameter after data, dims and perm.cat	
args = c(as.integer(vs), args)


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
	st <- Sys.time()
	
	if( presort & !stepwise ){
		Z <- xtabs(Freq~.,data = data)
		storage.mode(Z) <- "integer"
		dims <- as.integer( dim(Z) )
		cumdims <- c(0,cumsum(dims))
		
		res0 <- .Call("preclass",Z,dims,perm.cat)
		
		for( i in 1:nd ){
			orders[[i]] <- res0[ (1+cumdims[i]):cumdims[i+1]]+1
			data[,i] <- factor(data[,i], levels = levels(data[,i])[ orders[[i]] ])	
		}
		# Z = xtabs(Freq~.,data = data)
		#print(Z)
#cat("presort time = ",Sys.time()-st)
		st <- Sys.time()
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
#print("start")
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
	
	scaled.crit <- NULL	
		
	if( fun == "mvclass" ){
		if( vs == 1 & nd == 2){
			scaled.crit	 <- crit / ihcrit(Z)
		}
		if( vs == 1 ){
			scaled.crit <- crit/ iccrit(Z)	
		}
	}
	if( fun == "quicktile" ){
		scaled.crit <- crit/ iccrit(Z)	
	}
		
		
		
	if( return.data ){
		if(nd > 2 & return.type == "matrix"){
			return.type <- "array"
		}
		ret <- switch(return.type, 
				data.frame 	=  as.data.frame(data)			,
				matrix 		=  as.matrix( xtabs(Freq~.,data=data)) 	,
				table 		=  xtabs(Freq~.,data=data) 			 	,
				ftable 		=  ftable( xtabs(Freq~.,data=data)) 	 	,
				array 		=  array(xtabs(Freq~.,data=data), dim = dims) )
				if(return.type == "matrix" & nd < 3){
					class(ret) <- "matrix"	
				}
		if(!is.null(tree)){
			tree <- attr(res, "tree")
#			for(k in seq_along(tree)){
#				if("order" %in% names(tree[[k]])){
#						tree[[k]]$order <- tree[[k]]$order[orders[[k]]]
#				}	
#			}
			attr(ret, "tree") <- tree	
		}
		attr(ret,"criterion") <- crit
		if(!is.null(scaled.crit)){
			attr(ret,"scaled.criterion") <- scaled.crit
		}
		attr(ret, "orders") <- orders
		return(ret)
	}else{
		return( list(crit, scaled.crit, orders) )
		
	}
	
} 



optile.matrix = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "matrix", vs = 0, tree = NULL, ...){
	
	dx <- dim(x)
	
	x <- as.data.frame(as.table(x))
	
	NextMethod("optile",object = x, fun = fun, presort = presort, foreign = foreign, args = args, perm.cat = perm.cat, method = method, iter = iter, past = past, scale = scale, 
    freqvar = freqvar, return.data = return.data, return.type = return.type, vs = vs , tree = tree)
} 

optile.array = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "array", vs = 0, tree = NULL, ...){
	
	dx <- dim(x)
	
	x <- as.data.frame(as.table(x))
	
	NextMethod("optile",object = x, fun = fun, presort = presort, foreign = foreign, args = args, perm.cat = perm.cat, method = method, iter = iter, past = past, scale = scale, 
    freqvar = freqvar, return.data = return.data, return.type = return.type, vs = vs, tree = tree )
} 

optile.table = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "table", vs = 0, tree = NULL, ...){
#print("spectabular")
	dx <- dim(x)
	
	x <- as.data.frame(x)
	
	NextMethod("optile",object = x, fun = fun, presort = presort, foreign = foreign, args = args, perm.cat = perm.cat, method = method, iter = iter, past = past, scale = scale, 
    freqvar = freqvar, return.data = return.data, return.type = return.type, vs = vs , tree = tree)
} 

optile.ftable = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "ftable", vs = 0, tree = NULL, ...){
	
	dx <- dim(x)
	
	x <- as.data.frame(x)
	
	NextMethod("optile",object = x, fun = fun, presort = presort, foreign = foreign, args = args, perm.cat = perm.cat, method = method, iter = iter, past = past, scale = scale, 
    freqvar = freqvar, return.data = return.data, return.type = return.type, vs = vs, tree = tree )
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
	
	if(! "ca" %in% .packages(all.available = TRUE) ){
		install.packages("ca")	
	}
	require(ca)
	
	nd <- length(dims)
		
	ca1 <- ca(data)
	cd2 <- cumsum(dims)
	cd1 <- c(1,cd2[-nd]+1)
	orders <- list()
	
	for( k in 1:nd ){
		#orders[[k]] <- 	order( ca1$colcoord[cd1[k]:cd2[k],1])-1
		orders[[k]] <- do.call(order, as.data.frame( round( ca1$colcoord[cd1[k]:cd2[k],], 12) ))-1 
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
		##print(dim(S$u))
		##print(dim(S$v))
		##ord = order(S$v[,1])
		##if(ord[1] > ord[length(ord)]){
		##	coords[[i]] = rev(coords[[i]])	
		##}
		##cat("rows ",order(S$u[,1]))
		##cat("cols ",order(S$v[,1]))
		
		#neworder <- order(coords[[i]])
		neworder <- do.call(order, as.data.frame(round(S$v, 12)))
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


# different ways of including (cluster) trees:
# >> standard input, fun = "treeclass", tree <- list(tree1, tree2, ... ) 
# 	--> tree1, tree2, ... is 
# 		- a tree object (contains merge and order)
# 		- a function which generates a tree for the categories based on the mv table
# 		- "hc" with optional arguments method, link and dist
# 		- 0, NULL, FALSE for no tree
# 	--> the return object has an additional attribute "tree" containing the generated or reordered trees
# 	--> the number of clusters is set to the number of factor levels
 
# >> a list of cluster results and/or variables
# 		- a simple factor
# 		- a cluster object with a tree (merge and order)
# 	--> the return object is a data.frame (or as specified by return.type) 
# 		with an additional attribute "treelist" containing the generated or reordered trees
# 	--> for each cluster tree the number of clusters has to be specified:
# 		> in args( ncl <- c(1,2,3,4,...) )
# 		> in tree1$ncl, tree2$ncl, ...

treeclass <- function( data, dims, perm.cat, vs, treelist = "hc", hcargs = list(), ... ){

	nd <- length(dims)
	#print(data)
	treelist <- as.list(treelist)
	if(length(treelist) == 1){
		treelist <- as.list( rep(treelist[[1]], nd) )	
	}
	treelengths <- rep(0,nd)
	treevec <- list()
	orders <- list()
	
	tt <- data
	data <- as.data.frame(data)
	
	#print("treelist")
	#print(treelist)
	
	for(i in 1:nd){
		tr <- treelist[[i]]
		istree <- all( c("merge","order") %in% names(tr) )# height?
		
		if(perm.cat[i] > 0){
			
		if( istree ){
			# "hclust" indicates a tree from hclust, amap or flashClust
			# "twins" is a tree from the diana or agnes function in the cluster package
			if( dims[i] == length(tr$order) ){
				stree <- tr
				data[,i] <- factor(data[,i], levels = levels(data[,i])[stree$order])
			}else{
				stree <- subtree(tr, k = dims[i])
				treelist[[i]] <- stree
				# CHECK ORDER CHECK !!!!
			}
			treevec <- c(treevec,untree(stree$merge,stree$order))
			orders[[i]] <- stree$order
			#data[,i] = factor(data[,i],levels = levels(data[,i])[ tr$order ])
			treelengths[i] <- nrow(stree$merge)
			#data[,i] <- factor(data[,i], levels = levels(data[,i])[stree$order])
		}else{
			# no tree is defined unless there is both obj$merge and obj$order which we use
			# TODO: obj$dist for plotting? => better in separate function
			if( toString(tr) == "hc" ){
					# TODO
						mth <- ifelse(is.null(hcargs$method), "complete", hcargs$method)
						dist.mth <- ifelse(is.null(hcargs$dist.method), "euclidean", hcargs$dist.method)
						M <- matrix( apply(tt, i, function(s){
							(s/sum(s))
							}), ncol= dims[i])
						
						if( "amap" %in% .packages(all.available = TRUE) ){
							require(amap)
							stree <- hcluster(t(M), method = dist.mth, link = mth)
						}else{
							DM <- dist(t(M))
#if( "flashClust" %in% .packages(all.available = TRUE) ){
#require(flashClust)
#	stree <- flashClust(DM, method = mth)
#}else{
								stree <- hclust(DM, method = mth)
#}
						}
						treevec <- c(treevec,untree(stree$merge,stree$order))
						orders[[i]] <- stree$order
						treelist[[i]] <- stree
						#data[,i] = factor(data[,i],levels = levels(data[,i])[ tr$order ])
						treelengths[i] <- nrow(stree$merge)
						data[,i] <- factor(data[,i], levels = levels(data[,i])[stree$order])
			}else{
			#if( toString(tr) %in% c("0","none","FALSE","NULL","NA","") )
					# TODO
					treelengths[i] <- 0
			}
		}
		
		}else{
			treelengths[i] <- 0
		}
	
	}
	# mv table with possibly updaten category orders
	tt <- xtabs(Freq~.,data=data)
	# trees etc. are ready
	# CHECK: do we have to rearrange the levels acc. to tr$order ??
	#tt <- xtabs(Freq~.,data=data)
storage.mode(tt) <- "integer"
storage.mode(treevec) <- "integer"
storage.mode(treelengths) <- "integer"
	minc <- as.integer(0)
#print("----- CHECK PARAMETERS IN R -----")
#cat("treevec = ", treevec)
#print("----------")
#cat("treelengths = ", treelengths)
#print("----------")
#print(tt)
#print("----------")
	try(
		res <- .Call("quickmv", tt, as.integer(dims), as.integer(perm.cat), as.integer(vs), minc, treevec, treelengths)
	)
	cumdims <- c(0,cumsum(dims))
	#print(res)
	for(k in 1:nd){	
		if(treelengths[k] > 0){
			treelist[[k]]$order <- treelist[[k]]$order[res[ (1+cumdims[k]):cumdims[k+1] ] + 1]
			res[ (1+cumdims[k]):cumdims[k+1] ] <- orders[[k]][ res[ (1+cumdims[k]):cumdims[k+1] ] + 1 ] - 1
		}
	}
	attr(res, "tree") <- treelist
	return(res)
}


untree <- function(x,y, ind = FALSE){
	stopifnot(ncol(x) == 2)
	ng <- nrow(x)
	ind1 <- vector(mode="list",length=ng)
	ind2 <- vector(mode="list",length=ng)
	for(i in 1:ng){
		if(x[i,1] < 0){
			ind1[[i]] <- c( ind1[[i]] , -x[i,1] )
		}else{
			ind1[[i]] <- c( ind1[[i]] , ind1[[ x[i,1] ]], ind2[[ x[i,1] ]])
		}
		if(x[i,2] < 0){
			ind2[[i]] <- c( ind2[[i]] , -x[i,2] )
		}else{
			ind2[[i]] <- c( ind2[[i]] , ind1[[ x[i,2] ]], ind2[[ x[i,2] ]] )
		}
	}	
	ind1 <- lapply(ind1, function(z) match(z,y))
	ind2 <- lapply(ind2, function(z) match(z,y))
	if(!ind){
		ret <- as.vector( rbind( sapply( ind1, function(z) z[1]), sapply(ind2, range )) )
	}else{
		ret <- mapply(function(x,y) list( x,y ), x = ind1,y=ind2, SIMPLIFY =FALSE)	
	}
	return(ret)
}

optile.list = function (x, fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
    freqvar = NULL, return.data = TRUE, return.type = "table", vs = 0, tree = NULL, k = NULL, h = NULL, ...){
	
	#stopifnot( all(labels(ncl) %in% c("h","k", "")) )
	#if(!is.null(k) & !is.null(h)){
	#	warning("k is used before h")	
	#}
	
	if(is.null(h)){
			h <- rep(NA,length(x))
		}
	if(is.null(k)){
			k <- rep(NA,length(x))
		}
		if(length(k) == 1){
			k = rep(k, length(x))	
		}
	if(length(h) == 1){
			h = rep(h, length(x))	
		}
	ncl <- mapply( function(k,h)(
		if(is.na(k)){
			s <- h
			names(s) <- "h"
			return(s)	
		}else{
			s <- k
			names(s) <- "k"
			return(s)
		}
	), k = k, h = h ,SIMPLIFY = FALSE)
	
 
		 		mp <-	mapply( function(z,kh){
					if( is.vector(z) ){
						#factor(z)
						z
					}else{
						lab <- names(kh)
						if(is.na(kh)){
							if( !is.null(z$k) ){
								kh <- z$k
								lab <- "k"
							}
							if( !is.null(z$h) ){
								kh <- z$h
								lab <- "h"
							}	
						}
						stopifnot( names(kh) %in% c("k","h") )
						stopifnot( "merge" %in% names(z) )
						stopifnot( "order" %in% names(z) )
			
						if( lab == "k" ){
							tree <- subtree(z, k = kh, h = NULL)
						}else{
							stopifnot( "height" %in% names(z) )
							tree <- subtree(z, k = NULL, h = kh)	
						}
						ret <- factor(tree$data, levels = tree$order)
						attr(ret,"tree") <- tree
						return(ret)
					}
				}, z = x, kh = ncl, SIMPLIFY = FALSE) 
	#print(mp)						
	treelist <- lapply(mp, function(z) attr(z,"tree"))
	data <- do.call(cbind,mp)
# how to get the trees out?
	#print(sapply(data,levels))
	x <- subtable(data, 1:ncol(data), allfactor=TRUE)
	#print(x)
	#args[[ length(args)+1 ]] <- treelist
	#print(treelist)
	#names(args)[[ length(args) ]] <- "treelist"
	
	NextMethod("optile",object = x, fun = "treeclass", presort = FALSE, foreign = NULL, args = args, perm.cat = perm.cat, method = NULL, iter = iter, past = past, scale = scale, 
    freqvar = freqvar, return.data = return.data, return.type = return.type, vs = vs, tree = treelist )
}


subtree = function(tree, k = NULL, h = NULL){
	
	data <- tree$order 
	itr <- untree(tree$merge, tree$order, ind = TRUE)
	
	if( is.null(k) & is.null(h) ){
		stop("Neither k nor h specified")	
	}
	if( !is.null(k) & !is.null(h) ){
		stop("Unclear specification: found both k and h.")	
	}
	if(!is.null(h)){
		k <- sum(tree$height > h)	
	}
	
	nit <- length(itr)
	
	itr <- itr[ (nit-k+2):nit ]
	mrg <- tree$merge[(nit-k+2):nit,,drop=FALSE]
	
	#compute level order
	
	tmrg <- t(mrg)
	mins <- sapply(itr,function(x) sapply(x, min))
	#mins[tmrg < (nit-k+2)] <- rank(mins[tmrg < (nit-k+2)])
	#tmrg[tmrg < (nit-k+2)] <- -c(1:k)
	#ord <- rank(mins[tmrg < (nit-k+2)])
	tmrg[tmrg < (nit-k+2)] <- -rank(mins[tmrg < (nit-k+2)])
	tmrg[tmrg >= (nit-k+2)] <- tmrg[tmrg >= (nit-k+2)] - nit + k - 1
	mrg <- t(tmrg)
	rm(tmrg)
	
	# TODO: better solution than a loop?
	#s <- 1
	for( i in 1:(k-1) ){
		i1 <- mrg[i,1]
		if( i1 < 0 ){
			data[ itr[[i]][[1]] ] <- -i1 #-s
			#s <- s+1
		}
		i2 <- mrg[i,2]
		if( i2 < 0 ){
			data[ itr[[i]][[2]] ] <- -i2 #-s
			#s <- s+1
		}
	}
	# CHECK
	labs <- as.vector(t(mrg))
	labs <- labs[labs < 0]
	labs <- match(1:k, -labs)
	# CHECK
	
	ret <- new("list")
	ret$merge <- mrg
	ret$order <- 1:k
	ret$labels <- labs
	#ret$data <- factor(data[ order(tree$order) ], levels <- ord)
	ret$data <- factor(data[ order(tree$order) ])
	ret$height <- tree$height[ (nit-k+2):nit ]
	class(ret) <- "hclust"
	return(ret)
}




noderange <- function(mrg,ord, xval = NULL, ind = FALSE){
	stopifnot(ncol(mrg) == 2)
	ng <- nrow(mrg)
	
	ind1 <- lapply(1:ng,function(z) c(ng+2,0))
	ind2 <- ind1
	for(i in 1:ng){
		if(mrg[i,1] < 0){
			ind1[[i]][2] <- max( match(abs(mrg[i,1]),ord), ind1[[i]][2] )
			ind1[[i]][1] <- min( match(abs(mrg[i,1]),ord), ind1[[i]][1] )
		}else{
			ind1[[i]][1] <- min( ind1[[i]][1] , ind1[[ mrg[i,1] ]][1], ind2[[ mrg[i,1] ]][1])
			ind1[[i]][2] <- max( ind1[[i]][2] , ind1[[ mrg[i,1] ]][2], ind2[[ mrg[i,1] ]][2])
		}
		if(mrg[i,2] < 0){
			ind2[[i]][2] <- max( match(abs(mrg[i,2]),ord), ind2[[i]][2] )
			ind2[[i]][1] <- min( match(abs(mrg[i,2]),ord), ind2[[i]][1] )
		}else{
			ind2[[i]][1] <- min( ind2[[i]][1] , ind1[[ mrg[i,2] ]][1], ind2[[ mrg[i,2] ]][1])
			ind2[[i]][2] <- max( ind2[[i]][2] , ind1[[ mrg[i,2] ]][2], ind2[[ mrg[i,2] ]][2])
		}
	}	
	
	ret <- mapply(function(x,y) list( x,y ), x = ind1,y=ind2, SIMPLIFY =FALSE)	
	return(ret)
}

nodexc <- function(mrg,ord, xval = NULL, ind = FALSE){
	stopifnot(ncol(mrg) == 2)
	ng <- nrow(mrg)

	cnt <- as.list(rep(0,ng))
	
	for(i in 1:ng){
		if(mrg[i,1] < 0){
			a <-  (match(abs(mrg[i,1]),ord)-0.5)/(ng+1)
		}else{
			a <- cnt[[mrg[i,1]]][3]	
		}
		if(mrg[i,2] < 0){
			b <-  (match(abs(mrg[i,2]),ord)-0.5)/(ng+1)
		}else{
			b <- cnt[[mrg[i,2]]][3]	
		}
		cnt[[i]] <- c(a,b,(a+b)/2)
	}	

	return(cnt)
}



nodeheight <- function(mrg,ord,ht){
	stopifnot(ncol(mrg) == 2)
	ng <- nrow(mrg)
	#ind1 <- vector(mode="list",length=ng)
	#ind2 <- vector(mode="list",length=ng)
	LH <- as.list(ht)
	RH <- LH
	for(i in 1:ng){
		if(mrg[i,1] > 0){
			LH[[i]] <- LH[[i]] - ht[[mrg[i,1]]]
		}
		if(mrg[i,2] > 0){
			RH[[i]] <- RH[[i]] - ht[[mrg[i,2]]]
		}
	}	

	ret <- mapply(function(x,y) list( x,y ), x = LH,y=RH, SIMPLIFY =FALSE)	
	return(ret)
}

drawt = function(tree, vertical = FALSE, vp = NULL, ... ){
	if(is.null(tree)){
		return(invisible(TRUE))	
	}
	ng <- length(tree$order)
	maxht <- max(tree$height)
	nht <- nodeheight( tree$merge, tree$order, tree$height/maxht	)
	#nrg <- noderange( tree$merge, tree$order )
	#nleft <- lapply(nrg, function(z) mean(z[[1]])/ng )
	#nright <- lapply(nrg, function(z) mean(z[[2]])/ng )
	nxc <- nodexc(tree$merge, tree$order)
	if( !vertical ){
	mapply( function(x,hlr, h){
		grid.lines( x = c(x[1],x[2]), y = c(h,h), vp = vp )
		grid.lines( x = c(x[1],x[1]), y = c(h-hlr[[1]],h), vp = vp )
		grid.lines( x = c(x[2],x[2]), y = c(h-hlr[[2]],h) , vp = vp)
		}, x = nxc , hlr = nht, h = as.list(tree$height/maxht))
	}else{
		mapply( function(x,hlr, h){
		grid.lines( y = 1-c(x[1],x[2]), x = 1-c(h,h) , vp = vp)
		grid.lines( y = 1-c(x[1],x[1]), x = 1-c(h-hlr[[1]],h) , vp = vp)
		grid.lines( y = 1-c(x[2],x[2]), x = 1-c(h-hlr[[2]],h) , vp = vp)
		}, x = nxc , hlr = nht, h = as.list(tree$height/maxht))
	}
invisible(TRUE)
	
}


tfluctile = function(x, tree = NULL, dims = c(1,2), tw = 0.2, border = NULL, tile.col= hsv(0.1,0.1,0.1,alpha=0.6),
 bg.col = NULL, vp = NULL, ... ){
	stopifnot( !is.null( attr(x,"tree") ) | !is.null(tree) )
	
	
	if(is.null(border)){
		border <- 0.05	
	}
	if( !is.null(tree) ){
		tree1 <- tree[[ dims[1] ]]	
		tree2 <- tree[[ dims[2] ]]	
	}else{
		tree1 <- attr(x,"tree")[[ dims[1] ]]	
		tree2 <- attr(x,"tree")[[ dims[2] ]]
	}
	if(length(dim(x)) > 2){
		x <- as.data.frame(as.table(x))
		x <- xtabs(x$Freq~x[,dims[1]]+x[,dims[2]])
	}
	lefts <- ifelse( all(is.na(tree1)), 0, tw )
	tops <- ifelse( all(is.na(tree2)), 0, tw )
	
	if(is.null(vp)){
		grid.newpage()	
	}else{
		pushViewport(vp)	
	}
	
	vp00 <- viewport(x=lefts,y=1-tops, just=c("left","top"), width = 1-lefts, height = 1-tops)
	
	if(!all(is.na(tree1))){
	vp0L <- viewport(x=0,y=1-tops, just=c("left","top"), width = lefts, height = (1-tops))
	pushViewport(vp0L)
	pushViewport( viewport(0.5,0.5,0.96,1-2*border) )
	drawt(tree1,vertical=TRUE, vp = NULL)
	upViewport(2)
	}
	if(!all(is.na(tree2))){
	vp0T <- viewport(x=lefts,y=1, just=c("left","top"), width = (1-lefts), height = tops)	
	pushViewport(vp0T)
	pushViewport( viewport(0.5,0.5,1-2*border,0.96) )
	drawt(tree2,vertical=FALSE, vp = NULL)
	upViewport(2)
	}
	pushViewport(vp00)
	fluctile(x, add = TRUE, tile.col = tile.col, bg.col = bg.col)
	upViewport()
	return(invisible(TRUE))
}
