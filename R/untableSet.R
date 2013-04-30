untableSet2 <-
function(data, freqvar = NULL){
	
	data = as.data.frame(data)
	if("Freq" %in% names(data) & is.null(freqvar)){ freqvar = "Freq" }
	stopifnot(freqvar %in% names(data))
	
	ind = which(names(data) != freqvar)
	fi = which(names(data) == freqvar)
	names(data)[fi] = "Freq"
	
	n = ncol(data)
	m = nrow(data)
	X = data.frame(matrix(ncol=n-1,nrow=0))
	zind = which(data$Freq > 1)
	zero = which(data$Freq == 0)

	X = sapply(zind,function(x) spread(data[x,ind],nrow = data[x,fi]))
	X = do.call("rbind",X)
	vn = names(data)[ind]
	names(X)=vn
	X = rbind(X,as.matrix(data[c(1:m)[-c(zero,zind)],ind]))
	return(suppressWarnings(data.frame(X)))
}
untableSet <- function(data, freqvar = NULL){
	data <- as.data.frame(data)
	if("Freq" %in% names(data) & is.null(freqvar)){ freqvar <- "Freq" }
	
	stopifnot(freqvar %in% names(data))
	ind = which(names(data) != freqvar)
	fi = which(names(data) == freqvar)
	names(data)[fi] <- "Freq"
	n <- ncol(data)
	m <- nrow(data)
	lvls <- lapply(data[,-fi],function(z) levels(as.factor(z)))
	
	lol <- apply(data,1,function(z){
		 matrix(rep(z[-fi],z[fi]),nrow = as.integer(z[fi]),byrow=T)	
	})
	X <- as.data.frame(do.call(rbind,lol))
	names(X) <- names(data)[-fi]
	for(i in ind){
		X[,i] <- factor(X[,i], levels = lvls[[i]]	)
	}
	return(X)
}
