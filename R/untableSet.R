untableSet <-
function(dset, freqvar = NULL){
	dset = data.frame(dset)
	if("Freq" %in% names(dset) & is.null(freqvar)){ freqvar = "Freq" }
	dset = data.frame(dset)
	stopifnot(freqvar %in% names(dset))
	ind = which(names(dset) != freqvar)
	fi = which(names(dset) == freqvar)
	names(dset)[fi] = "Freq"
	n = ncol(dset)
	m = nrow(dset)
	X = data.frame(matrix(ncol=n-1,nrow=0))
	zind = which(dset$Freq > 1)
	zero = which(dset$Freq == 0)

	X = sapply(zind,function(x) spread(dset[x,ind],nrow = dset[x,fi]))
	X = do.call("rbind",X)
	vn = names(dset)[ind]
	names(X)=vn
	X = rbind(X,as.matrix(dset[which( !(c(1:m) %in% c(zero,zind)) ),ind]))
	return(suppressWarnings(data.frame(X)))
}

