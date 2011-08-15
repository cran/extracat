#subtable <-
#function(dset,cols,keep.zero=FALSE,allfactor=TRUE, freqvar = NULL){
#	dset = data.frame(dset)
#	if( "Freq" %in% names(dset) & is.null(freqvar) ){ freqvar = "Freq"}
#	
#	n = length(cols)
#	if(allfactor){
#		int = which(sapply(dset[,cols],is.integer))
#		for( i in int ){
#			dset[,cols[i]] = as.factor(dset[,cols[i]])
#		}
#	}
#	indK = which( sapply(dset[,cols],function(x) inherits(x,"factor") ) )
#	orig.labs = lapply(dset[,cols],function(x) levels(as.factor(x)))
#
#	if( !is.null(freqvar) ){
#		stopifnot( freqvar %in% names(dset) )
#		names(dset)[ which(names(dset) == freqvar) ] = "Freq"
#		if(keep.zero){
#			S = do.call("paste", as.list(c("Freq~",names(dset[cols]),sep="+")) )
#			subtable = data.frame(xtabs(formula(S),data=dset))
#			return(subtable)
#		}
#		if(!keep.zero){
#			Tt = do.call("paste",c(dset[,cols],sep="::"))
#			fTt = xtabs(dset$Freq~Tt)
#			ftab = colsplit(c(names(fTt)),pattern="::",names=names(dset)[cols])
#			ftab=data.frame(cbind(ftab,fTt[][][1:nrow(fTt)]))
#			names(ftab)[ncol(ftab)]="Freq"
#			subtable=ftab
#			subtable$Freq = as.integer(subtable$Freq)
#			if(allfactor){
#				subtable = data.frame(sapply(subtable[,c(1:n)],factor),subtable$Freq)
#			}
#			for( i in indK ){
#				subtable[,i] = factor(subtable[,i],levels = orig.labs[[i]])	
#			}
#			names(subtable)[ncol(subtable)]="Freq"
#			
#			return(subtable)
#		}
#	}else{
#		Tt = do.call("paste",c(dset[,cols],sep="::"))
#		TT = data.frame( cbind( Tt,rep( 1,length(Tt) ) ) )
#		names(TT) = c("Tt","Freq")
#		TT$Freq=as.numeric(TT$Freq)
#		fTt = xtabs(TT$Freq~Tt)
#		ftab = colsplit(c(names(fTt)),pattern="::",names=names(dset)[cols])
#		ftab=data.frame(cbind(ftab,fTt[][][1:nrow(fTt)]))
#		names(ftab)[ncol(ftab)]="Freq"
#			for( i in indK ){
#				ftab[,i] = factor(ftab[,i],levels = orig.labs[[i]])	
#			}
#		if(keep.zero){
#			S = do.call("paste", as.list(c("Freq~",names(ftab)[1:n],sep="+")) )
#			subtable=data.frame(xtabs(formula(S),data=ftab))
#			if(allfactor){
#				subtable = data.frame(sapply(subtable[,c(1:n)],factor),as.numeric(subtable$Freq))
#				names(subtable)[n+1] = "Freq"	
#			}
#			for( i in indK ){
#				subtable[,i] = factor(subtable[,i],levels = orig.labs[[i]])	
#			}
#			return(subtable)
#		}else{
#			if(allfactor){
#				ftab = data.frame(sapply(ftab[,1:n],factor),as.numeric(ftab$Freq))
#				names(ftab)[n+1] = "Freq"	
#			}
#			for( i in indK ){
#				ftab[,i] = factor(ftab[,i],levels = orig.labs[[i]])	
#			}
#			return(ftab)
#		}
#	}
#}
subtable <- function(data, cols, freqvar = NULL, keep.zero=FALSE, allfactor =FALSE){
	data=as.data.frame(data)
	names(data)[which(names(data)==freqvar)] = "Freq"
	
	if(allfactor){
		int = which(sapply(data[,cols], function(v) is.integer(v)|is.numeric(v)))
		for( i in int ){
			data[,cols[i]] = as.factor(data[,cols[i]])
		}
	}

	if(!keep.zero){
		if("Freq" %in% names(data)){
			fid = which(names(data)=="Freq")
			ss = count2(data, cols, weights = data[,fid])
		}else{
			ss = count2(data, cols)
		}
		ss <- ss[ss$Freq > 0,]
	}else{
		if("Freq" %in% names(data)){
			fid = which(names(data)=="Freq")
			ss = as.data.frame(xtabs(Freq~.,data=data[,c(cols,fid)]))
		}else{
			ss = as.data.frame(table(data[,c(cols)]))
		}
	}
#	if(allfactor){
#			ss = data.frame(sapply(ss[,-ncol(ss)],factor),ss$Freq)
#			names(ss)[ncol(ss)] = "Freq"
#	}
	return(ss)
}
	
count2 = function (df, vars = NULL, weights = NULL) 
{
    if (is.vector(df)) {
        df <- data.frame(x = df)
    }
    if (!is.null(vars)) {
        vars <- as.quoted(vars)
        df <- quickdf(eval.quoted(vars, df))
    }
if(!is.null(weights)){
    id <- plyr:::ninteraction(df, drop = TRUE)
    u_id <- !duplicated(id)
    labels <- df[u_id, , drop = FALSE]
    labels <- labels[order(id[u_id]), , drop = FALSE]
    Freq <- xtabs(weights~id)
	class(Freq) = "vector"
    unrowname(data.frame(labels, Freq))
}else{
	id <- plyr:::ninteraction(df, drop = TRUE)
    u_id <- !duplicated(id)
    labels <- df[u_id, , drop = FALSE]
    labels <- labels[order(id[u_id]), , drop = FALSE]
    Freq <- tabulate(id, attr(id, "n"))
    unrowname(data.frame(labels, Freq))
}
}




