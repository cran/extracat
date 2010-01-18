subtable <-
function(dset,cols,keep.zero=TRUE,allfactor=TRUE, freqvar = NULL){
	dset = data.frame(dset)
	if( "Freq" %in% names(dset) & is.null(freqvar) ){ freqvar = "Freq"}
	
	n = length(cols)
	if(allfactor){
		int = which(sapply(dset[,cols],is.integer))
		for( i in int ){
			dset[,cols[i]] = as.factor(dset[,cols[i]])
		}
	}
	indK = which( sapply(dset[,cols],function(x) inherits(x,"factor") ) )
	orig.labs = lapply(dset[,cols],function(x) levels(as.factor(x)))

	if( !is.null(freqvar) ){
		stopifnot( freqvar %in% names(dset) )
		names(dset)[ which(names(dset) == freqvar) ] = "Freq"
		if(keep.zero){
			S = do.call("paste", as.list(c("Freq~",names(dset[cols]),sep="+")) )
			subtable = data.frame(xtabs(formula(S),data=dset))
			return(subtable)
		}
		if(!keep.zero){
			Tt = do.call("paste",c(dset[,cols],sep=":"))
			fTt = xtabs(dset$Freq~Tt)
			ftab = colsplit(c(names(fTt)),split=":",names=names(dset)[cols])
			ftab=data.frame(cbind(ftab,fTt[][][1:nrow(fTt)]))
			names(ftab)[ncol(ftab)]="Freq"
			subtable=ftab
			subtable$Freq = as.integer(subtable$Freq)
			if(allfactor){
				subtable = data.frame(sapply(subtable[,c(1:n)],factor),subtable$Freq)
			}
			for( i in indK ){
				subtable[,i] = factor(subtable[,i],levels = orig.labs[[i]])	
			}
			names(subtable)[ncol(subtable)]="Freq"
			
			return(subtable)
		}
	}else{
		Tt = do.call("paste",c(dset[,cols],sep=":"))
		TT = data.frame( cbind( Tt,rep( 1,length(Tt) ) ) )
		names(TT) = c("Tt","Freq")
		TT$Freq=as.numeric(TT$Freq)
		fTt = xtabs(TT$Freq~Tt)
		ftab = colsplit(c(names(fTt)),split=":",names=names(dset)[cols])
		ftab=data.frame(cbind(ftab,fTt[][][1:nrow(fTt)]))
		names(ftab)[ncol(ftab)]="Freq"
			for( i in indK ){
				ftab[,i] = factor(ftab[,i],levels = orig.labs[[i]])	
			}
		if(keep.zero){
			S = do.call("paste", as.list(c("Freq~",names(ftab)[1:n],sep="+")) )
			subtable=data.frame(xtabs(formula(S),data=ftab))
			if(allfactor){
				subtable = data.frame(sapply(subtable[,c(1:n)],factor),as.numeric(subtable$Freq))
				names(subtable)[n+1] = "Freq"	
			}
			for( i in indK ){
				subtable[,i] = factor(subtable[,i],levels = orig.labs[[i]])	
			}
			return(subtable)
		}else{
			if(allfactor){
				ftab = data.frame(sapply(ftab[,1:n],factor),as.numeric(ftab$Freq))
				names(ftab)[n+1] = "Freq"	
			}
			for( i in indK ){
				ftab[,i] = factor(ftab[,i],levels = orig.labs[[i]])	
			}
			return(ftab)
		}
	}
}

