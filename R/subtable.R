
subtable <- function(data, cols, freqvar = NULL, keep.zero=FALSE, allfactor =FALSE, return.type = class(data)){
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
			#ss = count2(data, cols, weights = data[,fid])
			ss = count(data, cols, wt_var = fid)
		}else{
			#ss = count2(data, cols)
			ss = count(data, cols)
		}
		names(ss)[names(ss)=="freq"] <- "Freq"
		
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
	if("table" %in% return.type )	ss <- xtabs(Freq~.,data=ss)
		
	return(ss)
}
	
# count2 = function (df, vars = NULL, weights = NULL) 
# {
    # if (is.vector(df)) {
        # df <- data.frame(x = df)
    # }
    # if (!is.null(vars)) {
        # vars <- as.quoted(vars)
        # df <- quickdf(eval.quoted(vars, df))
    # }
# if(!is.null(weights)){
    # id <- plyr:::ninteraction(df, drop = TRUE)
    # u_id <- !duplicated(id)
    # labels <- df[u_id, , drop = FALSE]
    # labels <- labels[order(id[u_id]), , drop = FALSE]
    # Freq <- xtabs(weights~id)
	# class(Freq) = "vector"
    # unrowname(data.frame(labels, Freq))
# }else{
	# id <- plyr:::ninteraction(df, drop = TRUE)
    # u_id <- !duplicated(id)
    # labels <- df[u_id, , drop = FALSE]
    # labels <- labels[order(id[u_id]), , drop = FALSE]
    # Freq <- tabulate(id, attr(id, "n"))
    # unrowname(data.frame(labels, Freq))
# }
# }

subtable.table = function(x,cols){
	x <- as.table(x)
	x2 <- apply(x,cols,sum)
	dim(x2) <- dim(x)[cols]
	return(as.table(x2))
	
}



