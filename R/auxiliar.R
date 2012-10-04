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
		
		
		
Burt = function(x){
	if(inherits(x,"table")){
        return(Burt.table(x))
    }
	if(!"Freq" %in% names(x)){
	   x <- subtable(x, 1:ncol(x))
	}
	   fi <- which(names(x) == "Freq")
	   nd <- ncol(x)-1
	Z <- as.data.frame(do.call(cbind, sapply(x[, 
											 1:nd], imat, simplify = FALSE)))
	Z <- t(Z * x[, fi]) %*% as.matrix(Z)
return(Z)
}


idat = function(x, allcat = FALSE){
	if("Freq" %in% names(x)){
		fi <- which(names(x) =="Freq")
		s <- x[,fi]
		x <- x[,-fi]	
	}else{
		s <- NULL	
	}
	if(allcat){
		ret<-sapply(x,imat)
	}else{
		ret<-sapply(x,function(z){
			y<-imat(z)
			y <- y[,-ncol(y),drop=FALSE]
		})
	}
	ret <- as.data.frame(do.call(cbind,ret))
	if(!is.null(s)){
		ret$Freq <- s	
	}
	return(ret)
}
Burt.table <- function(x){
	stopifnot(inherits(x,"table"))
	nd <- length(dim(x))
B2 <- NULL
for(i in 1:nd){
	B1 <- NULL
	for(j in 1:nd){
		B0<-apply(x,unique(c(i,j)),sum)
        if(i == j){
            B0 <- diag(B0)
        }
        colnames(B0) <- dimnames(x)[[j]]
        rownames(B0) <- dimnames(x)[[i]]
		 B1 <- cbind(B1,B0)
		}
		B2 <- rbind(B2,B1)
} 
return(B2)
}