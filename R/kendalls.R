kendalls <- function(x){
	cs <- colSums(x)
	rs <- rowSums(x)
	css <- rev(cumsum(rev(cs)))
	rss <- rev(cumsum(rev(rs)))
	n <- nrow(x)
	m <- ncol(x)
	N <- sum(x)
	xx <- sum(cs*(css-cs)[1:m])
	y <- sum(rs*(rss-rs)[1:n])
	
	
	storage.mode(x) <- "integer"
	x2 <- x[,m:1]
	
	dims <- as.integer(c(n,m))
	crt1 <- .Call("classcrit",x,dims,as.integer(0))
	crt2 <- .Call("classcrit",x2,dims,as.integer(0))
	
	tau <- (crt1-crt2)/sqrt(max(xx,1))/sqrt(max(1,y))
#scrt <- crt1/x/y*N^2
#return(c(tau,crt1,scrt,n,m,N))
	return(tau)
}


