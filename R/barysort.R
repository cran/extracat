barysort = function(x, vs = 1){
	n <- nrow(x)
	m <- ncol(x)
	
	

	rs <- rowSums(x)
	cs <- colSums(x)
	optimal <- FALSE
	crt <- classcrit(x)
	while(!optimal){
		
		b0 <- ( x %*% 1:m )/rs
		ord <- order(b0)
		x <- x[ord,]
		rs <- rs[ord]
		x2 <- x
		crt2 <- classcrit(x2)
		
		b1 <- ( t(1:n) %*% x )/cs
		ord <- order(b1)
		x <- x[,ord]
		cs <- cs[ord]
		
		crt0 <- classcrit(x)
		
		if(crt0 <= crt & crt2 <= crt){
			optimal <- TRUE	
			if(crt2 <= crt){
				crt <- crt2
				x <- x2
			}
		}else{
			crt <- crt0	
		}
		
		
	}
	return(x)
}