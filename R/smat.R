smat = function(n,m,N = 2*n*m,k = 1,mix = TRUE){
	v1 = rep(1:n,m)
	v2 = rep(1:m,each=n)
	
	dmax = max(n,m)-1
	dv =  (1 - abs(v1/n-v2/m) )^k
	
	N = max(1,N-n*m)
	
	Freq = sapply(dv, function(z){
			return(rpois(1,lambda=z*(N-m*n)/m*n+1))
		})
	SM = xtabs(Freq~v1+v2)
	if(mix){
			return(SM[sample(1:n),sample(1:m)])
		}else{
			return(SM)
		}
}


smat2 = function(n,m,N){
	pv = runif(m)
	M=sapply(pv, function(x) table( c(0:n,rbinom(N,n,x)) ))
	
	return(M[sample(1:n),sample(1:m)]-1)
}