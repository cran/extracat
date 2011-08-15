getsets = function(x){
	s = x[,1]
	s2 = x[,2]
	ind = seq_along(s)
	ind2 = seq_along(s2)
	while(length(ind) > 0){
		
		s[ind] = sapply(s[ind],function(z){
			c(z[z<0], x[ z[z>0], ])
		},simplify=FALSE)
		ind = which(sapply(s, function(v){
				any(v > 0)
			}))
		}
	while(length(ind2) > 0){	
		s2[ind2] = sapply(s2[ind2],function(z){
			c(z[z<0], x[ z[z>0], ])
		},simplify=FALSE)
		ind2 = which(sapply(s2, function(v){
				any(v > 0)
			}))
		}
		return(list(s,s2))
}


optilehc = function(M, d.method = "euclidean", h.method = "centroid", 
	critfun = "class", perm.cat = c(1,2), symmetric = FALSE, return.data = TRUE){



critfunz = c("class","hamming","euclidean","ME")

stopifnot(critfun %in% critfunz)

cm = match(critfun, critfunz)

pv = as.integer(1:2 %in% perm.cat)
sm = as.integer(symmetric)
st0 = Sys.time()


	
	n = nrow(M)
	m = ncol(M)
	
	#prop.table?
	row.mat = M/rowSums(M)
	col.mat = t(M)/colSums(M)
	
	row.dm = Dist(row.mat, method = d.method)
	col.dm = Dist(col.mat, method = d.method)
	hcr = flashClust(row.dm, method = h.method)
	hcc = flashClust(col.dm, method = h.method)
	
M = M[hcr$order,hcc$order]

col.sets = getsets(hcc$merge)
row.sets = getsets(hcr$merge)
c1 = sapply(col.sets[[1]], function(z){
		as.integer(min(match(-z,hcc$order))-1)
	})
c2 = sapply(col.sets[[2]], function(z){
		as.integer(min(match(-z,hcc$order))-1)
	})
c3 = sapply(col.sets[[2]], function(z){
		as.integer(max(match(-z,hcc$order))-1)
	})
r1 = sapply(row.sets[[1]], function(z){
		as.integer(min(match(-z,hcr$order))-1)
	})
r2 = sapply(row.sets[[2]], function(z){
		as.integer(min(match(-z,hcr$order))-1)
	})
r3 = sapply(row.sets[[2]], function(z){
		as.integer(max(match(-z,hcr$order))-1)
	})
nrp = as.integer(length(r1))
ncp = as.integer(length(c1))
dims = as.integer(c(n,m))
storage.mode(M)="integer"

hcopt = .Call("opthclust", M, dims, as.integer(cm), pv,sm, r1, r2, r3, c1, c2, c3, nrp, ncp, as.integer(0:(n-1)), as.integer(0:(m-1))) 

rind = hcopt[1:n]
cind = hcopt[(n+1):(n+m)]
val = hcopt[n+m+1]
M2 = M[rind,cind]


cat("optimization complete after ",Sys.time()-st0)

if(return.data){
	attr(M2,"var.orders") = list(rind,cind)
	attr(M2,"criterion") = val
	return(M2)
}else{
	return(list(rind,cind,val))
}
}
