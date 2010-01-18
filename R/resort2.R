resort <-
function(V = iset(which=iset.cur()),sel = iset.selected(), listen = FALSE){

nv = names(V)
i1 = which(sapply( strsplit(nv,"\\."), function(x) "C" %in% x))

n = length(i1)
i0 = 1:n
i3 = which(sapply( strsplit(nv,"\\."), function(x) "S" %in% x))
i4 = which(sapply( strsplit(nv,"\\."), function(x) "I" %in% x))
V0 = V[,i0]
V3 = V[,i3]
V4 = V[,i4]

indK  = which( sapply(V0,class) != "numeric" )
nK = length(indK)
TT = lapply(V0[,indK],function(x){
tapply(1:nrow(V0),x,I)
})
tmp = mapply( function(w,y){
unlist(sapply( w, function(x){
sel.ind = which(x %in% sel)
not.ind = which( !(x %in% sel) )
n = length(x)
nsel = length(sel.ind)
t0 = vector(length=n,mode="numeric")
t0[sel.ind] =  sort(y[x],decreasing=T)[nsel:1]
t0[not.ind] =  sort(y[x],decreasing=T)[n:(nsel+1)]
return(t0)
}))
},TT,data.frame(V3[,indK]),SIMPLIFY =F)
ind = lapply(TT,function(x) unlist(x))
V3[,indK] = mapply(function(x,y) x[order(y,decreasing=F)],tmp,ind)


if(.GlobalEnv$sort.individual){
			e1 = new.env()
				e1$M3 = V3
				e1$M2 = V3+V4
				for( z in 2:nK ){
					i = indK[z]
						sapply(TT[[z]], function(s){
							sel.ind = which(s %in% sel)
							not.ind = which( !(s %in% sel) )
							e1$M3[s[sel.ind],i] = e1$M3[s[sel.ind],i][rank(e1$M2[s[sel.ind],i-1])]
							e1$M2[s[sel.ind],i] = e1$M2[s[sel.ind],i][rank(e1$M2[s[sel.ind],i-1])]
							e1$M3[s[not.ind],i] = e1$M3[s[not.ind],i][rank(e1$M2[s[not.ind],i-1])]
							e1$M2[s[not.ind],i] = e1$M2[s[not.ind],i][rank(e1$M2[s[not.ind],i-1])]
						})
				}
			V3 = e1$M3
}
Vnew = V3+V4

for( i in 1:n ){
	ivar.update(V[[ i1[i] ]],Vnew[,i],batch = T)
}
iset.updateVars()
if(listen){
listen()
}

}

