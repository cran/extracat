cpcp = function(V,ord = NULL, freqvar = NULL,numerics = NULL, gap.type = "equal.tot",na.rule = "omit",
spread=0.3,gap.space = 0.2
,sort.individual=FALSE, jitter = FALSE,
plot = TRUE, return.df = !plot, ... ){
	ieXt = FALSE
	scr.res=c(1280,1024)
	if( plot & .Platform$OS.type == "unix" ){
		if(!("JGR" %in% .packages(all.available=FALSE))){
			cat("Please use the 'JGR' console for iplots. To download it visit\n 'http://www.rosuda.org/software/'")
			return(invisible(TRUE))
			}
	}
			
	if(plot & !( "iplots" %in% .packages(all.available = TRUE) ) ){
		cat("Installing required package 'iplots'...")
		install.packages('iplots')   
		}
	if(plot){ 	
			require(iplots)
	}
	
	
	V = data.frame(V)
	m = nrow(V)
	
	if(!is.null(numerics)){
		stopifnot(is.null(freqvar))
		for(j in numerics){
			if(is.integer(V[,j]) & jitter){
				V[,j] = as.numeric(V[,j]) + runif(m,-0.5,0.5) # +- 1?
			}else{
				V[,j] = as.numeric(V[,j])
			}
		}
	}

	stopifnot(gap.type %in% c("equal.gaps","equal.tot","spread"))
	stopifnot(spread <= 1)
	stopifnot(gap.space < 1)
	stopifnot(dim(V) > 1)

	if(na.rule == "omit"){
		V = na.omit(V)	
	}else{ # how handle NAs in numeric variables ??
		V = sapply(V,function(x){
				if(!class(x) == "numeric"){
				w = x
				w[is.na(x)] = "N/A"
					return(w)
				}else{
					return(x)	
				}
			})
	}
	
	if( "Freq" %in% names(V) & is.null(freqvar) ){ freqvar = "Freq" }
		
	is.ft = !is.null(freqvar)
	

	if(is.ft){
		stopifnot(freqvar %in% names(V))
		if(is.null(ord)){ ord = which(names(V) != freqvar)	}
		n = length(ord)
		
		fi = which(names(V) == freqvar)
		#not.ord = which( !( 1:ncol(V) %in% c(ord,fi) ) )
		not.ord = c(1:ncol(V))[-c(ord,fi)]
		
		VS = subtable(V,c(ord,not.ord),keep.zero=F,allfactor=F,freqvar=freqvar)
		lvls = lapply(VS[,1:n],function(x){
				 if(class(x) != "numeric") levels(as.factor(x)) 
				})
		VS = untableSet(VS, freqvar = freqvar)
	}else{
		if(is.null(ord)){ ord = 1:ncol(V)	}
		
		n = length(ord)
		not.ord = c(1:ncol(V))[-ord]
		
		VS = V[,c(ord,not.ord)]
		lvls = lapply(VS[,1:n],function(x){
				 if(!inherits(x, "numeric")) levels(as.factor(x)) 
				})
	}
	

	V = VS		
	m = nrow(V)
	classes = sapply(V[,1:n],class)
	indK = which(classes != "numeric")
	indS = which(classes == "numeric")
	nK = length(indK)
	for(j in indK){
		V[,j] = factor(V[,j],levels = lvls[[j]])
	}

	m = nrow(V)

	V2 = matrix(ncol=n,nrow=m)
	V3 = matrix(0,ncol=n,nrow=m)
	V4 = matrix(0,ncol=n,nrow=m)
	
	S = vector(mode="list",length=nK)
	S = sapply(V[,indK],function(x) table(x),simplify=FALSE)

	V4 = V[,1:n]
	V4[,indK] = sapply(V4[,indK], as.integer)
	
	ord = do.call(order,c(V4,decreasing=F))
	V = V[ord,]
	V4 = V4[ord,]	
			if( gap.type == "equal.gaps" ){
				nmax = max(sapply(V[,indK],nlevels))
			#	gap.space = min(gap.space,0.9/(maxn-1))
			}
			#------------------------------------------------------------------------------#
			
			curmax = sapply(S,max)
			TT = sapply(V[,indK],function(x){
					tapply(1:nrow(V),x,I)
				},simplify=FALSE)
			
			
			nlvl = sapply(V[,indK],nlevels)
	
			if( gap.type == "equal.tot" ){
					eta = sapply(S,function(x) x/sum(x),simplify=FALSE)
					sc = sapply(nlvl,function(x) c(0:(x-1)) *( 2*gap.space/( (x - 1)*(1 - gap.space))-1 ),simplify=FALSE) # includes the integer value
					eta2 = sapply(eta,function(x){
							nj = length(x)
							return(cumsum(c(0,x[1:(nj-1)])+ c(0,x[2:(nj)])))
						}, simplify=FALSE)
					shift = mapply(function(x,y) x + y,	eta2,  sc, SIMPLIFY = FALSE)
				}
			if( gap.type == "equal.gaps" ){
					eta = sapply(S,function(x) x/sum(x),simplify=FALSE)
					sc = sapply(nlvl,function(x) c(0:(x-1)) *( 2*gap.space/( (nmax - 1)*(1 - gap.space))-1),simplify=FALSE)
					eta2 = sapply(eta,function(x){
							nj = length(x)
							return(cumsum(c(0,x[1:(nj-1)])+ c(0,x[2:(nj)])))
						},simplify=FALSE)
					shift = mapply(function(x,y) x + y,	eta2,  sc,SIMPLIFY=FALSE)
			}
			if(gap.type == "spread"){
					eta = sapply(S,function(x) x/max(x)*spread/2,simplify=FALSE)
					shift = lapply(nlvl,function(x) rep(0,x))
			}
			
			seqs = lapply(1:nK,function(y){

					unlist(
						sapply(c(1:nlvl[y]),function(x){
							return( seq(shift[[y]][x]-eta[[y]][x],shift[[y]][x]+eta[[y]][x],2*c(eta[[y]]/S[[y]])[x])[1:S[[y]][x]] )
						})
					)
				})

			ind = lapply(TT,function(x) unlist(x))

			ord.seqs = mapply(function(x,y) return(x[order(y)]),seqs,ind)


			V3[,indK] = ord.seqs
			V2 = V3+V4
		
			.GlobalEnv$sort.individual = FALSE
			if( sort.individual ){
				e1 = new.env()
				e1$M3 = V3
				e1$M2 = V2
				for( z in 2:nK ){
					i = indK[z]
						sapply(TT[[z]], function(s){
							e1$M3[s,i] = e1$M3[s,i][rank(e1$M2[s,i-1])]
							e1$M2[s,i] = e1$M2[s,i][rank(e1$M2[s,i-1])]
						})
				}
			V3 = e1$M3
			V2 = V3+V4
			.GlobalEnv$sort.individual = TRUE
			}
	
			#.GlobalEnv$sort.individual = FALSE
			#if( sort.individual ){
			#	e1 = new.env()
			#	e1$M3 = V3
			#	e1$M2 = V2
			#	for( z in (nK-1):1 ){
			#		i = indK[z]
			#			sapply(TT[[z]], function(s){
			#				e1$M3[s,i] = e1$M3[s,i][rank(e1$M2[s,i+1])]
			#				e1$M2[s,i] = e1$M2[s,i][rank(e1$M2[s,i+1])]
			#			})
			#	}
			#V3 = e1$M3
			#V2 = V3+V4
			#.GlobalEnv$sort.individual = TRUE
			#}
			
			
			
			vn = names(V)[1:n]
			cpcp=data.frame(V,V2,V3,V4)

			names(cpcp)=c(names(V),paste("C",vn,sep="."),paste("S",vn,sep="."),paste("I",vn,sep="."))	
		if(plot){
			s = iset.new("icpcp",cpcp)
			ipcp(s[,(1:n) + ncol(V) ])
			# setting the size and location as the upper half of the screen
			iplot.location(x=10, y=10, relative=FALSE, plot=iplot.cur())
			iplot.size(width=(scr.res[1])*(n-1)/n, height=scr.res[2]/2.2, plot=iplot.cur())
		
			
		}
		#itext <- function(x, y=NULL, labels=seq(along=x), ax=NULL, ay=NULL, ..., plot=iplot.cur())
			
		if(return.df){
			return(cpcp)	
		}else{
			return(invisible(TRUE))	
		}

}
