

fluctile <- function(tab, gap.prop = 0.1, hsplit = TRUE, env = NULL, abbrev = 3, border = NULL, label = TRUE,...  ){
	
	dm <- dim(tab)
	maxv <- max(tab)
	nd <- length(dm)
	
	data <- as.data.frame(as.table(tab))
	
	
# logical vector for the split directions
	if(length(hsplit) == 1){
		hsp <- rep( c(hsplit,!hsplit),ncol(data) )[1:(ncol(data)-1)] 
		hsplit <- hsp
	}else{
		hsp <- hsplit
	}
	for(i in which(!hsp) ){
		data[,i] <- factor(data[,i], levels = rev(levels(data[,i])))	
	}
	tab <- xtabs(Freq~., data = data)
	
# prepare for border and labs
	nx <- min(sum(hsp*label),2)
	ny <- min(sum(!hsp*label),2)
	if(is.null(border)){
		border <- 	0.05* c( (nx+1), (ny+1) )
	}
	
	wph <- dm[2]/dm[1]
#dev.new( width = 1000*(1-gap.prop)*wph + gap.prop*1000, height = 1000 )
	dev.new()
	
	vp0 = viewport(x = border[1]*nx/(nx+1) + (1-border[1])/2 , y = border[2]/(ny+1) + (1-border[2])/2, width = 1-border[1], height = 1-border[2],name="base")
	if(length(dm)==2){
		
		pushViewport(vp0)
		grid.rect(gp = gpar(fill=rgb(0,0,0,alpha=0.1),col=NA))
		
# handle the last 2 dimensions
		if( hsp[1] ){
			if( hsp[2] ){
				dim(tab) <- c(1, prod(dim(tab)))	
			}else{
				tab <- as.table(t(tab))	
				dm <- dim(tab)
			}
		}else{
			if( !hsp[2] ){
				dim(tab) <- c( prod(dim(tab)), 1)	
			}
		}
		gridfluc(tab,gap.prop, maxv=maxv)
		popViewport()
		
# add labels
		if(label){
			rn <- abbreviate(dimnames(tab)[[1]],abbrev)
			cn <- abbreviate(dimnames(tab)[[2]],abbrev)
			
			m <- dm[2]
			xc <- seq(0.5,m-0.5) / m
			xc <- xc - gap.prop/(m-1)/2
			xc <- xc/(1 - gap.prop/(m-1))
			
			
			n <- dm[1]
			yc <- seq(0.5,n-0.5) / n
			yc <- yc - gap.prop/(n-1)/2
			yc <- yc/(1 - gap.prop/(n-1))
			
			
			vpR <- viewport(x = border[1]/4, y = 0.5, width = border[1]/2, height = 1-border[2],name="rowlabs")
			pushViewport(vpR)
			grid.text(rn, x = 0.5, y = yc, gp = gpar(cex=1.2), rot = 65)
			popViewport()
			
			vpC <- viewport(x = border[1]/2 + (1-border[2])/2, y = 1-border[2]/4, width = 1-border[1], height = border[2]/2,name="collabs")
			pushViewport(vpC)
			grid.text(cn, x = xc, y = 0.5,rot = 25, gp = gpar(cex=1.2))
			
			popViewport()
		}
		return(vp0)
	}else{
		e1 <- new.env()
		ltrs <- expand.grid(letters,letters,letters)
		
		e1$vpn <- paste(ltrs[,1],ltrs[,2],ltrs[,3],sep="")
		e1$k <- 0
		hsplit <- !hsplit
		ss <- fluctree(dm,parent=vp0, hsplit=hsp, gap.prop=gap.prop, env=e1)
		grid.newpage()
		pushViewport(ss)
		seekViewport("base")
		
		e1$k <- 0	
		
#flucplot creates a viewport tree in the e1 environment
		flucplot(tab = tab, gap.prop = gap.prop, hsplit = hsp, env=e1)
		
# use the tree in e1 to plot the 2-dimensional fluctuation diagrams
		mapply(function(x,y,hs,gp) {
#pushViewport(ss)
			   seekViewport(y)
			   
#			if(hsp[1] & hsp[2]) dim(x) <- c(1, prod(dim(x)))
#			if(!hsp[1] & !hsp[2]) dim(x) <- c( prod(dim(x)), 1)
#			if(hsp[1] & !hsp[2]) x <- as.table(t(x))
			   if(hsp[length(hsp)-1] & hsp[length(hsp)]) dim(x) <- c(1, prod(dim(x)))
			   if(!hsp[length(hsp)-1] & !hsp[length(hsp)]) dim(x) <- c( prod(dim(x)), 1)
			   if(hsp[length(hsp)-1] & !hsp[length(hsp)]) x <- as.table(t(x))
			   
#TODO: do something similar for the labels
			   
			   gridfluc(x,gp, maxv=maxv)
			   }, x = e1$tablist, y = e1$namlist, hs = e1$hslist, gp = e1$gplist)	
		
		upViewport()
# go back to surface
		for(i in 1:(nd-2)){
			upViewport()
		}
		
#######################################################################################################	
# -------------------------------------------- LABELING  -------------------------------------------- #
#
		labs <- lapply(data[,-(nd+1)],function(s) abbreviate(levels(as.factor(s)),abbrev))
		
		vp1 <- viewport(x = border[1]*nx/(nx+1)/2, y = border[2]/(ny+1) + (1-border[2])/2, width = border[1]*nx/(nx+1), height = 1-border[2],name="ylab")
		pushViewport(vp1)
#grid.rect(0.5,0.5,1,1,gp=gpar(fill=rgb(0,0,0,alpha=0.1)))
		
		suppressWarnings( label <- label & rep( TRUE, nd ))
		
# create labels for the y-axis
		
		ind <- which(label & !hsp)
		rpt <- c(1,cumprod( dim(tab)[ind] ))
		
		nlvl <- dim(tab)[ind]
		gaps <- 0
		blocksize <- 1
		
		for( i in 1:ny ){
# viewport for the label column
			currvp <- viewport( x = 1/ny/2+ (i-1)/ny, y = 0.5, width = 1/ny, height = 1)
			pushViewport(currvp)
			
			labs <- rep(levels(data[, ind[i]]), rpt[i])
			nl <- length(labs)
			
# the gaps for the new dimension
			newgaps <- c(0:(nlvl[i]-1)) * blocksize * gap.prop / (nlvl[i]-1)
			newgaps <- rep(newgaps, rpt[i])
			
			gaps <- rep(gaps, each = nlvl[i])
			
			x <- 0.5
			if( i == 1 ){
				y <- seq(1/nl/2,1-1/nl/2,1/nl)*((1-gap.prop)^i) + newgaps 
			}else{
# the old coordinates +- the new gaps and cellwidths
				y <- rep(y,each = nlvl[i]) + newgaps - max(newgaps)/2 + (1-gap.prop)*blocksize * rep(seq(1/nlvl[i]/2,1-1/nlvl[i]/2,1/nlvl[i])-0.5, rpt[i])
			}
			
			grid.text(labs, x = x, y = y , rot = 65, gp = gpar(cex=1.5))
#grid.points( x = rep(x,length(y)), y = y , gp = gpar(col="red"))
			popViewport()
			
			blocksize <- blocksize * (1-gap.prop) / nlvl[i]
		}
		popViewport()
		
		vp2 <- viewport(x = border[1]*nx/(nx+1) + (1-border[1])/2, y = 1 - border[2]*ny/(ny+1)/2, width = 1-border[1], height = border[2]*ny/(ny+1),name="xlab")
		pushViewport(vp2)
#grid.rect(0.5,0.5,1,1,gp=gpar(fill=rgb(0,0,0,alpha=0.1)))
		
		ind <- which(label & hsp)
		rpt <- c(1,cumprod( dim(tab)[ind] ))
		
		nlvl <- dim(tab)[ind]
		gaps <- 0
		blocksize <- 1
		
		for( i in 1:nx ){
# viewport for the label column
			
			currvp <- viewport( x = 0.5, y = 1 - 1/nx/2 - (i-1)/nx, width = 1, height = 1/nx)
			pushViewport(currvp)
			
			labs <- rep(levels(data[, ind[i]]), rpt[i])
			nl <- length(labs)
			
			
			newgaps <- c(0:(nlvl[i]-1)) * blocksize * gap.prop / (nlvl[i]-1)
			newgaps <- rep(newgaps, rpt[i]) 
			
#cat("labeling variable", names(data)[ind[i]], " with nlvl = ", nlvl[i], " and levels = ", labs," and rpt = ", rpt[i]) 
			
			y <- 0.5
			if( i == 1 ){
				x <- seq(1/nl/2,1-1/nl/2,1/nl)*((1-gap.prop)^i) + newgaps 
			}else{
				x <- rep(x,each = nlvl[i]) + newgaps - max(newgaps)/2 + (1-gap.prop)*blocksize * rep(seq(1/nlvl[i]/2,1-1/nlvl[i]/2,1/nlvl[i])-0.5, rpt[i])
			}
#cat( "lab.x = ",x)
			grid.text(labs, x = x, y = y , rot = 25, gp = gpar(cex=1.5))
#grid.points( x = x, y = rep(y,length(x)) , gp = gpar(col="red"))
			popViewport()
			
			blocksize <- blocksize * (1-gap.prop) / nlvl[i]
		}
#
# -------------------------------------------- LABELING  -------------------------------------------- #
#######################################################################################################	
		
		return(ss)
	}
}	

flucplot <- function(tab, gap.prop, hsplit, env, ...){
	if(! "tablist" %in% ls(env) ){
		env$tablist = list()
		env$namlist = list()
		env$hslist = list()
		env$gplist = list()
		env$k2 = 0
	}
	
	dm <- dim(tab)
	
#	if(length(hsplit) == 1){
#		hsplit = !hsplit	
#	}else{
		hsplit = hsplit[-1]	#hsplit[-length(hsplit)]#	
#	}
	
#print(env$k)
#print(env$vpn[env$k])
#print(dm)
	
	if(length(dm) == 2){
		
		env$k2 <- env$k2+1
		k2 <- env$k2
		env$tablist[[k2]] <- tab
		env$namlist[[k2]] <- env$vpn[env$k]
		env$hslist[[k2]] <- hsplit
		env$gplist[[k2]] <- gap.prop
		env$k <- env$k+1
	}else{
		env$k <- env$k+1
		apply(tab,1,function(z){
			  flucplot(z,gap.prop, hsplit, env)
			  })
	}
#
	return(invisible(TRUE))
}

fluctree <- function(dims,parent, hsplit, gap.prop, env, ...){
	nv <- dims[1]
	dims <- dims[-1]
	
	
	if(hsplit[1]){
		w <- rep((1-gap.prop)/nv,nv)
		x <-  w/2 + 0:(nv-1)*(w + gap.prop/(nv-1)) 
		h <- rep(1,nv)
		y <- rep(0.5,nv)
	}else{
		h <- rep((1-gap.prop)/nv,nv)
		y <- h/2 + 0:(nv-1)*(h + gap.prop/(nv-1))
		w <- rep(1,nv)
		x <- rep(0.5,nv)
	}
#	if(length(hsplit) == 1){
#		hsplit <- !hsplit	
#}else{
		hsplit <- hsplit[-1]	#hsplit[-length(hsplit)]#	
#}
	
	
	if(length(dims) > 2){
		children <- vpList()
		for(i in 1:nv){
			env$k <- env$k+1
			tmp <- viewport( x[i],y[i],w[i],h[i],just="centre" , name = env$vpn[env$k])
			children[[i]] <- fluctree(dims,parent = tmp,hsplit, gap.prop, env)
			
		}
		return(vpTree(parent,children))
		
	}else{
		children <- vpList()
		for(i in 1:nv){
			env$k <- env$k+1
			children[[i]] <- viewport( x[i],y[i],w[i],h[i],just="centre" , name = env$vpn[env$k])
		}
		return(vpTree(parent,children))
	}
	
	
}

gridfluc <- function(tab, gap.prop = 0.1, maxv = NULL,vp = NULL, ...){
	
	n <- nrow(tab)
	m <- ncol(tab)
	w <- (1-gap.prop)/m
	h <- (1-gap.prop)/n
	x <-  w/2  
	y <- h/2 
	if(n > 1){
		y <- y + 0:(n-1)*(h+gap.prop/(n-1))	
	} 
	if(m > 1){
		x <- x + 0:(m-1)*(w+gap.prop/(m-1))
	}  
	if(is.null(maxv)){
		maxv <- max(tab)	
	}
	tab <- sqrt(tab/maxv)
	
	draw2( h, w,  t(replicate(n,x)) , replicate(m,y), border = NA,bg="lightgrey",vp=vp)
	ht <- as.matrix(h*tab)
	wt <- as.matrix(w*tab)
	
	draw2(ht, wt,  t(replicate(n,x)), replicate(m,y),border = NA,bg=rgb(0,0,0,alpha=0.7), vp=vp)
}


draw2 <- function (H, W, X, Y, alpha = 1, border = "black", bg = "white", 
vp = NULL) 
{
    grid.rect(x = unit(X, "npc"), y = unit(Y, "npc"), width = unit(W, 
																   "npc"), height = unit(H, "npc"), just = "centre", 
			  default.units = "npc", name = NULL, gp = gpar(col = border, 
															fill = bg, alpha = alpha), draw = TRUE, vp = vp)
}

addrect = function( vp ,breaks, col = "red", gap.prop = 0, rev.y = FALSE){
	yc <- breaks[[1]]
	xc <- breaks[[2]]
	
	nyc <- length(yc)
	nxc <- length(xc)
	stopifnot( nxc > 1 & nxc==nyc)
	n <- yc[nyc]
	m <- xc[nxc]
	pushViewport(vp)
	
	xc <- (xc-1)/(m-1)
	xc <- xc - gap.prop/(m-1)/2
	xc <- xc/(1 - gap.prop/(m-1))
	xc[1] <- 0
	xc[nxc] <- 1
	
	yc <- (yc-1)/(n-1)
	yc <- yc - gap.prop/(n-1)/2
	yc <- yc/(1 - gap.prop/(n-1))
	yc[1] <- 0
	yc[nyc] <- 1
	
	dyc <- diff(yc)
	dxc <- diff(xc)
	if(rev.y){
		mapply( function(x,y,w,h){
			   grid.rect(x,y,w,h,gp=gpar(fill=NA,col=col,lwd=2),just=c("left","top"))
			   }, x = as.list( xc[-nxc] ), y = as.list( 1-yc[-nyc] ),w = as.list(dxc),h = as.list(dyc))
		
	}else{
		mapply( function(x,y,w,h){
			   grid.rect(x,y,w,h,gp=gpar(fill=NA,col=col,lwd=2),just=c("left","bottom"))
			   }, x = as.list( xc[-nxc] ), y = as.list( yc[-nyc] ),w = as.list(dxc),h = as.list(dyc))
	}
	
	upViewport()
	return(invisible(TRUE))
	
}
kendalls <- function(M){
	cs <- colSums(M)
	rs <- rowSums(M)
	css <- rev(cumsum(rev(cs)))
	rss <- rev(cumsum(rev(rs)))
	n <- nrow(M)
	m <- ncol(M)
	N <- sum(M)
	x <- sum(cs*(css-cs)[1:m])
	y <- sum(rs*(rss-rs)[1:n])
	
	
	storage.mode(M) <- "integer"
	M2 <- M[,m:1]
	
	dims <- as.integer(c(n,m,1))
	crt1 <- .Call("simplecrit",M,dims,as.integer(1))
	crt2 <- .Call("simplecrit",M2,dims,as.integer(1))
	
	tau <- (crt2-crt1)/sqrt(x)/sqrt(y)
	scrt <- crt1/x/y*N^2
#return(c(tau,crt1,scrt,n,m,N))
	return(tau)
}




#optileplot = function(x, tau0 = NULL, col ="red",iter = 100, gap.prop = 0.2,presort=TRUE, floor = 0, critfun = "class", rev.y = TRUE,...){
#	stopifnot(inherits(x,"table"))
#	stopifnot( length(dim(x)) == 2 )
#	x = optile(x,iter=iter,presort=presort, critfun=critfun)
#	n = nrow(x)
#	m = ncol(x)
#	storage.mode(x) = "integer"
#	if(is.null(tau0)){
#		tau0 = kendalls(x)
#	}
#	if( floor > 0 ){
#		x2 = apply(x,1:2, function(z){
#				   ifelse(z < floor, 0, z)
#				   })
#		storage.mode(x2)="integer"
#		cuts = .Call("getclust",x2,as.integer(dim(x)),tau0)
#	}else{
#		cuts = .Call("getclust",x,as.integer(dim(x)),tau0)
#	}
#	
#	r = length(cuts)/2
#	if(rev.y){
#		x = x[nrow(x):1,]	
#	}
#	breaks = list(c(1,cuts[1:r]+1),c(1,cuts[(r+1):(r+r)]+1))
#		
#	bs = fluctile(x, gap.prop=gap.prop)
#	addrect(bs,breaks,col, gap.prop, rev.y = rev.y)
#	return(invisible(TRUE))
#}

#optileplot = function(x, tau0 = NULL, col ="red", gap.prop = 0.2, floor = 0, rev.y = TRUE,fun = "class", presort = FALSE, foreign = NULL, args = list(), perm.cat = TRUE, method = NULL, iter = 1, past = NULL, scale = FALSE, 
#    freqvar = NULL, return.data = TRUE, return.type = "matrix", vs = 0,...){
#	stopifnot(inherits(x,"table") | inherits(x,"matrix") )
#	stopifnot( length(dim(x)) == 2 )
#	x = optileXYZ(x, fun = fun, presort = presort, foreign = foreign, args = args, perm.cat = perm.cat, method = method, iter = iter, past = past, scale = scale, 
#    freqvar = freqvar, return.data = return.data, return.type = "matrix", vs = vs)
#	n = nrow(x)
#	m = ncol(x)
#	print("x opt = ")
#	print(x)
#	storage.mode(x) = "integer"
#	if(is.null(tau0)){
#		tau0 = kendalls(x)
#	}
#	if( floor > 0 ){
#		x2 = apply(x,1:2, function(z){
#				   ifelse(z < floor, 0, z)
#				   })
#		storage.mode(x2)="integer"
#		cuts = .Call("getclust",x2,as.integer(dim(x)),tau0)
#	}else{
#		cuts = .Call("getclust",x,as.integer(dim(x)),tau0)
#	}
#	
#	r = length(cuts)/2
#	if(rev.y){
#		x = x[nrow(x):1,]	
#	}
#	breaks = list(c(1,cuts[1:r]+1),c(1,cuts[(r+1):(r+r)]+1))
#		
#	bs = fluctile(x, gap.prop=gap.prop)
#	addrect(bs,breaks,col, gap.prop, rev.y = rev.y)
#	return(invisible(TRUE))
#}

cfluctile <- function(x, tau0 = NULL, col ="red", gap.prop = 0.2, floor = 0, rev.y = FALSE,...){
	stopifnot(inherits(x,"table") | inherits(x,"matrix") )
	stopifnot( length(dim(x)) == 2 )
	
	n <- nrow(x)
	m <- ncol(x)
	storage.mode(x) <- "integer"
	if(is.null(tau0)){
		tau0 <- kendalls(x)
	}
	if( floor > 0 ){
		x2 <- apply(x,1:2, function(z){
				   ifelse(z < floor, 0, z)
				   })
		storage.mode(x2) <- "integer"
		cuts <- .Call("getclust",x2,as.integer(dim(x)),tau0)
	}else{
		cuts <- .Call("getclust",x,as.integer(dim(x)),tau0)
	}
	
	r <- length(cuts)/2
	if(rev.y){
		x <- as.table(x[nrow(x):1,]	)
	}
	breaks <- list(c(1,cuts[(r+1):(r+r)]+1),c(1,cuts[1:r]+1))
	
	bs <- fluctile(x, gap.prop=gap.prop)
	addrect(bs,breaks,col, gap.prop, rev.y = !rev.y)
	return(invisible(TRUE))
}



bloma <- function(k, ncr = 12, ar = 1){
	Mlist = list()
	N=0
	M=0
	for(i in 1:k){
		n <- max(2,rpois(1,ncr))
		m <- max(2,rpois(1,ar*ncr))
		Mlist[[i]] <- extracat:::smat2(n,m,rpois(1,500))
		N <- N+n
		M <- M+m
	}
	S <- matrix(0,ncol=M,nrow=N)
	xind <- 1:M
	yind <- 1:N
	
	for(i in 1:k){
		MS <- Mlist[[i]]
		
		
		xi <- sample(xind, size = ncol(MS))	
		yi <- sample(yind, size = nrow(MS))
		S[yi,xi] <- MS
		xind <- xind[-which(xind %in% xi)]
		yind <- yind[-which(yind %in% yi)]		
	}
	return(S)
}


	
scaledclass <- function(M, vs = 1, ... ){
	dims <- as.integer(c(dim(M),1))
	stopifnot( length(dims) == 3 )
	
	MI <- outer(rowSums(M),colSums(M))/sum(M)
	storage.mode(M) <- "integer"
	storage.mode(MI) <- "integer"
	crt <- .Call("simplecrit",M, dims, as.integer(vs))	
	crt0 <- .Call("simplecrit",MI, dims, as.integer(vs))
	return(crt/crt0)
}
