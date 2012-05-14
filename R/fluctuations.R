

fluctile <- function(tab, dir = "b", just = "c", hsplit = FALSE,shape ="r", gap.prop = 0.1, border = NULL, label = TRUE, lab.opt = list(abbrev = 3, lab.cex = 1.2), add = FALSE, tile.col = hsv(0.1,0.1,0.1,alpha=0.6), bg.col = NULL, ...  ){
	#tab <- t(tab)
	dm <- dim(tab)
	maxv <- max(tab)
	nd <- length(dm)
	
    if( dir %in% c("vh","hv","b","both")){
        dir <- "b"
    }
    if( dir %in% c("h","horizontal")){
        dir <- "h"
    }
    if( dir %in% c("v","vertical")){
        dir <- "v"
    }
    if(nchar(just)[1] == 2){
        just <- c(substr(just[1],1,1),substr(just[1],2,2))
    }
    just <- sapply(just, function(j){
        switch(j, t="top", b = "bottom", c = "center", r = "right", l = "left", NULL)
    })
    if("abbrev" %in% names(lab.opt)){
        abbrev <- lab.opt$abbrev
    }else{
        abbrev <- 3
    }
    if("lab.cex" %in% names(lab.opt)){
        lab.cex <- lab.opt$lab.cex
    }else{
        lab.cex = 1.2
    }
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
	if(!add){
		#dev.new()
		grid.newpage()
	}
	vp0 = viewport(x = border[1]*nx/(nx+1) + (1-border[1])/2 , y = border[2]/(ny+1) + (1-border[2])/2, width = 1-border[1], height = 1-border[2],name="base")
	if(length(dm)==2){
		
		pushViewport(vp0)
		if(!add){
		#grid.newpage()
			grid.rect(gp = gpar(fill=rgb(0,0,0,alpha=0.1),col=NA))
		}
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
		gridfluc(tab,dir,just,shape,gap.prop, maxv=maxv, bg = bg.col, col = tile.col)
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
			grid.text(rn, x = 0.5, y = yc, gp = gpar(cex=lab.cex), rot = 65)
			popViewport()
			
			vpC <- viewport(x = border[1]/2 + (1-border[2])/2, y = 1-border[2]/4, width = 1-border[1], height = border[2]/2,name="collabs")
			pushViewport(vpC)
			grid.text(cn, x = xc, y = 0.5,rot = 25, gp = gpar(cex=lab.cex))
			
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
			   
			   gridfluc(x,dir,just,shape,gp, maxv=maxv, bg = bg.col, col = tile.col)
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
            if(abbrev != FALSE){
                labs <- sapply(labs, function(l) abbreviate(l,abbrev))
            }
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

gridfluc <- function(tab,dir = "b", just = "c", shape = "r", gap.prop = 0.1, maxv = NULL,vp = NULL, col = NULL, bg = NULL, ...){
	
	if( is.null(bg) ){
		bg <- "lightgrey"
	}
	if( is.null(col) ){
		col <- hsv(0,0,0,alpha=0.7)
	}
	
	n <- nrow(tab)
	m <- ncol(tab)
	w <- (1-gap.prop)/m
	h <- (1-gap.prop)/n
#centered x- and y-coords:
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
	
	
	draw2( h, w,  t(replicate(n,x)) , replicate(m,y), border = NA,bg=bg,vp=vp, just = "centre")
	
    if("right" %in% just){
        x <- x+w/2
    }
    if("left" %in% just){
        x <- x-w/2
    }
    if("top" %in% just){
        y <- y+h/2
    }
    if("bottom" %in% just){
        y <- y-h/2
    }
    if(dir == "b"){
        tab <- sqrt(tab/maxv)
        ht <- as.matrix(h*tab)
        wt <- as.matrix(w*tab)
	}
    if(dir == "h"){
        tab <- tab/maxv
        wt <- as.matrix(w*tab^0)
        ht <- as.matrix(h*tab)
    }
    if(dir == "v"){
        tab <- tab/maxv
        ht <- as.matrix(h*tab^0)
        wt <- as.matrix(w*tab)
    }
if(dir == "n"){
    tab <- tab/maxv
    ht <- matrix(0,ncol=ncol(tab), nrow=nrow(tab))
    wt <- ht
    ht[tab > 0] <- h
    wt[tab > 0] <- w
}
	if(shape == "r"){
		draw2(ht, wt,  t(replicate(n,x)), replicate(m,y), border = NA,bg=col, vp=vp, just = just)
	}else{
		if(shape == "c"){
			angles <- seq(0,360,0.5)/180*pi
		}
		if(shape == "o"){
			angles <- seq(22.5,360,45)/180*pi
		}
		if(shape == "d"){
			angles <- seq(0,360,90)/180*pi
		}
		corners <-  cbind( cos(angles), sin(angles))
		mapply(function(x,y,h,w){
			   grid.polygon(x = x+corners[,1]*w/2, y = y+corners[,2]*h/2, gp = gpar(fill= col, col = NA))
			   }, x = t(replicate(n,x)), y = replicate(m,y), h = ht, w = wt)
	
	}
	
}


draw2 <- function (H, W, X, Y, alpha = 1, border = "black", bg = "white", 
vp = NULL, just = "centre") 
{
grid.rect(x = unit(X, "npc"), y = unit(Y, "npc"), width = unit(W, 
"npc"), height = unit(H, "npc"), just = just, default.units = "npc", 
name = NULL, gp = gpar(col = border, fill = bg, alpha = alpha), 
draw = TRUE, vp = vp)
}
addrect = function( vp ,breaks, col = "red", lwd = 2, lty = 1, gap.prop = 0, rev.y = FALSE){
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
			   grid.rect(x,y,w,h,gp=gpar(fill=NA,col=col,lwd=lwd, lty = lty),just=c("left","top"))
			   }, x = as.list( xc[-nxc] ), y = as.list( 1-yc[-nyc] ),w = as.list(dxc),h = as.list(dyc))
		
	}else{
		mapply( function(x,y,w,h){
			   grid.rect(x,y,w,h,gp=gpar(fill=NA,col=col,lwd=lwd, lty = lty),just=c("left","bottom"))
			   }, x = as.list( xc[-nxc] ), y = as.list( yc[-nyc] ),w = as.list(dxc),h = as.list(dyc))
	}
	
	upViewport()
	return(invisible(TRUE))
	
}



cfluctile <- function(x, tau0 = NULL, method="Kendall", col ="red", lwd = 2, lty = 1, gap.prop = 0.2, 
floor = 0, rev.y = FALSE, add = FALSE, shape = "r", just = "c", dir = "b", ...){
	stopifnot(inherits(x,"table") | inherits(x,"matrix") )
	stopifnot( length(dim(x)) == 2 )
    
    if(kendalls(x) < 0){
        x <- x[,ncol(x):1]
        print("Reversed column category order...")
    }

	if( method %in% c("Kendall","kendall","tau",1) ){
		method <- as.integer(1)	
	}
	if( method %in% c("kappa","Cohen",2) ){
		method <- as.integer(2)	
	}
	
	n <- nrow(x)
	m <- ncol(x)
	storage.mode(x) <- "integer"
	if(is.null(tau0)){
		if(method == 1){
			tau0 <- kendalls(x)
		}else{
			tau0 <- 0	
		}
	}
	if( floor > 0 ){
		x2 <- apply(x,1:2, function(z){
				   ifelse(z < floor, 0, z)
				   })
		storage.mode(x2) <- "integer"
		cuts <- .Call("getclust",x2,as.integer(dim(x)),tau0,method)
	}else{
		cuts <- .Call("getclust",x,as.integer(dim(x)),tau0,method)
	}

	r <- length(cuts)/2
	if(rev.y){
		x <- as.table(x[nrow(x):1,]	)
	}
	breaks <- list(c(1,cuts[1:r]+1),c(1,cuts[(r+1):(r+r)]+1))
	if(!add){
		bs <- fluctile(x, gap.prop=gap.prop, shape = shape, just = just, dir = dir)
	}else{
		bs <- fluctile(x, gap.prop=gap.prop, bg.col=rgb(0,0,0,alpha=0),tile.col = rgb(0,0,0,alpha=0), add = TRUE, shape = shape, just = just, dir = dir)	
	}
	addrect(bs,breaks,col, gap.prop, rev.y = !rev.y, lwd = lwd, lty = lty)
	return(invisible(TRUE))
}



fluctile3d = function(x, shape = "cube", col = "darkgrey", alpha = 0.8,...){
	if(!"rgl" %in% .packages(all.available = TRUE)){
		cat("Please install the package 'rgl' to run this function")
		return(invisible(TRUE))
	}
	check <- tryCatch(require(rgl), error = function(e) FALSE)
	if(!check){ 
		cat("Problems with package 'rgl' occured...")
		return(invisible(TRUE))
	}
	if(shape %in% c("o","oct","octahedron")){
        shape <- "octagon"
	}
    if(!is.data.frame(x)){
		x2 <- subtable(as.data.frame(as.table(x)),1:3)	
		dim <- dim(x)
	}else{
		x2 <- subtable(x,1:3)
		dim <- sapply(x,nlevels)[1:3]	
	}
	maxfreq <- max(x2$Freq)
	lvls <- lapply(x2, levels)
	x2 <- sapply(x2, as.integer)
	
	#MX <- identityMatrix()
	MX <- matrix(0,ncol=4,nrow=4)
	diag(MX) <- 1
	open3d()
#	lines3d(x = c(0.5,dim[1]+0.5), y = c(0.5,0.5), z =  c(0.5,0.5))
#	lines3d(x = c(0.5,0.5), y = c(0.5,dim[2]+0.5), z =  c(0.5,0.5))
#	lines3d(x = c(0.5,0.5), y = c(0.5,0.5), z =  c(0.5,dim[3]+0.5))
#	lines3d(x = c(dim[1]+0.5,0.5), y = c(dim[2]+0.5,dim[2]+0.5), z =  c(dim[3]+0.5,dim[3]+0.5))
#	lines3d(x = c(dim[1]+0.5,dim[1]+0.5), y = c(dim[2]+0.5,0.5), z =  c(dim[3]+0.5,dim[3]+0.5))
#	lines3d(x = c(dim[1]+0.5,dim[1]+0.5), y = c(dim[2]+0.5,dim[2]+0.5), z =  c(dim[3]+0.5,0.5))
	wire3d( translate3d( cube3d( trans = scaleMatrix(dim[1]/2,dim[2]/2,dim[3]/2)), (dim[1]+1)/2,(dim[2]+1)/2,(dim[3]+1)/2)) 
	dot3d( translate3d( cube3d( trans = scaleMatrix(dim[1]/2,dim[2]/2,dim[3]/2)), (dim[1]+1)/2,(dim[2]+1)/2,(dim[3]+1)/2))
	
	
	text3d( x = 1:dim[1], y = rep(0, dim[1]), z = rep(0, dim[1]), texts = lvls[[1]])
	text3d( x = rep(0, dim[2]), y = 1:dim[2], z = rep(0, dim[2]), texts = lvls[[2]])
	text3d( x = rep(0, dim[3]), y = rep(0, dim[3]), z = 1:dim[3], texts = lvls[[3]])
	
	text3d( x = 1:dim[1], y = rep(dim[2]+1, dim[1]), z = rep(dim[3]+1, dim[1]), texts = lvls[[1]])
	text3d( x = rep(dim[1]+1, dim[2]), y = 1:dim[2], z = rep(dim[3]+1, dim[2]), texts = lvls[[2]])
	text3d( x = rep(dim[1]+1, dim[3]), y = rep(dim[2]+1, dim[3]), z = 1:dim[3], texts = lvls[[3]])
	
	apply(x2,1,function(z){
		s <- ((z[4]/maxfreq)^(1/3))/2
		  if(shape == "cube"){
			shade3d( translate3d( cube3d(col=col, trans = scaleMatrix(s,s,s)), z[1], z[2], z[3]), alpha = alpha )	
		  }
		  if(shape == "octagon"){
			 shade3d( translate3d( octahedron3d(col=col, trans = scaleMatrix(s,s,s)), z[1], z[2], z[3]), alpha = alpha )	
		  }

	})
	
	
	return(invisible(TRUE))
	
}
cfluctile3d = function(x, xc = NULL, yc = NULL, zc = NULL, shape ="cube", col = c("darkgrey","red"), alpha = c(0.8,0.2),...){
															
if(shape %in% c("o","oct","octahedron")){
	shape <- "octagon"
}
	if(is.null(xc)){
		xc <- attr(x,"xc")
		if(is.null(xc)){
			xc <- as.integer(attr(x,"grps")[[1]])
		}	
	}
	if(is.null(yc)){
		yc <- attr(x,"yc")	
		if(is.null(yc)){
			yc <- as.integer(attr(x,"grps")[[2]])
		}
	}
	if(is.null(zc)){
		zc <- attr(x,"zc")
		if(is.null(zc)){
			zc <- as.integer(attr(x,"grps")[[3]])
		}	
	}
	stopifnot( all( c(is.null(xc),is.null(yc),is.null(zc)) == FALSE))
	
	colA <- col[1]
	colB <- ifelse(length(col) > 1, col[2], "red")
	alphaA <- alpha[1]
	alphaB <- ifelse(length(alpha) > 1, alpha[2], 0.2)
	
	if(!is.data.frame(x)){
		dim <- dim(x)
		x <- subtable(as.data.frame(as.table(x)),1:3,allfactor=TRUE)	
		
	}else{
		dim <- sapply(x,nlevels)[1:3]
		x <- subtable(x,1:3,allfactor=TRUE)
			
	}
	ncl <- nlevels(as.factor(xc))
	xc <- lapply(1:ncl, function(v){
			which(xc == v)
		})
	yc <- lapply(1:ncl, function(v){
			which(yc == v)
		})
	zc <- lapply(1:ncl, function(v){
			which(zc == v)
		})
		
	stopifnot(length(xc) == length(yc) & length(yc) == length(zc))
	fluctile3d(x, col = colA, shape = shape, alpha = alphaA)
	mapply(function(s1,s2,s3){
		r1 <- diff(range(s1))+1
		c1 <- mean(range(s1))
		r2 <- diff(range(s2))+1
		c2 <- mean(range(s2))
		r3 <- diff(range(s3))+1
		c3 <- mean(range(s3))
		  shade3d( translate3d( cube3d(col=colB, trans = scaleMatrix(r1/1.98,r2/1.98,r3/1.98)), c1, c2, c3), alpha = alphaB )	
		
	}, s1 = xc, s2 = yc, s3 = zc)
		return(invisible(TRUE))
}
