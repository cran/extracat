
 
facetshade <- function( data, mapping = NULL, f, bg.all = TRUE, keep.orig = FALSE, ... ){

body = ggplot

	# if(inherits(data,"ggplot")){
		# p <- data
		# if(is.null(mapping)){
			# mapping <- p$mapping
			# if(length(mapping) == 0) mapping <- p$layer[[1]]$mapping
		# }
		# data <- data$data
	#}else{
		p <- body(mapping = mapping, ... )
	#}
	gnames <- colnames(attr(terms(f),"factors"))
	
	ind <- which(names(data) %in% gnames)
	

	ord <- do.call(order,data[,ind,drop=FALSE])
	data <- data[ord,]
	
	mdata <- subtable(data,ind)
	gs <- mdata$Freq
	ng <- nrow(mdata)
	mdata <- mdata[rep(1:ng, each=nrow(data)),]
	
	
	xn <- toString(mapping$x)
	yn <- toString(mapping$y)
	
	
	mdata[xn] <- data[xn]
	
	if(yn != ""){
		mdata[yn] <- data[yn]
	}	

	if(!is.null(mapping$colour)){
		cn <- toString(mapping$colour)
		if(!cn %in% names(mdata)){
			mdata[cn] <- rep(unlist(data[cn]),ng)
		}
	}
	if(!is.null(mapping$group)){
		gn <- toString(mapping$group)
		if(!gn %in% names(mdata)){
			mdata[gn] <- rep(unlist(data[gn]),ng)
		}
	}
	if(!is.null(mapping$fill)){
		fn <- toString(mapping$fill)
		if(!fn %in% names(mdata)){
			mdata[fn] <- rep(unlist(data[fn]),ng)
		}
	}
	if(keep.orig){
		for(i in gnames){
			mdata[paste("orig",i,sep=".")] <- rep(unlist(data[i]),ng)
		}
	}


	if(!bg.all){
		#mdata <- subset(mdata,gv == 0)
		cgs <- cumsum(gs)
		cgs <- cbind(c(1,cgs[-ng]+1),cgs) + (0:(ng-1))*nrow(data)
		rm.ind <- unlist( apply(cgs,1,function(z) z[1]:z[2]) )
		mdata <- mdata[-rm.ind,]
	}
	

p <- p %+% mdata + facet_grid(f) + guides(colour = guide_legend(title=NULL))

return(p)

}
