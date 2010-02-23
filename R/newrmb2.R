
# -------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------- #

draw = function(H,W,X,Y,alpha = 1,border = "black",bg = "white",vp=NULL){
	grid.rect(x =  unit(X, "npc"), y = unit(Y, "npc"),
          width = unit(W, "npc"), height = unit(H, "npc"),
          just = c("left","bottom"),
          default.units = "npc", name = NULL,
          gp=gpar(col=border,fill = bg, alpha = alpha), draw = TRUE, vp = vp
	)
}


rmb = function(f, dset, hsplit = NULL, spine = FALSE,  hlcat = 1,  eqwidth = FALSE,
                tfreq = "Id", max.scale = 1,  use.na = FALSE, expected = NULL,
                mod.type = "poisson", resid.type = "pearson",
                use.expected.values = FALSE, resid.max = NULL, cut.rv = TRUE, cut.rs = 5,                
                base = 0.2, mult = 1.5, colv = NULL, lab = TRUE, yaxis = TRUE, label = NULL,
                min.alpha = 0.1, base.alpha = 0.75, boxes = TRUE, lab.tv = FALSE,
                varnames = TRUE, abbr = FALSE, lab.cex = 1.2,...){

# ----- parameter check 1 ----- #
		
		dset = data.frame(dset)
		f = as.formula(f)
		
		if(max.scale > 1 ){
			stop(simpleError("Wrong max.scale specification!"))
		}
				
		if(!( (base < 1) & (base > 0) ) ){
			stop(simpleError("Wrong base specification!"))
		}
		if( !((min.alpha <= 1) & (min.alpha >= 0)) ){
			stop(simpleError("Wrong alpha specification!"))
		}
		if(!(   (eqwidth == T)|(tfreq %in% c("Id","sqrt","log")) )){
			stop(simpleError("Wrong eqwidth/tfreq specification!"))
		}
		

# ----- terms extraction ----- #
	fterms = attr(terms(f),"term.labels")
	nv = length(fterms)
		
	nv0 = nv
	if(is.null(hsplit)){
		hsplit = rep(c(T,F),nv)[1:nv]	
	}
		
	hsplit[nv] = T
	tv = fterms[nv]
	
	ind =  match(fterms,names(dset))
# 	>>> a dummy-variable is used to handle the case of no row-variables
	if( sum(!hsplit) < 1 ){
		dset = data.frame(probability = rep("prob",nrow(dset)),dset)
		nv=nv+1
		hsplit = c(F,hsplit)
		ind = c(1,ind+1)
	}
	if( sum(hsplit) < 2 ){
		dset = data.frame(distribution = rep(tv,nrow(dset)),dset)
		nv=nv+1
		hsplit = c(T,hsplit)
		ind = c(1,ind+1)
		if("probability" %in% names(dset)){ expected = NULL }
	}

# 	>>> preparing the labels for plotting including the abbreviation option	
	if(!abbr){
		rclabs = lapply(dset[,ind],function(x) levels(as.factor(x)))
	}else{
		rclabs = lapply(dset[,ind],function(x) abbreviate(levels(as.factor(x)),3))
	}
	lab.tv = ifelse(spine,FALSE,lab.tv)

	if(is.null(label)){
		label = rep(T,nv)	
	}
	label[nv] = lab.tv
	if( sum(label) == 0 ){
		lab = FALSE	
	}
	rlabs = lapply(which( (!hsplit)*label > 0),function(x) rclabs[[x]])
	clabs = lapply(which(hsplit*label > 0),function(x) rclabs[[x]])
	
#	orig.labs = lapply(dset[,ind],function(x) levels(as.factor(x)))
	if( length(terms(f)) > 2 ){
		names(dset)[which( suppressWarnings(names(dset) == terms(f)[[2]]))] = "Freq"
	}
	dset = subtable(dset,ind,keep.zero=F)

	if(length(hlcat) > 1 & !spine){
		#levels(dset[,nv])[which(! 1:nv %in% hlcat )] = NA 
		levels(dset[,nv])[-hlcat] = NA
		
		dset[,nv] = factor(dset[,nv],levels=levels(dset[,nv])[rank(hlcat)])
		if(lab.tv){
			clabs[[length(clabs)]] = levels(dset[,nv])	
		}
	}
	
	if(!use.na){
		dset = na.omit(dset)
	}else{
		for( i in 1:(ncol(dset)-1) ){
			ind = is.na(dset[,i])
			if(any(ind)){
				levels(dset[,i]) = c(levels(dset[,i]),"N/A")	
				dset[which(ind),i] = "N/A"
			}	
		}	
	}
	
	
# ----- descriptive parameters ----- #	
	ntc = nlevels(dset[,nv])
	ntc0 = ntc
	nlvl = sapply(dset[,(sapply(dset,class)=="factor")], nlevels)

	col.nlvl = nlvl[hsplit]
	row.nlvl = nlvl[!hsplit]								

	nc = prod(col.nlvl)
	nr = prod(row.nlvl)

	ncv = length(col.nlvl)
	nrv = length(row.nlvl)

# ----- parameter check 2 ----- #

		#if(max(hlcat) > ntc){
		#	stop(simpleError("Wrong hlcat specification!"))
		#}
		if(suppressWarnings(max(unlist(expected))) > (nv - (mod.type=="logit")) ){
			stop(simpleError("Wrong expected specification!"))
		}
# ----- restoring original variable orders ----- #			
#	for( i in 1:nv ){
#		dset[,i] = factor(dset[,i],levels = orig.labs[[i]])	
#	}

# ----- computation of the underlying (relative) frequencies ----- #		
	tt1 = ftable(tapply(dset$Freq,as.list(dset[,1:nv]),sum),col.vars=which(hsplit))
	tt2 = spread(ftable(tapply(dset$Freq,as.list(dset[,1:(nv-1)]),sum),col.vars=which(hsplit[1:(nv-1)])),ncol=ntc)

	tt1[is.na(tt1)] = 0 
	tt2[is.na(tt2)] = 1 

	H0 = tt1/tt2
	H = H0

# ----- a few more auxiliary variables ----- #	
	nlvl = sapply(dset[,(sapply(dset,class)=="factor")], nlevels) # muss nicht sein
	
	col.nlvl = nlvl[hsplit]
	row.nlvl = nlvl[!hsplit]								

	nc = prod(col.nlvl)
	nr = prod(row.nlvl)

	row.base = ifelse(nr > 1, base * min(1,nrv/ncv), 0) 
	col.base = ifelse(nc > ntc, base * min(1,ncv/nrv), 0) 

	ncv = length(col.nlvl)
	nrv = length(row.nlvl)
	
	rind = which( (!hsplit)*label > 0)
	cind = which(hsplit*label > 0)
	nrl = length(rind)*lab
	ncl = length(cind)*lab

# ----- model computation in expected mode ----- #	
	if(!is.null(expected) ){
		int.dset = lapply(dset,as.integer)
		dset = dset[do.call("order",int.dset),]
				
		if(!(mod.type %in% c("poisson","polr"))){
			stop(simpleError("Wrong mod.type specification!"))
		}
		# ------ logit response residuals ----- #
		if(mod.type == "polr"){ 
			if(resid.type != "response"){
				print(simpleWarning("Argument resid.type ignored. Only response residuals are implemented for polr models."))
				resid.type="response"	
			}
			nameslist = names(dset)[1:(nv-1)]		
			nameslist = nameslist[which(  !(nameslist %in% c("probability","distribution"))  )]
			single.terms = do.call("paste",c(as.list( nameslist ),sep="+"))
			interaction.terms = do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(nameslist[x]),sep="*"))),sep="+"))
			full.terms = do.call("paste",c(as.list(nameslist),sep="*"))
			
			#single.terms = do.call("paste",c(as.list(names(dset)[1:(nv-1)]),sep="+"))
			#full.terms = do.call("paste",c(as.list(names(dset)[1:(nv-1)]),sep="*"))
			#interaction.terms = do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(names(dset)[x]),sep="*"))),sep="+"))
			mod.formula = as.formula(paste(fterms[nv0],"~",single.terms,"+",interaction.terms,sep=""))
			modS.formula = as.formula(paste(fterms[nv0],"~",full.terms,sep=""))	
			print(head(dset))
			print(mod.formula)
			print(modS.formula)
			modS = polr( formula = as.formula(modS.formula),data = dset, weights = dset$Freq, method ="logistic")
			mod = polr( formula = as.formula(mod.formula),data = dset, weights = dset$Freq, method ="logistic")
			pred = predict(mod,newdata = subtable(dset,c(1:(nv-1)),keep.zero=F),type="probs")
			fit.values = as.vector(t(pred))
			tt1r = ftable(tapply(fit.values,as.list(dset[,1:nv]),sum),col.vars=which(hsplit))
			resid.mat = (H-tt1r)
			p.value=anova(mod,modS)[2,7]
			
			if( use.expected.values ){
				H0 = tt1r/tt2
				H = H0
			}
			
		}
		# ----- poisson model ----- #	
		if(mod.type == "poisson"){
			nameslist = names(dset)[1:(nv)]		
			nameslist = nameslist[which(  !(nameslist %in% c("probability","distribution"))  )]
			single.terms = do.call("paste",c(as.list( nameslist ),sep="+"))
			interaction.terms = do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(nameslist[x]),sep="*"))),sep="+"))
						
			#single.terms = do.call("paste",c(as.list(names(dset)[1:(nv-1)]),sep="+"))
			#interaction.terms = do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(names(dset)[x]),sep="*"))),sep="+"))
			mod.formula = paste("Freq~",single.terms,"+",interaction.terms,sep="")
			mod = glm( as.formula(mod.formula),data = dset,  family = "poisson")
		    mod.residuals = residuals(mod,type=resid.type)
			resid.mat =  ftable(tapply(mod.residuals,as.list(dset[,1:nv]),sum),col.vars=which(hsplit))
			
			Df  = mod$df.residual
			Dev = mod$deviance
			p.value = 1-pchisq(Dev,Df)
			
			if( use.expected.values ){
				pred = predict(mod,type="response")
				pred.mat =  ftable(tapply(pred,as.list(dset[,1:nv]),sum),col.vars=which(hsplit))
				H0 = pred.mat/tt2
				H = H0
			}
			
		}
	}

# ----- modify parameters for rmbplot mode ----- #	
	if( spine ){
		nlvl[nv] = 1
		nc = nc/ntc
		ntc = 1
	}

# ----- computation of rectangle widths, heights and x/y coordinates in matrix form ----- #	
	W = matrix(1,ncol=nc,nrow=nr)
	S = space(nlvl,hsplit,base=base,mult=mult,last.col.zero = T,last.row.zero = F)

	x = ( seq(0,1-1/nc,1/nc)   )*(1-col.base)
	if(nc > nlevels(dset[,nv])){
		x = x + c(0,cumsum(S[[1]]))*col.base/sum(S[[1]])
	}
	if( suppressWarnings(any(S[[2]])) ){
		y = (seq(0,1-1/nr,1/nr))*(1-row.base) + c(0,cumsum(S[[2]]))*row.base/sum(S[[2]])
	}else{
		y = 0	
	}
		
	X = spread(t(matrix( x )),nrow=nr)
	Y0 = spread(matrix( y ),ncol=nc)
	Y = Y0
	X2 = X[,seq(1,nc,ntc)]
	Y2 = Y[,seq(1,nc,ntc)]


	H2 = matrix(1,ncol=nc/ntc,nrow=nr)
	tt3 = ftable(tapply(dset$Freq,as.list(dset[,1:(nv-1)]),sum),col.vars=which(hsplit[1:(nv-1)]))
	
	tt3[is.na(tt3)] = 0
	W2 = tt3 / max(tt3)

	if(nc > nlevels(dset[,nv])){
		width.cor = (ntc-1)*S[[1]][1]/sum(S[[1]])*col.base
	}else{
		width.cor = 0
	}
# ----- ---------------- ----- #	
# ----- plotting section ----- #	
# ----- ---------------- ----- #

# ----- opening the basic viewport ----- #	
	#dev.new()
	#vp0 <- viewport(x = 0.02+yaxis*0.02, y = 0.5, w = 0.96-yaxis*0.04 - 0.1*(!is.null(expected)), h = 0.96, just = c("left","centre"), name = "vp0")
	#pushViewport(vp0)
	#if(lab){
	#	vp1 <- viewport(x = 1-yaxis*0.02, y = 0, w = 1-(nrl+yaxis*2/3)*0.06, h = 1-ncl*0.06, just = c("right", "bottom"), name = "vp1")
	#	pushViewport(vp1)
	#}
	s0 = 0.02 # border
	s1 = 0.04 # yaxis
	s2 = 0.06 # labs
	s3 = 0.1  # expected scale
	dev.new()
	vp0 <- viewport(x = 0.5, y = 0.5, w = 1- 2*s0 , h = 1 - 2*s0, just = "centre", name = "vp0")
	pushViewport(vp0)
	vp1 <- viewport(x = 1-yaxis*s1 - s3*(!is.null(expected)) , y = 0, w = 1- nrl*s2 - yaxis*2*s1 - s3*(!is.null(expected)), h = 1-ncl*s2, just = c("right", "bottom"), name = "vp1")
	pushViewport(vp1)
	
	alpha = spread(W2,ncol=ntc)*base.alpha
	alpha[alpha < min.alpha] = min.alpha
	if(!is.null(expected)){
		alpha[alpha < 1] = 1	
	}

# ----- preparing the color matrix ----- #	
			if(is.null(expected)){
				if( is.null(colv)|length(colv) < ntc0 ){
					colv = rainbow(ntc0,alpha=base.alpha)
					if((!spine) & eqwidth){
						colv = rep(rainbow(1,alpha=base.alpha),ntc)
					}
				}else{
					colv = colv[1:ntc0]	
				}					
				col.Mat = spread( t(rep(colv , nc/ntc)),nrow=nr)
			}else{
				resid.mat[is.nan(resid.mat)] = 0
				resid.mat[is.na(resid.mat)] = 0
				if(is.null(resid.max)){
					resid.max = max(abs(resid.mat),na.rm=T)	
				}
				resid.mat[resid.mat < -resid.max] = -resid.max
				resid.mat[resid.mat > resid.max]  =  resid.max
				rmax = ifelse(abs(resid.max) < 10^{-6},1,resid.max)
				col.Mat = apply(resid.mat/rmax,2,function(x){
						qx = abs(x)
						if(cut.rv){ qx =  ceiling(qx * cut.rs)/cut.rs }
					rgb(sign(x) < 0,0,sign(x) > 0,alpha = qx*base.alpha )
					})
			}

# ----- plotting itself ----- #
# 	>>> divided into eqwidth = T | eqwidth = F >>> tfreq = Id/sqrt/log
# 	>>> in spine mode, the target categories will be plotted one after another (see the for-loop)
# 	>>> the 3 diff. draw-function add the rectangles for the relative frequencies, the weights and the background respectively
	if(eqwidth){
		draw(H2/nrow(H2)*(1-row.base),W2/ncol(W2)*(1-col.base) + (round(W2,digits=7) > 0)*width.cor,X2,Y2,alpha = 0.25, bg = "grey")
		draw(H2/nrow(H2)*(1-row.base),H2/ncol(H2)*(1-col.base) + width.cor,X2,Y2,alpha = 0.1, bg = "grey")
		if( spine ){
			for( i in 1:length(hlcat) ){
				if(i > 1){ 
					Y = Y + H/nrow(H)*(1-row.base)
				}
				xind = seq( hlcat[i], nc*ntc0, ntc0 )
				H = H0[,xind]
			 	draw(H/nrow(H)*(1-row.base),W/ncol(W)*(1-col.base),X,Y,alpha = 1, bg = col.Mat[,xind])
			}
			
			Y = Y0
		}else{
			H = apply(H,c(1,2),function(x) return(min(x/max.scale,1)))
			draw(H/nrow(H)*(1-row.base),W/ncol(W)*(1-col.base),X,Y,alpha = alpha, bg = col.Mat)
		}
	}else{
		
		if(tfreq == "Id"){
			W3 = spread(W2,ncol=ntc)

			if( nr > 1 ){
				X3 = spread(X[,seq(1,nc,ntc)],ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.base) 
			}else{
				X3 = spread(t(X[,seq(1,nc,ntc)]),ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.base) 
			}
			draw(H2/nrow(H2)*(1-row.base),W2/ncol(W2)*(1-col.base) + (round(W2,digits=7) > 0)*width.cor,X2,Y2,alpha = 0.25, bg = "grey")
			draw(H2/nrow(H2)*(1-row.base),H2/ncol(H2)*(1-col.base) + width.cor,X2,Y2,alpha = 0.1, bg = "grey")
			if( spine ){
				for( i in 1:length(hlcat) ){
					if(i > 1){ 
						Y = Y + H/nrow(H)*(1-row.base)
					}
					xind = seq( hlcat[i], nc*ntc0, ntc0 )
					H = as.matrix(H0[,xind])
					if( nr < 2 ){
						H = t(H)	
					}
					draw(H/nrow(H)*(1-row.base),W3/ncol(W3)*(1-col.base),X3,Y,alpha = 1, bg = col.Mat[,xind])
				}
				Y = Y0
			}else{
				H = apply(H,c(1,2),function(x) return(min(x/max.scale,1)))
				draw(H/nrow(H)*(1-row.base),W3/ncol(W3)*(1-col.base),X3,Y,alpha = 1, bg = col.Mat)
			}
		}
		if(tfreq == "log"){
			W2log = log(tt3+1)/max(log(tt3+1))
			W3 = spread(W2log,ncol=ntc)
			if( nr > 1 ){
				X3 = spread(X[,seq(1,nc,ntc)],ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.base) 
			}else{
				X3 = spread(t(X[,seq(1,nc,ntc)]),ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.base) 
			}
			draw(H2/nrow(H2)*(1-row.base),W2log/ncol(W2log)*(1-col.base) + (round(W2log,digits=7) > 0)*width.cor,X2,Y2,alpha = 0.25, bg = "grey")
			draw(H2/nrow(H2)*(1-row.base),H2/ncol(H2)*(1-col.base) + width.cor,X2,Y2,alpha = 0.1, bg = "grey")
			if( spine ){
				for( i in 1:length(hlcat) ){
					if(i > 1){ 
						Y = Y + H/nrow(H)*(1-row.base)
					}
					xind = seq( hlcat[i], nc*ntc0, ntc0 )
					H = as.matrix(H0[,xind])
					if( nr < 2 ){
						H = t(H)	
					}
					draw(H/nrow(H)*(1-row.base),W3/ncol(W3)*(1-col.base),X3,Y,alpha = 1, bg = col.Mat[,xind])
				}
				Y = Y0
			}else{
				H = apply(H,c(1,2),function(x) return(min(x/max.scale,1)))
				draw(H/nrow(H)*(1-row.base),W3/ncol(W3)*(1-col.base),X3,Y,alpha = 1, bg = col.Mat)
			}
		}
		if(tfreq == "sqrt"){
			W2sqrt = sqrt(tt3)/max(sqrt(tt3))
			W3 = spread(W2sqrt,ncol=ntc)
			if( nr > 1 ){
				X3 = spread(X[,seq(1,nc,ntc)],ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.base) 
			}else{
				X3 = spread(t(X[,seq(1,nc,ntc)]),ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.base) 
			}
			if( spine ){
				for( i in 1:length(hlcat) ){
					if(i > 1){ 
						Y = Y + H/nrow(H)*(1-row.base)
					}
					xind = seq( hlcat[i], nc*ntc0, ntc0 )
					H = as.matrix(H0[,xind])
					if( nr < 2 ){
						H = t(H)	
					}
					draw(H/nrow(H)*(1-row.base),W3/ncol(W3)*(1-col.base),X3,Y,alpha = 1, bg = col.Mat[,xind])
				}
				Y = Y0
			}else{
				H = apply(H,c(1,2),function(x) return(min(x/max.scale,1)))
				draw(H/nrow(H)*(1-row.base),W3/ncol(W3)*(1-col.base),X3,Y,alpha = 1, bg = col.Mat)
			}
			draw(H2/nrow(H2)*(1-row.base),W2sqrt/ncol(W2sqrt)*(1-col.base) + (round(W2sqrt,digits=7) > 0)*width.cor,X2,Y2,alpha = 0.25, bg = "grey")
			draw(H2/nrow(H2)*(1-row.base),H2/ncol(H2)*(1-col.base) + width.cor,X2,Y2,alpha = 0.1, bg = "grey")
		
		}
	}
		
	if(yaxis){
		sapply(y,function(z){
			grid.yaxis(at =seq(z,z+1/length(y)*(1-row.base),1/length(y)*(1-row.base)/5), label = round(seq(0,max.scale,max.scale/5),digits = 3), main = TRUE,
			edits = NULL, name = NULL,
			gp = gpar(), draw = TRUE, vp = NULL)
		})
		sapply(y,function(z){
			grid.yaxis(at =seq(z,z+1/length(y)*(1-row.base),1/length(y)*(1-row.base)/5), label = round(seq(0,max.scale,max.scale/5),digits = 3), main = FALSE,
			edits = NULL, name = NULL,
			gp = gpar(), draw = TRUE, vp = NULL)
		})
	}
	upViewport()
# ----- ---------------- ----- #	
# ----- labeling section ----- #	
# ----- ---------------- ----- #		
	

	
	
	if(lab){
			if(any(hsplit*label)){
	# ----- labeling the x-axis ----- #	

	vpX <- viewport(x = 1-yaxis*s1 - s3*(!is.null(expected)), y = 1-ncl*s2, w = 1-nrl*s2 - yaxis*2*s1 - s3*(!is.null(expected)), h = ncl*s2, just = c("right", "bottom"), name = "vpX")
	pushViewport(vpX)
	cur = 1
	zc = cumprod(col.nlvl)[ which(label[which(hsplit)]) ]
    for(i in 1:ncl){
		vpXX <- viewport(x = 1, y = 1 - i*1/ncl, w = 1, h = 1/ncl, just = c("right", "bottom"), name = "vpXX")
		pushViewport(vpXX)
			if( varnames ){
				grid.text(  names(dset)[cind[i]],x = 0.5 , y = 5/6, just = "centre",gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
			}
			#z = (nlvl[cind])[i]*cur
			z = zc[i]
			grid.text(  rep( clabs[[i]] ,cur), x = X[1,seq(1,nc,nc/z)] + (X[1,seq(nc/z,nc,nc/z)] - X[1,seq(1,nc,nc/z)] +1/nc*(1-col.base))/2, y = 1/3, just = "centre",gp=gpar(cex=lab.cex)) 
			if( boxes ){
				grid.rect(  y = 1/3 , x = X[1,seq(1,nc,nc/z)] , just = c("left","centre"), width = X[1,seq(nc/z,nc,nc/z)] - X[1,seq(1,nc,nc/z)] +1/nc*(1-col.base), height = 1/2, gp = gpar(fill="transparent"))
			}
			cur = z
		popViewport()
	}
		
	upViewport()
	}
		if(any( (!hsplit)*label)){
# ----- labeling the y-axis ----- #			
	vpY <- viewport(x = nrl*s2, y = 0, w = nrl*s2, h = 1-ncl*s2, just = c("right", "bottom"), name = "vpY")
	pushViewport(vpY)
	cur = 1
	zr = cumprod(row.nlvl)[ which(label[which(!hsplit)]) ]
	for(i in 1:nrl){
		
		vpYY <- viewport(x = 0 + (i-1)*1/nrl, y = 0, w = 1/nrl, h = 1, just = c("left", "bottom"), name = "vpYY")
		pushViewport(vpYY)
			if( varnames & (nr > 1) ){
				grid.text(  names(dset)[rind[i]],x = 1/6 , y = 0.5, just = "centre",rot=90,gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
			}			
			#z = (nlvl[rind])[i]*cur
			z = zr[i]
			grid.text(  rep( rlabs[[i]] ,cur), y = Y[seq(1,nr,nr/z),1] + (Y[seq(nr/z,nr,nr/z),1] - Y[seq(1,nr,nr/z),1] +1/nr*(1-row.base))/2, x = 2/3, just = "centre",rot=90,gp=gpar(cex=lab.cex)) 
			if( boxes ){
				grid.rect(  x = 2/3 , y = Y[seq(1,nr,nr/z),1] , just = c("centre","bottom"), height = Y[seq(nr/z,nr,nr/z),1] - Y[seq(1,nr,nr/z),1] +1/nr*(1-row.base), width = 1/2, gp = gpar(fill="transparent"))
			}
					
			cur = z
		popViewport()
	}
	upViewport()
	}
	}
	# ----- scale for expected mode ----- #
	if(!is.null(expected)){

		vpS <- viewport(x = 1, y = 0 , w = s3, h = 0.8, just = c("right","bottom"), name = "vpS")
		pushViewport(vpS)
		colv = rgb(rep(c(1,0),each=cut.rs),0,rep(c(0,1),each=cut.rs),alpha = c(seq(1,1/cut.rs,-1/cut.rs),seq(1/cut.rs,1,1/cut.rs))*base.alpha)
	
		grid.rect(  x = 0.1 , y = seq(0,1-1/2/cut.rs,1/2/cut.rs) , just = c("left","bottom"), height = 1/2/cut.rs, width = 0.3 , gp = gpar(fill=colv))
		grid.text( label=round(seq(-resid.max,resid.max,resid.max/min(10,cut.rs)),digits=2), y = seq(0,1,1/2/min(10,cut.rs)), x = 0.45, just = c("left","centre"),rot=0,gp=gpar(cex=lab.cex))
		upViewport()	
		vpS2 <- viewport(x = 1, y = 0.8  , w = 0.1, h = 0.2 , just = c("right","bottom"), name = "vpS2")
		pushViewport(vpS2)
		
		grid.text( label=resid.type, y = 0.5, x = 0.2, just = "centre",rot=90,gp=gpar(cex=lab.cex*1.2))
		grid.text( label="residuals", y = 0.5, x = 0.4, just = "centre",rot=90,gp=gpar(cex=lab.cex*1.2))
		grid.text( label=paste("p.value =",format(p.value,digits=2)), y = 0.5, x = 0.7, just = "centre",rot=90,gp=gpar(cex=lab.cex*1.1))
	}
}

# ------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------ #


space = function(nlvl,hsplit,mult = 1.1,base = 0.1,last.col.zero = F, last.row.zero = F){
	# 	>>> function to compute the gaps widths for the rmb plot
	# 	>>> might be extended in further versions of the package
	
	# ----- parameter check ----- #
		if( !(  (base > 0)&(base < 1) ) ){
			stop(simpleError("Wrong base specification!"))
		}
	
	col.nlvl = nlvl[hsplit]
	row.nlvl = nlvl[!hsplit]								

	nc = prod(col.nlvl)
	nr = prod(row.nlvl)

	ncv = length(col.nlvl)
	nrv = length(row.nlvl)

	col.spaces = rep(base,nc-1)
	row.spaces = rep(base,nr-1)

	if( ncv > 1 ){
		cind = 1:(nc-1)
		cur = 1
		for( i in (ncv-1):1 ){ 
			cur = cur*col.nlvl[i+1]
			curind = which( (cind %% cur) == 0)
			col.spaces[ curind  ] = col.spaces[ curind  ] * mult
		}
	}
	if( nrv > 1 ){
		rind = 1:(nr-1)
		cur=1
		for( i in (nrv-1):1 ){ 
			cur = cur*row.nlvl[i+1]
			curind = which( (rind %% cur) == 0)
			row.spaces[ curind  ] = row.spaces[ curind  ] * mult
		}
	}
	if(last.col.zero){
		col.spaces[which(col.spaces == base)] = 0
		col.spaces = col.spaces/base
	}
	if(last.row.zero){
		row.spaces[which(row.spaces == base)] = 0
		row.spaces = row.spaces/base
	}
	return( list(col.spaces,row.spaces))
}

# ------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------ #


spread = function(M,ncol = 1,nrow = 1){
	# 	>>> this function turns each matrix entry into a nrow x ncol - matrix
	M = as.matrix(M)
	M2 = apply(M,2,function(x) rep(x,each = nrow)  )
	if( (nrow(M) == 1)&(nrow == 1) ){
		M2 = matrix(M2,ncol = ncol(M))
	}
	M3 = t(apply(t(M2),2,function(x) rep(x,each = ncol)  ))
	if( (ncol(M) == 1)&(ncol == 1) ){
		M3 = matrix(M3,nrow = nrow(M2))
	}
	return(M3)
}

