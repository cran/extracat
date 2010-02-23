
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


rmbX = function(env, init = T, f,dset,spine = F,hlcat = 1,  hsplit = NULL,expected = NULL,resid.max = NULL,
eqwidth=F,tfreq=1,max.scale = 1,...){

	# ----- read global parameters from e1 ----- #
	mod.type = env$mod.type
	resid.type = env$resid.type
	base = env$base
	mult = env$mult
	varnames = env$varnames
	min.alpha = env$min.alpha
	abbr = env$abbr
	lab.cex = env$lab.cex
	lab.tv = env$lab.tv
	yaxis = env$yaxis
	boxes = env$boxes
	dset = env$dset
	colv = env$colv
	base.alpha = env$base.alpha
	use.expected.values = env$use.expected.values
	cut.rv = env$cut.rv
	cut.rs = env$cut.rs
	lab = env$lab
	# NOT INCLUDE: lab, eqwidth
	
	# ----- read special parameters from e1 ----- #
	resid.max = env$resid.max
	max.scale = env$max.scale
if(!init){
	X = env$X
	X2 = env$X2
	Y = env$Y
	y = env$y
	Y2 = env$Y2
	H = env$H
	H2 = env$H2
	W = env$W
	W2 = env$W2
	S = env$S
	tfreq = env$tfreq

	H0 = H
	Y0 = Y

	eqwidth = env$eqwidth
	spine = env$spine
	hsplit = env$hsplit
	expected = env$expected
	resid.max = env$resid.max
	lab = env$lab
	if(env$spinechange){
		max.scale = 1
	}else{
		max.scale = env$max.scale
	}
	p.value = env$p.value
	hlcat = env$hlcat
}
if(!suppressWarnings(any(hlcat))){
	hlcat = 1	
}
# ----- parameter check 1 ----- #
		if(!is.formula(f) ){
			stop(simpleError("Wrong formula specification!"))
		}
		if(!is.data.frame(dset) ){
			stop(simpleError("Wrong dset specification!"))
		}
		if(!( (base < 1) & (base > 0) ) ){
			stop(simpleError("Wrong base specification!"))
		}
		if( !((min.alpha <= 1) & (min.alpha >= 0)) ){
			stop(simpleError("Wrong alpha specification!"))
		}
	

# ----- terms extraction ----- #
	fterms = attr(terms(f),"term.labels")
	nv = length(fterms)
	if(is.null(hsplit)){
		hsplit = rep(c(T,F),nv)[1:nv]
		env$hsplit = hsplit	
	}
	tv = fterms[nv]
	hsplit[nv] = T
	
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
	label = rep(T,nv) # if is null label....
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

if(init | env$spinechange){
# ----- computation of the underlying (relative) frequencies ----- #		
	tt1 = ftable(tapply(dset$Freq,as.list(dset[,1:nv]),sum),col.vars=which(hsplit))
	tt2 = spread(ftable(tapply(dset$Freq,as.list(dset[,1:(nv-1)]),sum),col.vars=which(hsplit[1:(nv-1)])),ncol=ntc)

	tt1[is.na(tt1)] = 0 
	tt2[is.na(tt2)] = 1 

	H0 = tt1/tt2
	H = H0
}
# ----- a few more auxiliary variables ----- #	
	nlvl = sapply(dset[,(sapply(dset,class)=="factor")], nlevels) # muss nicht sein
	
	col.nlvl = nlvl[hsplit]
	row.nlvl = nlvl[!hsplit]								

	nc = prod(col.nlvl)
	nr = prod(row.nlvl)

	row.base = ifelse(nr > 1, base * min(1,nrv/ncv), 0) 
	col.base = ifelse(nc > ntc0, base * min(1,ncv/nrv), 0) 

	ncv = length(col.nlvl)
	nrv = length(row.nlvl)


	rind = which( (!hsplit)*label > 0)
	cind = which(hsplit*label > 0)
	nrl = length(rind)*lab
	ncl = length(cind)*lab
# ----- model computation in expected mode ----- #	

	if(!is.null(expected) & (init|env$modchange)){
			
		int.dset = lapply(dset,as.integer)
		dset = dset[do.call("order",int.dset),]
		if(!(mod.type %in% c("poisson","polr"))){
			stop(simpleError("Wrong mod.type specification!"))
		}
		# ----- standardized or deviance residuals in poisson model ----- #	
			nameslist = names(dset)[1:(nv)]		
			nameslist = nameslist[which(  !(nameslist %in% c("probability","distribution"))  )]
			single.terms = do.call("paste",c(as.list( nameslist ),sep="+"))
			if(init){
				interaction.terms = do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(nameslist[x]),sep="*"))),sep="+"))
			}else{
				interaction.terms = env$interaction.terms
			}
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
		

# ----- modify parameters for spineplot mode ----- #	
	if( spine ){
		nlvl[nv] = 1
		nc = nc/ntc
		ntc = 1
	}

if(init | env$spinechange){
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
}else{
	tt3 = W2 *sum(dset$Freq)/sum(W2)
}
	if(nc > nlevels(dset[,nv])){
		width.cor = (ntc-1)*S[[1]][1]/sum(S[[1]])*col.base
	}else{
		width.cor = 0
	}
# ----- ---------------- ----- #	
# ----- plotting section ----- #	
# ----- ---------------- ----- #

# ----- opening the basic viewport ----- #	
	s0 = 0.02 # border
	s1 = 0.04 # yaxis
	s2 = 0.06 # labs
	s3 = 0.1  # expected scale
grid.rect()
	vp0 <- viewport(x = 0.5, y = 0.5, w = 1- 2*s0 , h = 1 - 2*s0, just = "centre", name = "vp0")
	pushViewport(vp0)
	vp1 <- viewport(x = 1-yaxis*s1 - s3*(!is.null(expected)) , y = 0, w = 1- nrl*s2 - yaxis*2*s1 - s3*(!is.null(expected)), h = 1-ncl*s2, just = c("right", "bottom"), name = "vp1")
	pushViewport(vp1)
	
	alpha = spread(W2,ncol=ntc)
	alpha[alpha < min.alpha] = min.alpha
	if(!is.null(expected)){
		alpha[alpha < 1] = 1	
	}

# ----- preparing the color matrix ----- #	
if(init | env$modchange){
			if(is.null(expected)){

				if(is.null(env$colv) | length(env$colv) < ntc0){	
					#if( (!spine) & eqwidth ){
					#	env$colv = rep(rainbow(1,alpha = base.alpha),ntc0)
					#}else{
					#	env$colv = rainbow(ntc0,alpha = base.alpha)
					#}
					env$colv = rainbow(ntc0,alpha=base.alpha)
					if((!spine) & eqwidth){
						env$colv = rep(rainbow(1,alpha=base.alpha),ntc)
					}
				}				
				colv = env$colv[1:ntc0]
				col.Mat = spread( t(rep(colv , nc/ntc)),nrow=nr)
			}else{
				resid.mat[is.nan(resid.mat)] = 0
				resid.mat[is.na(resid.mat)] = 0
				if(is.null(env$max.resid.scale)){
					resid.max = max(abs(resid.mat),na.rm=T)	
				}	
			
				resid.mat[resid.mat < -resid.max] = -resid.max
				resid.mat[resid.mat > resid.max]  =  resid.max
				rmax = ifelse(resid.max == 0, 1, resid.max)
				#col.Mat = apply(resid.mat/rmax,2,function(x) rgb(sign(x) < 0,0,sign(x) > 0,alpha = abs(x)  ))	
				col.Mat = apply(resid.mat/rmax,2,function(x){
						qx = abs(x)
						if(cut.rv){ qx =  ceiling(qx * cut.rs)/cut.rs }
					rgb(sign(x) < 0,0,sign(x) > 0,alpha = qx*base.alpha )
					})
			}
			env$col.Mat = col.Mat
}else{
			col.Mat = env$col.Mat	
}
# ----- plotting itself ----- #
# 	>>> divided into eqwidth = T | eqwidth = F >>> tfreq = Id/sqrt/log
# 	>>> in spineplot mode, the target categories will be plotted one after another (see the for-loop)
# 	>>> the 3 diff. draw-function add the rectangles for the relative frequencies, the weights and the background respectively
	if(eqwidth){
		if(spine){
			W2 = apply(W2,c(1,2),function(x) return(min(x/tfreq,1)))
		}
			draw(H2/nrow(H2)*(1-row.base),W2/ncol(W2)*(1-col.base) + (round(W2,digits=7) > 0)*width.cor,X2,Y2,alpha = 0.25, bg = "grey")
			draw(H2/nrow(H2)*(1-row.base),H2/ncol(H2)*(1-col.base) + width.cor,X2,Y2,alpha = 0.1, bg = "grey")
	
		if( spine ){
			for( i in 1:length(hlcat) ){
				if(i > 1){ 
					Y = Y + H/nrow(H)*(1-row.base)
				}
				xind = seq( hlcat[i], nc*ntc0, ntc0 )
				H = H0[,xind]
			 	draw(H/nrow(H)*(1-row.base)*max.scale,W/ncol(W)*(1-col.base),X,Y,alpha = ifelse(is.null(expected),0.9,1), bg = col.Mat[,xind])
			}
			H = H0
			Y = Y0
		}else{
			H = apply(H,c(1,2),function(x) return(min(x/max.scale,1)))
			draw(H/nrow(H)*(1-row.base),W/ncol(W)*(1-col.base),X,Y,alpha = alpha, bg = col.Mat)
		}
	}else{
		if(is.numeric(tfreq)){
			W2num = apply(tt3,c(1,2),function(x) return(min(x/(max(tt3)*tfreq),1)))
			W3 = spread(W2num,ncol=ntc)
			if( nr > 1 ){
				X3 = spread(X[,seq(1,nc,ntc)],ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.base) 
			}else{
				X3 = spread(t(X[,seq(1,nc,ntc)]),ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.base) 
			}
			draw(H2/nrow(H2)*(1-row.base),W2num/ncol(W2num)*(1-col.base) + (round(W2num,digits=7) > 0)*width.cor,X2,Y2,alpha = 0.25, bg = "grey")
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
				H = H0
				Y = Y0
			}else{
				H = apply(H,c(1,2),function(x) return(min(x/max.scale,1)))
				draw(H/nrow(H)*(1-row.base),W3/ncol(W3)*(1-col.base),X3,Y,alpha = 1, bg = col.Mat)
			}
		}
	}
	if(yaxis & max.scale > 0){
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
		if(any(!hsplit * label)){
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
		upViewport(1)
	}
	upViewport(1)
if(init | env$spinechange){	
	env$X =X
	env$X2 = X2
	env$Y = Y
	env$y = y
	env$Y2 = Y2
	env$H = H
	env$H2 = H2
	env$W = W
	env$W2 = W2
	env$S = S
	env$tfreq = tfreq	
	env$eqwidth = eqwidth
	env$spine = spine
	env$hlcat = hlcat
	env$expected = expected
	env$resid.max = resid.max
	env$max.scale = max.scale
	env$lab = lab
	env$fterms = fterms
	if(!is.null(expected)){
		env$p.value = p.value
	}
}
}

 irmb = function(f, dset, Z=100, use.na = FALSE, expected = TRUE,
               resid.type = "pearson", use.expected.values = FALSE,
               max.resid.scale = NULL, cut.rv = TRUE, cut.rs = 5, base = 0.2,
               mult = 1.5, colv = NULL, yaxis = TRUE, min.alpha = 0.1,
               base.alpha = 0.75, boxes = TRUE, lab.tv = FALSE, varnames = TRUE,
               abbr = FALSE, lab.cex = 1.2,...){

	if(!("iWidgets" %in% .packages(all.available=TRUE) ) ){
		cat("The required package iWidgets is not available.\n Please use \n install.packages('iWidgets',,'http://rforge.net/',type='source') \n for installation.")
		return(invisible(TRUE))
		}
	if(!("JGR" %in% .packages(all.available=FALSE) ) ){
		cat("iWidgets requires 'JGR' in order to run correctly.\n Please visit \n http://www.rosuda.org/software/ \n to download it.")
		return(invisible(TRUE))
		}
	library(iWidgets)
	
	e1 = new.env()
	w = iwindow()
	# ----- save global parameters to e1 ----- #
	e1$mod.type = "poisson"
	e1$resid.type = resid.type
	e1$use.expected.values =  use.expected.values
	e1$cut.rv = cut.rv
	e1$cut.rs = cut.rs
	e1$base = base
	e1$mult = mult
	e1$varnames = varnames
	e1$min.alpha = min.alpha
	e1$abbr = abbr
	
	e1$lab.cex = lab.cex
	e1$lab.tv = lab.tv
	e1$yaxis = yaxis
	e1$boxes = boxes
	e1$colv = colv
	e1$base.alpha = base.alpha
	e1$fterms = attr(terms(f),"term.labels")
	e1$dset = subtable(dset,cols = which(names(dset) %in% e1$fterms),keep.zero=F)
	
	for( i in 1:(ncol(e1$dset)-1) ){
		ind = is.na(e1$dset[,i])
		if(any(ind)){
			levels(e1$dset[,i]) = c(levels(e1$dset[,i]),"N/A")	
			e1$dset[which(ind),i] = "N/A"
		}	
	}
	# NOT INCLUDE: lab, eqwidth
	
	# ----- save important initial parameter settings to e1 ----- #
	e1$max.scale = 1
	e1$max.resid.scale = max.resid.scale
	e1$resid.max = max.resid.scale
	e1$tfreq = 1
	e1$f = f
	
	if(expected){
		e1$expected = as.list(1:length(attr(terms(f),"term.labels")))
	}else{
		e1$expected = NULL	
	}
	e1$eqwidth = FALSE
	e1$spine = F
	e1$spinechange = F
	e1$hlcat = 1
	e1$updateplot = T
	e1$modchange = T
	e1$hsplit = NULL
	e1$lab = TRUE
	g = igroup()
	g2 = igroup(window=g,horizontal = F)
	add(g2,igraphics(width=1024,height=768))
	rmbX(env = e1,init =T,f,dset = e1$dset,spine = F,hlcat = 1,expected = e1$expected,  hsplit = e1$hsplit,max.scale = 1,eqwidth=e1$eqwidth,lab = e1$lab,tfreq = 1)
	mslider = islider(min=0,max=Z,horizontal = F,value = Z, window=g, handler=function(h,...){
		e1$max.scale = get.value(h$obj)/Z
		if((!e1$spine) & e1$updateplot){
			rmbX(env = e1 ,init = F, f = e1$f,dset = e1$dset)
		}
		})
	tslider = islider(min=0,max=Z,horizontal = T,value = Z, window=g2, handler=function(h,...){
		e1$tfreq = get.value(h$obj)/Z
		if(e1$updateplot){
			rmbX(env = e1 ,init = F, f = e1$f,dset = e1$dset)
		}
		})
	g3 = igroup(window=g2,horizontal = F)
	g3a = igroup(window=g3,horizontal = T)
	icheckbox("spineplot", checked = F, window = g3a,handler=function(h,...){
			e1$spine = !e1$spine
			e1$spinechange = T
			ntc = nlevels(e1$dset[,which(names(e1$dset) == e1$fterms[length(e1$fterms)])])
			if(max(e1$hlcat) > ntc){
				e1$hlcat = 1:ntc
				set.value(hledit,do.call(paste,c(as.list(1:ntc),sep=",")))
			}
			rmbX(env = e1 ,init = F, f = e1$f,dset = e1$dset)
			e1$spinechange = F
			})
	icheckbox("eqwidth", checked = e1$eqwidth, window = g3a,handler=function(h,...){
			e1$eqwidth = !e1$eqwidth
			
			rmbX(env = e1 ,init = F, f = e1$f,dset = e1$dset)
			
			})
	icheckbox("lab", checked = e1$lab, window = g3a,handler=function(h,...){
			e1$lab = !e1$lab
			
			rmbX(env = e1 ,init = F, f = e1$f,dset = e1$dset)
			
			})
	g4 = igroup(window=g3,horizontal=T)
	add.space(g4,32)
	ilabel("highlighting categories:",window=g4)
	add.space(g4,10)
	hledit = iedit(window=g4,text = e1$hlcat, handler=function(h,...){
		tmp = e1$hlcat
		op=options()
		e1$hlcat <- lapply(strsplit(get.value(h$obj),"\\,"),as.integer  )[[1]]
		options(error = function(){
				e1$hlcat = tmp
				set.value(h$obj,do.call(paste,c(as.list(tmp),sep=",")))
			})
			
		if(e1$spine & e1$updateplot){
			ntc = nlevels(e1$dset[,which(names(e1$dset) == e1$fterms[length(e1$fterms)])])
			if(max(e1$hlcat) > ntc){
				e1$hlcat = 1:ntc
				set.value(hledit,do.call(paste,c(as.list(1:ntc),sep=",")))
			}
			rmbX(env = e1 ,init = F, f = e1$f,dset = e1$dset)
			options(op)
			}
		})
		add.space(g4,20)
	if(expected){
		g5 = igroup(window=g3,hroizontal=T)
		add.space(g5,10)
		e1$interaction.terms = do.call("paste",c(lapply(e1$expected,function(x) do.call("paste",c(as.list(names(e1$dset)[x]),sep="*"))),sep="+"))
		modlab = ilabel(text = "expected interaction terms:",window = g5)
		add.space(g5,10)
		modedit = iedit(window = g5,text =  e1$interaction.terms, width = 40, handler = function(h,...){
				tmp = e1$interaction.terms
				op = options()
				e1$interaction.terms = get.value(h$obj)
				if(e1$updateplot){
					options(error = function(){
						e1$interaction.terms = tmp
						set.value(h$obj,tmp)
						})
						e1$modchange = T
						rmbX(env = e1 ,init = F, f = e1$f,dset = e1$dset)
						options(op)
						e1$modchange = F
				}
			})
			add.space(g5,20)
	}

	g5 = igroup(window = g2, horizontal = T)
	e1$nv = length(e1$fterms)
	e1$genhsplit = e1$hsplit
	e1$genfterms = e1$fterms
	add.space(g5,97)
	ilabel(text = "variables:",window = g5)
	add.space(g5,10)
	for(i in 1:e1$nv){
			icheckbox(e1$genfterms[i], checked = T, window = g5,handler=function(h,...){
			if(get.value(h$obj)){
				e1$fterms = c(e1$fterms,h$obj$text)
				e1$hsplit = c(e1$hsplit,e1$genhsplit[which(e1$genfterms == h$obj$text)])
				e1$nv = e1$nv+1
			}else{
				e1$fterms = e1$fterms[which(e1$fterms != h$obj$text)]
				e1$hsplit = e1$hsplit[which(e1$fterms != h$obj$text)]
				e1$nv = e1$nv-1
				set.value(hcheck[[which(e1$genfterms == e1$fterms[length(e1$fterms)])]],TRUE)
			}
			if(length(e1$fterms) > 0){
				tmp = paste("~",do.call(paste,c(as.list(e1$fterms),sep="+")),sep="")
				e1$f = as.formula(tmp)
				set.value(fx, paste("      Current formula: >> ",tmp," <<   "))
				e1$updateplot = F
				set.value(hledit,do.call(paste,c(as.list(1:nlevels(e1$dset[,which(names(e1$dset) == e1$fterms[e1$nv])])),sep=",")))
				set.value(tslider,Z)
				set.value(mslider,Z)
				e1$updateplot = T
				
				if(expected){
					e1$expected = as.list(1:length(e1$fterms))
					set.value(modedit,do.call("paste",c(as.list(e1$fterms),sep="+")))
					e1$modchange = T
				}
				set.value(hcheck[[which(e1$genfterms == e1$fterms[length(e1$fterms)])]],TRUE)
				rmbX(env = e1,init =T,e1$f,dset = e1$dset,spine = e1$spine,hlcat = e1$hlcat,expected = e1$expected,  hsplit = e1$hsplit,eqwidth = e1$eqwidth,lab = e1$lab,tfreq = e1$tfreq,max.scale = e1$max.scale)
				e1$modchange = F
			}
			})
	}
	add.space(g5,50)
	fx = ilabel(text = do.call(paste,as.list(as.character(e1$f))),window = g5)
	add.spring(g5)
	g6 = igroup(window = g2, horizontal = T)
	add.space(g6,116)
	ilabel(text = "splits:",window = g6)
	add.space(g6,10)
	hcheck = vector(length = e1$nv, mode = "list")
	for(i in 1:e1$nv){
			hcheck[[i]] = icheckbox(e1$genfterms[i], checked= e1$genhsplit[i], window = g6,handler=function(h,...){
				e1$hsplit[which(e1$fterms == h$obj$text)] = get.value(h$obj)
				e1$genhsplit[which(e1$genfterms == h$obj$text)] = get.value(h$obj)
				e1$updateplot = F
				set.value(tslider,Z)
				set.value(mslider,Z)
				e1$updateplot = T
				rmbX(env = e1,init =T,e1$f,dset = e1$dset,spine = e1$spine,hlcat = e1$hlcat,expected = e1$expected,  hsplit = e1$hsplit,eqwidth = e1$eqwidth,lab = e1$lab,tfreq = e1$tfreq,max.scale = e1$max.scale)
			})
	}	
	add.spring(g6)
	
	add(w,g)
	visible(w,TRUE)
}

# ------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------ #


space = function(nlvl,vsplit,mult = 1.1,base = 0.1,last.col.zero = F, last.row.zero = F){
	# 	>>> function to compute the gaps widths for the rmb plot
	# 	>>> might be extended in further versions of the package
	
	# ----- parameter check ----- #
		if( !(  (base > 0)&(base < 1) ) ){
			stop(simpleError("Wrong base specification!"))
		}
	
	col.nlvl = nlvl[vsplit]
	row.nlvl = nlvl[!vsplit]								

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




