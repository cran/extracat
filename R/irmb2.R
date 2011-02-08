
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


rmbX = function(env, init = TRUE,...){

	# ----- read global parameters from e1 ----- #
	mod.type = env$mod.type
	resid.type = env$resid.type
	gap.prop = env$gap.prop
	gap.mult = env$gap.mult
	varnames = env$varnames
	min.alpha = env$min.alpha
	abbrev = env$abbrev
	lab.cex = env$lab.cex
	lab.tv = env$lab.tv
	yaxis = env$yaxis
	boxes = env$boxes
	#dset = env$dset
	colv = env$colv
	col = env$col
	col.opt = env$col.opt
	base.alpha = env$base.alpha
	use.expected.values = env$use.expected.values
	cut.rv = TRUE
	cut.rs = env$cut.rs
	lab = env$lab
	# NOT INCLUDE: lab, eqwidth
	
	# ----- read special parameters from e1 ----- #
	resid.max = env$resid.max
	max.scale = env$max.scale
	
#if(!init){

	freq.trans = env$freq.trans

dset.changed = FALSE
	eqwidth = env$eqwidth
	spine = env$spine
	col.vars = env$col.vars
	expected = env$expected
	resid.max = env$resid.max
	lab = env$lab
	if(env$spinechange){
		max.scale = 1
	}else{
		max.scale = env$max.scale
	}
	p.value = env$p.value
	cat.ord = env$cat.ord
	col = env$col
#}

# ----- parameter check 1 ----- #

		if(!is.data.frame(env$dset) ){
			stop(simpleError("Wrong dset specification!"))
		}
		if(!( (gap.prop < 1) & (gap.prop > 0) ) ){
			stop(simpleError("Wrong gap.prop specification!"))
		}
		if( !((min.alpha <= 1) & (min.alpha >= 0)) ){
			stop(simpleError("Wrong alpha specification!"))
		}
	

# ----- terms extraction ----- #
	fterms = env$fterms#attr(terms(f),"term.labels")
	nv = length(fterms)
	tv = fterms[nv]

	hsplit = env$col.vars[ match(fterms, env$genfterms)  ]

	ind =  match(fterms,names(env$dset))
# 	>>> a dummy-variable is used to handle the case of no row-variables
	if( sum(!hsplit) < 1 ){
		env$dset = data.frame(probability = rep("prob",nrow(env$dset)),env$dset)
		nv=nv+1
		hsplit = c(F,hsplit)
		ind = c(1,ind+1)
		dset.changed = TRUE
	}
	if( sum(hsplit) < 2 ){
		env$dset = data.frame(distribution = rep(tv,nrow(env$dset)),env$dset)
		nv=nv+1
		hsplit = c(T,hsplit)
		ind = c(1,ind+1)
		if("probability" %in% names(env$dset)){ expected = NULL }
		dset.changed = TRUE
	}
# 	>>> preparing the labels for plotting including the abbreviation option	
	

		
				
	
	if(is.numeric(abbrev)){
			abbr = TRUE
	}else{
			abbr = abbrev	
	}
	
	if(!abbr){
		rclabs = lapply(env$dset[,ind],function(x) levels(as.factor(x)))
	}else{
		rclabs = lapply(env$dset[,ind],function(x) abbreviate(levels(as.factor(x)),abbrev))
	}
	
	lab.tv = ifelse(spine,FALSE,lab.tv)
	label = rep(T,nv) # if is null label....
	label[nv] = lab.tv
	if( sum(label) == 0 ){
		lab = FALSE	
	}
	
	rlabs = lapply(which( (!hsplit)*label > 0),function(x) rclabs[[x]])
	clabs = lapply(which(hsplit*label > 0),function(x) rclabs[[x]])

	env$dset = subtable(env$dset,ind,allfactor=TRUE)

	if(length(cat.ord) > 1 & !spine){
		#levels(env$dset[,nv])[which(! 1:nv %in% cat.ord )] = NA
		env$dset[,nv] = factor(env$dset[,nv],levels=levels(env$dset[,nv])[cat.ord])
		#env$dset[,nv] = factor(env$dset[,nv],levels=levels(env$dset[,nv])[rank(cat.ord)])
		if(lab.tv){
			clabs[[length(clabs)]] = levels(env$dset[,nv])	
		}
		env$dset = na.omit(env$dset)
		dset.changed = TRUE
	}
	
# ----- descriptive parameters ----- #	
	ntc = nlevels(env$dset[,nv])
	ntc0 = ntc
	orig.ntc = env$ntc
	#nlvl = sapply(env$dset[,(sapply(env$dset,class)=="factor")], nlevels)
	nlvl = sapply(env$dset[,-ncol(env$dset)], nlevels)

	col.nlvl = nlvl[hsplit]
	row.nlvl = nlvl[!hsplit]								

	nc = prod(col.nlvl)
	nr = prod(row.nlvl)

	ncv = length(col.nlvl)
	nrv = length(row.nlvl)

if(init | env$spinechange){
# ----- computation of the underlying (relative) frequencies ----- #		
	tt1 = ftable(tapply(env$dset$Freq,as.list(env$dset[,1:nv]),sum),col.vars=which(hsplit))
	tt2 = spread(ftable(tapply(env$dset$Freq,as.list(env$dset[,1:(nv-1)]),sum),col.vars=which(hsplit[1:(nv-1)])),ncol=ntc)

	tt1[is.na(tt1)] = 0 
	tt2[is.na(tt2)] = 1 

	env$H0 = tt1/tt2
	env$H = env$H0
}

# ----- a few more auxiliary variables ----- #	
	#nlvl = sapply(env$dset[,(sapply(env$dset,class)=="factor")], nlevels) # muss nicht sein
	nlvl = sapply(env$dset[,-ncol(env$dset)], nlevels)
	
	col.nlvl = nlvl[hsplit]
	row.nlvl = nlvl[!hsplit]								

	nc = prod(col.nlvl)
	nr = prod(row.nlvl)

	row.gap.prop = ifelse(nr > 1, gap.prop * min(1,nrv/ncv), 0) 
	col.gap.prop = ifelse(nc > ntc0, gap.prop * min(1,ncv/nrv), 0) 

	ncv = length(col.nlvl)
	nrv = length(row.nlvl)


	rind = which( (!hsplit)*label > 0)
	cind = which(hsplit*label > 0)
	nrl = length(rind)*lab
	ncl = length(cind)*lab
# ----- model computation in expected mode ----- #	

	if(!is.null(expected) & (init|env$modchange)){
			
		int.dset = lapply(env$dset,as.integer)
		env$dset = env$dset[do.call("order",int.dset),]
		if(!(mod.type %in% c("poisson","polr"))){
			stop(simpleError("Wrong mod.type specification!"))
		}
		# ----- standardized or deviance residuals in poisson model ----- #	
			nameslist = names(env$dset)[1:(nv)]		
			nameslist = nameslist[which(  !(nameslist %in% c("probability","distribution"))  )]
			single.terms = paste( nameslist ,collapse="+")
			#if(init){
			#	interaction.terms = do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(nameslist[x]),sep="*"))),sep="+"))
			#}else{
				interaction.terms = env$interaction.terms
			#}
			
			mod.formula = paste("Freq~",single.terms,"+",interaction.terms,sep="")
			mod = glm( as.formula(mod.formula),data = env$dset,  family = "poisson")
		    mod.residuals = residuals(mod,type=resid.type)
			resid.mat =  ftable(tapply(mod.residuals,as.list(env$dset[,1:nv]),sum),col.vars=which(hsplit))
			
			Df  = mod$df.residual
			Dev = mod$deviance
			p.value = 1-pchisq(Dev,Df)
			
			if( use.expected.values ){
				pred = predict(mod,type="response")
				pred.mat =  ftable(tapply(pred,as.list(env$dset[,1:nv]),sum),col.vars=which(hsplit))
				env$H0 = pred.mat/tt2
				env$H = env$H0
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
	env$W = matrix(1,ncol=nc,nrow=nr)
	env$S = space(nlvl,hsplit,gap.prop=gap.prop,gap.mult=gap.mult,last.col.zero = T,last.row.zero = F)
	env$x = ( seq(0,1-1/nc,1/nc)   )*(1-col.gap.prop)
	if(nc > nlevels(env$dset[,nv])){
		env$x = env$x + c(0,cumsum(env$S[[1]]))*col.gap.prop/sum(env$S[[1]])
	}
	if( suppressWarnings(any(env$S[[2]])) ){
		env$y = rev((seq(0,1-1/nr,1/nr))*(1-row.gap.prop) + c(0,cumsum(env$S[[2]]))*row.gap.prop/sum(env$S[[2]]))# rev test
	}else{
		env$y = 0	
	}
		
	env$X = spread(t(matrix( env$x )),nrow=nr)
	env$Y0 = spread(matrix( env$y ),ncol=nc)
	env$Y = env$Y0
	env$X2 = env$X[,seq(1,nc,ntc)]
	env$Y2 = env$Y[,seq(1,nc,ntc)]

	env$H2 = matrix(1,ncol=nc/ntc,nrow=nr)
	tt3 = ftable(tapply(env$dset$Freq,as.list(env$dset[,1:(nv-1)]),sum),col.vars=which(hsplit[1:(nv-1)]))
	
	tt3[is.na(tt3)] = 0
	env$W2 = tt3 / max(tt3)
	
}else{
	tt3 = env$W2 *sum(env$dset$Freq)/sum(env$W2)
}

	if(nc > nlevels(env$dset[,nv])){
		width.cor = (ntc-1)*env$S[[1]][1]/sum(env$S[[1]])*col.gap.prop
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
	
	#alpha = spread(env$W2,ncol=ntc)
	alpha = spread(exp(2*(env$W2-1)),ncol=ntc)*base.alpha
	alpha[alpha < min.alpha] = min.alpha
	if(!is.null(expected)){
		alpha[alpha < 1] = 1	
	}
	colv=col

# ----- preparing the color matrix ----- #	
if(init | env$modchange){
			if(is.null(expected)){
				#if(length(col) < ntc0){
				#	num.col = orig.ntc^(spine | !eqwidth) #ntc0

				#	if( col %in% c("hsv","rgb") ){
				#		colv = rainbow(num.col,s = 0.8, v = 0.8, start = 0, end = max(1, ntc0 - 1)/ntc0, gamma = 1, alpha = 0.9)
				#	}else{
						#colv = rainbow_hcl(num.col,c = 70, l = 70, start = 0, end = 360 * (ntc0 - 1)/ntc0, gamma = 2.4)
				#		colv = rev(sequential_hcl(num.col,h = 180, c. = c(80, 30), l = c(30, 90)))
					
				#	}
				if(any(c("hsv","hcl","rgb","seq","sequential","sqn","sqt","div","diverging","diverge") %in% col)){
				tvind = which(names(env$orig.dset) == env$fterms[env$nv])
				num.col = nlevels(env$orig.dset[,tvind])^(spine | !eqwidth)#ntc0
				if( col %in% c("hsv","rgb") ){
					col.def = formals(rainbow)
					if(!("s" %in% labels(col.opt))){
						col.opt$s = eval(col.def$s)
					}
					if(!("v" %in% labels(col.opt))){
						col.opt$v = eval(col.def$v)
					}
					if(!("start" %in% labels(col.opt))){
						col.opt$start = eval(col.def$start)
					}
					if(!("end" %in% labels(col.opt))){
						col.opt$end = max(num.col-1,1)/num.col
					}
					if(!("alpha" %in% labels(col.opt))){
						col.opt$alpha = eval(col.def$alpha)
					}
					colv = rainbow(num.col,s = col.opt$s, v = col.opt$v, start = col.opt$start, end = col.opt$end, gamma = 1, alpha = col.opt$alpha)
					}
				if( col == "hcl" ){
					col.def = formals(rainbow_hcl)
					if(!("c" %in% labels(col.opt))){
						col.opt$c = eval(col.def$c)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l = eval(col.def$l)
					}
					if(!("start" %in% labels(col.opt))){
						col.opt$start = eval(col.def$start)
					}
					if(!("end" %in% labels(col.opt))){
						col.opt$end = 360 * (num.col - 1)/num.col
					}
				colv = rainbow_hcl(num.col,c = col.opt$c, l = col.opt$l, start = col.opt$start, end = col.opt$end, gamma = 2.4)
				}
				if( col %in% c("seq","sqt","sqn","sequential") ){
					col.def = formals(sequential_hcl)
					if(!("h" %in% labels(col.opt))){
						col.opt$h = eval(col.def$h)
					}
					if(!("c" %in% labels(col.opt))){
						col.opt$c = eval(col.def$c.)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l = eval(col.def$l)
					}
					if(!("power" %in% labels(col.opt))){
						col.opt$power = eval(col.def$power)
					}
					colv = rev(sequential_hcl(num.col,h = col.opt$h, c. = col.opt$c, l = col.opt$l, power = col.opt$power))
				}
				if( col %in% c("div","diverging","diverge") ){
					col.def = formals(diverge_hcl)
					if(!("h" %in% labels(col.opt))){
						col.opt$h = eval(col.def$h)
					}
					if(!("c" %in% labels(col.opt))){
						col.opt$c = eval(col.def$c)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l = eval(col.def$l)
					}
					if(!("power" %in% labels(col.opt))){
						col.opt$power = eval(col.def$power)
					}
					colv = diverge_hcl(num.col,h = col.opt$h, c = col.opt$c, l = col.opt$l, power = col.opt$power)
					
				}
				}else{
					colv = col[1:ntc0]	
				}	
				#cat("cat.ord=",env$cat.ord)
					
				#cat("env$ntc=",env$ntc)
				if(!spine){
					colv = colv[env$cat.ord^(!eqwidth)]
				}	
				env$colv = colv		
				env$col.Mat = spread( t(rep(colv , nc/ntc)),nrow=nr)

				
			}else{
				resid.mat[is.nan(resid.mat)] = 0
				resid.mat[is.na(resid.mat)] = 0
				if(is.null(env$max.resid.scale)){
					resid.max = max(abs(resid.mat),na.rm=T)	
				}	
			
				resid.mat[resid.mat < -resid.max] = -resid.max
				resid.mat[resid.mat > resid.max]  =  resid.max
				rmax = ifelse(resid.max == 0, 1, resid.max)
				env$col.Mat = apply(resid.mat/rmax,2,function(x){
						qx = abs(x)
						#if(cut.rv){ qx =  ceiling(qx * cut.rs)/cut.rs }
					qx = qx - qx %% 1/(cut.rs-1)
					alf = (exp(qx - 1)-exp(-1))/(1-exp(-1))
					rgb(sign(x) < 0,0,sign(x) > 0,alpha = alf*base.alpha ) #qx # cut.rs/2
					
					})
			}
}

# ----- plotting itself ----- #
# 	>>> divided into eqwidth = T | eqwidth = F >>> freq.trans = Id/sqrt/log
# 	>>> in spineplot mode, the target categories will be plotted one after another (see the for-loop)

	if(eqwidth){
		if(spine){
			env$W2 = apply(env$W2,c(1,2),function(x) return(min(x/freq.trans,1)))
		}
		
			draw(env$H2/nrow(env$H2)*(1-row.gap.prop),env$W2/ncol(env$W2)*(1-col.gap.prop) + (round(env$W2,digits=7) > 0)*width.cor,env$X2,env$Y2,alpha = 0.25, bg = "grey")
			draw(env$H2/nrow(env$H2)*(1-row.gap.prop),env$H2/ncol(env$H2)*(1-col.gap.prop) + width.cor,env$X2,env$Y2,alpha = 0.1, bg = "grey")
	
		if( spine ){
			for( i in 1:length(cat.ord) ){
				if(i > 1){ 
					env$Y = env$Y + env$H/nrow(env$H)*(1-row.gap.prop)
				}
				xind = seq( cat.ord[i], nc*ntc0, ntc0 )
				env$H = env$H0[,xind]
			 	draw(env$H/nrow(env$H)*(1-row.gap.prop)*max.scale,env$W/ncol(env$W)*(1-col.gap.prop),env$X,env$Y,alpha = ifelse(is.null(expected),0.9,1), bg = env$col.Mat[,xind])
			}
			env$H = env$H0
			env$Y = env$Y0
		}else{
			env$H = apply(env$H,c(1,2),function(x) return(min(x/max.scale,1)))
			
			draw(env$H/nrow(env$H)*(1-row.gap.prop),env$W/ncol(env$W)*(1-col.gap.prop),env$X,env$Y,alpha = alpha, bg = env$col.Mat)
		}
	}else{
		if(is.numeric(freq.trans)){
			
			env$W2num = apply(tt3,c(1,2),function(x) return(min(x/(max(tt3)*freq.trans),1)))
			env$W3 = spread(env$W2num,ncol=ntc)
			env$X3 = spread(env$X[,seq(1,nc,ntc),drop=FALSE],ncol=ntc) + t( t(env$W3) * c(0:(ntc-1)) )/ncol(env$X)*(1-col.gap.prop) 
		
			draw(env$H2/nrow(env$H2)*(1-row.gap.prop),env$W2num/ncol(env$W2num)*(1-col.gap.prop) + (round(env$W2num,digits=7) > 0)*width.cor,env$X2,env$Y2,alpha = 0.25, bg = "grey")
			draw(env$H2/nrow(env$H2)*(1-row.gap.prop),env$H2/ncol(env$H2)*(1-col.gap.prop) + width.cor,env$X2,env$Y2,alpha = 0.1, bg = "grey")
			if( spine ){
				for( i in 1:length(cat.ord) ){
					if(i > 1){ 
						env$Y = env$Y + env$H/nrow(env$H)*(1-row.gap.prop)
					}
					xind = seq( cat.ord[i], nc*ntc0, ntc0 )
					env$H = as.matrix(env$H0[,xind])
					if( nr < 2 ){
						env$H = t(env$H)	
					}
					draw(env$H/nrow(env$H)*(1-row.gap.prop),env$W3/ncol(env$W3)*(1-col.gap.prop),env$X3,env$Y,alpha = 1, bg = env$col.Mat[,xind])
				}
				env$H = env$H0
				env$Y = env$Y0
			}else{
			
				env$H = apply(env$H,c(1,2),function(x) return(min(x/max.scale,1)))
			
				draw(env$H/nrow(env$H)*(1-row.gap.prop),env$W3/ncol(env$W3)*(1-col.gap.prop),env$X3,env$Y,alpha = 1, bg = env$col.Mat)
			}
		}
	}
	if(yaxis & max.scale > 0){
		sapply( env$y,function(z){
			grid.yaxis(at =seq(z,z+1/length( env$y)*(1-row.gap.prop),1/length( env$y)*(1-row.gap.prop)/5), label = round(seq(0,max.scale,max.scale/5),digits = 3), main = TRUE,
			edits = NULL, name = NULL,
			gp = gpar(), draw = TRUE, vp = NULL)
		})
		sapply( env$y,function(z){
			grid.yaxis(at =seq(z,z+1/length( env$y)*(1-row.gap.prop),1/length( env$y)*(1-row.gap.prop)/5), label = round(seq(0,max.scale,max.scale/5),digits = 3), main = FALSE,
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
				grid.text(  names(env$dset)[cind[i]],x = 0.5 , y = 5/6, just = "centre",gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
			}
			#z = (nlvl[cind])[i]*cur
			z = zc[i]
			grid.text(  rep( clabs[[i]] ,cur), x = env$X[1,seq(1,nc,nc/z)] + (env$X[1,seq(nc/z,nc,nc/z)] - env$X[1,seq(1,nc,nc/z)] +1/nc*(1-col.gap.prop))/2, y = 1/3, just = "centre",gp=gpar(cex=lab.cex)) 
			if( boxes ){
				grid.rect(  y = 1/3 , x = env$X[1,seq(1,nc,nc/z)] , just = c("left","centre"), width = env$X[1,seq(nc/z,nc,nc/z)] - env$X[1,seq(1,nc,nc/z)] +1/nc*(1-col.gap.prop), height = 1/2, gp = gpar(fill="transparent"))
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
				grid.text(  names(env$dset)[rind[i]],x = 1/6 , y = 0.5, just = "centre",rot=90,gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
			}			
			#z = (nlvl[rind])[i]*cur
			z = zr[i]
			grid.text(  rep( rlabs[[i]] ,cur), y = env$Y[seq(1,nr,nr/z),1] + (env$Y[seq(nr/z,nr,nr/z),1] - env$Y[seq(1,nr,nr/z),1] +1/nr*(1-row.gap.prop))/2, x = 2/3, just = "centre",rot=90,gp=gpar(cex=lab.cex)) 
			if( boxes ){
				grid.rect(  x = 2/3 , y = env$Y[seq(1,nr,nr/z),1] , just = c("centre","bottom"), height = env$Y[seq(nr/z,nr,nr/z),1] - env$Y[seq(1,nr,nr/z),1] +1/nr*(1-row.gap.prop), width = 1/2, gp = gpar(fill="transparent"))
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
		alfav = (exp(seq(0,-1,-1/(cut.rs-1)))-exp(-1))/(1-exp(-1))
		colv = rgb(rep(c(1,0),each=cut.rs),0,rep(c(0,1),each=cut.rs),alpha = c(alfav,rev(alfav))*base.alpha)
	
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
	if(dset.changed){
		env$dset = env$orig.dset
	}
if(init | env$spinechange){	

	env$freq.trans = freq.trans	
	env$eqwidth = eqwidth
	env$spine = spine
	env$cat.ord = cat.ord
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









irmb = function(x,...){
	UseMethod("irmb")
}





irmb.table = function(x,   gap.prop = 0.2,
               gap.mult = 1.5, col = "hcl", col.opt=list(), Z=100, abbrev = FALSE, use.na = FALSE,
				 expected = TRUE, resid.type = "pearson", 
               max.resid.scale = NULL, cut.rs = 7, yaxis = TRUE, min.alpha = 0.1,
               boxes = TRUE, lab.tv = FALSE, varnames = TRUE,
               lab.cex = 1.2,...){
	
	formula = as.formula(paste("~",do.call(paste,c(as.list(names(dimnames(x))),sep="+")),sep=""))
	data = as.data.frame(ftable(x))
	irmb.formula(formula, data, 
                 use.na=use.na ,expected= expected ,
                 resid.type=resid.type ,
               resid.max=resid.max , cut.rs= cut.rs ,                
                gap.prop = gap.prop , gap.mult = gap.mult , col = col, col.opt = col.opt ,yaxis = yaxis , 
               min.alpha =  min.alpha , boxes = boxes , lab.tv = lab.tv ,
               varnames = varnames , abbrev = abbrev ,lab.cex = lab.cex)
}

irmb.ftable = function(x , gap.prop = 0.2,
               gap.mult = 1.5, col = "hcl", col.opt=list(), Z=100, abbrev = FALSE, use.na = FALSE, 
				expected = TRUE,  resid.type = "pearson", 
               max.resid.scale = NULL, cut.rs = 7, yaxis = TRUE, min.alpha = 0.1,
               boxes = TRUE, lab.tv = FALSE, varnames = TRUE,
               lab.cex = 1.2,...){
	rv = names(attr(x,"row.vars"))
	cv = names(attr(x,"col.vars"))
	col.vars = NULL
	
		
	formula = as.formula(paste("~",do.call(paste,c(as.list(c(rv,cv)),sep="+")),sep=""))
	data = as.data.frame(x)
	irmb.formula(formula, data, 
                 use.na=use.na ,expected= expected ,
                 resid.type=resid.type ,
               resid.max=resid.max  , cut.rs= cut.rs ,                
                gap.prop = gap.prop , gap.mult = gap.mult , col = col, col.opt = col.opt,
				yaxis = yaxis, min.alpha =  min.alpha , boxes = boxes , lab.tv = lab.tv ,
               varnames = varnames , abbrev = abbrev ,lab.cex = lab.cex)
}




 irmb.formula = function(formula, data, gap.prop = 0.2,
               gap.mult = 1.5, col = "hcl", col.opt =list(), Z=100, abbrev = FALSE, use.na = FALSE, expected = TRUE,
               resid.type = "pearson", 
               max.resid.scale = NULL, cut.rs = 7, yaxis = TRUE, min.alpha = 0.1,
               boxes = TRUE, lab.tv = FALSE, varnames = TRUE,
               lab.cex = 1.2,...){

	if(!("iWidgets" %in% .packages(all.available=TRUE) ) ){
		cat("The required package iWidgets is not available.\n Please use \n install.packages('iWidgets',,'http://rforge.net/',type='source') \n for installation.")
		return(invisible(TRUE))
		}
	if(!("JGR" %in% .packages(all.available=FALSE) ) ){
		cat("iWidgets requires 'JGR' in order to run correctly.\n Please visit \n http://www.rosuda.org/software/ \n to download it.")
		return(invisible(TRUE))
		}
	require(iWidgets)
	
	
		
	e1 = new.env()
	w = iwindow()
	
	e1$dset = data.frame(data)
	e1$f = as.formula(formula)
	if(length(terms(e1$f))>2){
		names(e1$dset)[which( suppressWarnings(names(e1$dset) == terms(e1$f)[[2]]))] = "Freq"
	}
	
	# ----- save global parameters to e1 ----- #
	e1$mod.type = "poisson"
	e1$resid.type = resid.type
	e1$expected = expected
	e1$use.expected.values =  FALSE
	e1$cut.rv = TRUE
	e1$cut.rs = cut.rs
	e1$gap.prop = gap.prop
	e1$gap.mult = gap.mult
	e1$varnames = varnames
	e1$min.alpha = min.alpha
	e1$abbrev = abbrev
	
	e1$lab.cex = lab.cex
	e1$lab.tv = lab.tv
	e1$yaxis = yaxis
	e1$boxes = boxes
	e1$col = col
	e1$col.opt = col.opt
	if( !("base.alpha" %in% ls()) ){
			base.alpha = 1	
	}
	e1$base.alpha = base.alpha
	e1$fterms = attr(terms(e1$f),"term.labels")
	e1$genfterms = e1$fterms
	e1$nv = length(e1$fterms)

	e1$col.vars = rep(c(T,F),e1$nv)[1:e1$nv]
	e1$col.vars[which(e1$genfterms == e1$fterms[e1$nv])] = TRUE
		

	
	#e1$dset = subtable(dset,cols = which(names(dset) %in% e1$fterms),keep.zero=F)
	e1$dset = subtable(e1$dset,cols = match(e1$fterms,names(e1$dset)),allfactor=TRUE)
	e1$orig.dset=e1$dset
	
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
	e1$freq.trans = 1
	

	if(e1$expected){
		#e1$expected = as.list(1:length(attr(terms(e1$f),"term.labels")))
		e1$exp.ls = as.list( c(1:(e1$nv-1)) ,e1$nv)
		e1$interaction.terms = paste( sapply(e1$exp.ls,function(x) paste(names(e1$dset)[x],collapse="*") ),collapse="+")
		
	}else{
		e1$expected = NULL	
	}
	
	e1$eqwidth = FALSE
	e1$spine = FALSE
	e1$spinechange = FALSE
	
	e1$updateplot = TRUE
	e1$modchange = TRUE
	
	e1$lab = TRUE
	
	e1$ntc = nlevels(e1$orig.dset[,which(names(e1$orig.dset) == e1$fterms[length(e1$fterms)])])
	e1$cat.ord = 1:e1$ntc
	
	g = igroup()
	g2 = igroup(window=g,horizontal = FALSE)
	add(g2,igraphics(width=1024,height=768))
	rmbX(env = e1,init = TRUE)
		mslider = islider(min=0,max=Z,horizontal = FALSE,value = Z, window=g, handler=function(h,...){
		if(e1$updateplot){
		if(!e1$spine & !get.value(vcdcheck)){# do nothing for spine and mosaic mode
			if(exists("H",e1)){
				e1$H = e1$H*e1$max.scale
			}
			e1$max.scale = get.value(h$obj)/Z
			if((!e1$spine) & e1$updateplot){
				rmbX(env = e1 ,init = FALSE)
			}
		}
		}
		})
	tslider = islider(min=0,max=Z,horizontal = TRUE,value = Z, window=g2, handler=function(h,...){
		e1$freq.trans = get.value(h$obj)/Z
		if(e1$updateplot){
		if(exists("H",e1)){
			e1$H = e1$H*e1$max.scale
		}
		
		#if(e1$updateplot){
			rmbX(env = e1 ,init = FALSE)
		}
		})
	g3 = igroup(window=g2,horizontal = FALSE)
	g3a = igroup(window=g3,horizontal = TRUE)
	icheckbox("spineplot", checked = FALSE, window = g3a,handler=function(h,...){
			e1$spine = !e1$spine
			e1$spinechange = TRUE
			#ntc = nlevels(e1$orig.dset[,which(names(e1$orig.dset) == e1$fterms[length(e1$fterms)])])
			if(max(e1$cat.ord) > e1$ntc){
				e1$cat.ord = 1:e1$ntc
				set.value(hledit,paste(1:e1$ntc,collapse=","))
			}
			e1$updateplot = FALSE
			set.value(tslider,Z)
			set.value(mslider,Z)
			e1$max.scale = 1
			e1$freq.trans = 1
			e1$updateplot = TRUE
			rmbX(env = e1 ,init = FALSE)
			e1$spinechange = FALSE
			set.value(vcdcheck,FALSE)
			})
	icheckbox("eqwidth", checked = e1$eqwidth, window = g3a,handler=function(h,...){
			e1$eqwidth = !e1$eqwidth
			
			set.value(tslider,Z)
			
			e1$freq.trans = 1
			
			rmbX(env = e1 ,init = FALSE)
			set.value(vcdcheck,FALSE)
			})
	icheckbox("lab", checked = e1$lab, window = g3a,handler=function(h,...){
			e1$lab = !e1$lab
			if(get.value(vcdcheck)){
				vcdX(e1)
			}else{
				rmbX(env = e1 ,init = FALSE)
			}
			})
	vcdcheck = icheckbox("mosaic", checked = FALSE, window = g3a,handler=function(h,...){
			vcdmosaic = get.value(h$obj)
			if(!("vcd" %in% .packages(all.available=TRUE) ) ){
				cat("This option requires the package vcd. Please use\n install.packages('vcd')\n to install it and try again.")
				set.value(vcdcheck,FALSE)
			}else{
		require(vcd)
			if(get.value(vcdcheck)){
				
			if(length(e1$cat.ord) > 1){
				#reverse order of target categories in vcd
				tvind = which(names(e1$dset) == e1$fterms[e1$nv])
				e1$dset[,tvind] = factor(e1$dset[,tvind], 
					levels = levels(e1$orig.dset[,which(names(e1$orig.dset) == e1$fterms[e1$nv])])[rev(e1$cat.ord)])
				e1$colv=rev(e1$colv)
					vcdX(e1)
				#reset the target variable in dset and the color vector colv
				
#				levels(e1$dset[,e1$nv])[which(! 1:e1$nv %in% e1$cat.ord )] = NA
#				e1$dset[,e1$nv] = factor(e1$dset[,e1$nv],levels=levels(e1$dset[,e1$nv])[rank(rev(e1$cat.ord))])
#				e1$colv=rev(e1$colv)
				e1$dset[,tvind] = factor(e1$dset[,tvind], 
					levels = levels(e1$orig.dset[,which(names(e1$orig.dset) == e1$fterms[e1$nv])])[e1$cat.ord])
				e1$colv=rev(e1$colv)
				
				if(exists("H",e1)){
					e1$H = e1$H*e1$max.scale
				}
				e1$updateplot = FALSE
				set.value(tslider,Z)
				set.value(mslider,Z)
				e1$max.scale = 1
				e1$freq.trans = 1
				e1$updateplot = TRUE
				
			}else{
				set.value(vcdcheck,FALSE)
			}
			}else{
				rmbX(env = e1 ,init = FALSE)
			}
		}
			})
			
			
	if(!is.null(e1$expected)){
		expectedcheck1 = icheckbox("res. shading", checked = TRUE, window = g3a,handler=function(h,...){
			
			if(!get.value(h$obj)){
				 e1$expected = NULL
			}else{
				e1$expected = TRUE
			}
			rmbX(env = e1 ,init = TRUE)
			if(get.value(vcdcheck)){
				vcdX(e1)
			}
		
			})
		expectedcheck2 = icheckbox("exp. values", checked = FALSE, window = g3a,handler=function(h,...){
			expcheck = get.value(h$obj)
			e1$use.expected.values = expcheck
			rmbX(env = e1 ,init = TRUE)
			if(get.value(vcdcheck)){
				vcdX(e1)
			}
		
			})
		
	}
			
	g3label = igroup(window=g3,horizontal=TRUE)
	fx = ilabel(text = paste(c("Current formula: >> ",as.character(e1$f)," <<"),collapse=""),window = g3label)

	g4 = igroup(window=g3,horizontal=TRUE)
	add.space(g4,32)
	ilabel("target category order:     ",window=g4)
	add.space(g4,10)
	
	hledit = iedit(window=g4,text = paste(e1$cat.ord,collapse=","), handler=function(h,...){
		e1$dset = e1$orig.dset
		tmp = e1$cat.ord
		op=options()
		#e1$dset=e1$orig.dset #reset the dataset
		e1$cat.ord <- lapply(strsplit(get.value(h$obj),"\\,"),as.integer  )[[1]]
		options(error = function(){
				e1$cat.ord = tmp
				set.value(h$obj,paste(tmp,collapse=","))
			})
			
		#if(e1$spine & e1$updateplot){
		#if(e1$updateplot){
			
			#ntc = nlevels(e1$orig.dset[,which(names(e1$orig.dset) == e1$fterms[length(e1$fterms)])])
			if(max(e1$cat.ord) > e1$ntc){
				e1$cat.ord = 1:e1$ntc
				set.value(hledit,paste(1:e1$ntc,collapse=","))
			}
			# this is inefficient
			rmbX(env = e1 ,init = TRUE)
			if(get.value(vcdcheck)){
				vcdX(e1)
			}
			
			options(op)
			#}
		})
		add.space(g4,80)
		
	if(expected){
		g5 = igroup(window=g3,horizontal=T)
		add.space(g5,10)
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
						e1$modchange = TRUE
						rmbX(env = e1 ,init = FALSE)
						options(op)
						e1$modchange = FALSE
						if(get.value(vcdcheck)){
							vcdX(e1)
						}
				}
			})
			add.space(g5,20)
	}


g5 = igroup(window = g2, horizontal = TRUE)
g5lab = igroup(window = g5, horizontal = FALSE)
g5ord = igroup(window = g5, horizontal = FALSE)
g5sel = igroup(window = g5, horizontal = FALSE)
g5split = igroup(window = g5, horizontal = FALSE)

ilabel(text = "variable",window = g5lab)
ilabel(text = "position",window = g5ord)
ilabel(text = "sel",window = g5sel)
ilabel(text = "horiz",window = g5split)
add.space( g5lab,10)
add.space( g5ord,10)
add.space( g5sel,10)
add.space( g5split,10)

	
	e1$gencol.vars = e1$col.vars
	
	#add.space(g5,97)
	#ilabel(text = "variables:",window = g5)
	
	
nchar.max = max(sapply(e1$genfterms,nchar))
e1$var.labels = sapply(e1$genfterms,function(x){
		paste(c(x,rep(" ",nchar.max-nchar(x)),":"),collapse="")
	})

	#add.space(g5,20)
e1$var.pos = c(1:e1$nv)
e1$var.sel = rep(TRUE,e1$nv)
radiobutton = list()
sel.check = list()
split.check = list()
for(i in 1:e1$nv){
	add.space( g5lab,5)
	ilabel(text = e1$var.labels[i],window = g5lab)
	add.space( g5lab,10)
	radiobutton[[i]] =  iradio(1:e1$nv, window=g5ord,horizontal=TRUE, name = e1$genfterms[i], id = i,selected=i,
          handler=function(h,...){
				buttonid = h$obj$id
			newval = get.value(h$obj)
			oldval = e1$var.pos[buttonid]
			
			
			# compute new positions:
			othervalues = sapply(radiobutton,function(x){
						get.value(x)
					})
			if(newval > oldval){	
				for(s in which(e1$var.pos > oldval & e1$var.pos <= newval)){
					set.value(radiobutton[[s]],othervalues[s]-1)
					e1$var.pos[s] = e1$var.pos[s]-1
				}
			}else{
				for(s in which(e1$var.pos >= newval & e1$var.pos < oldval)){
					set.value(radiobutton[[s]],othervalues[s]+1)
					e1$var.pos[s] = e1$var.pos[s]+1
				}
			}
			e1$var.pos[buttonid] = newval
			
			
			
			# change fterms, cat.ord, col.vars, ntc and redraw
			
			# get new fterms and formula from changed var.pos
			curr.ord = order(e1$var.pos)
			e1$fterms = e1$genfterms[curr.ord][which(e1$var.sel[curr.ord])]
			tmp = paste("~",paste(e1$fterms,collapse="+"),sep="")
			e1$f = as.formula(tmp)
			
			e1$ntc = nlevels(e1$orig.dset[,which(names(e1$orig.dset) == e1$fterms[length(e1$fterms)])])
			
			#change hledit and cat.ord according to orig.dset and reset dset
			tvind = which(names(e1$orig.dset) == e1$fterms[e1$nv])
			ntv = nlevels(e1$orig.dset[,tvind])
			e1$cat.ord = 1:ntv
			set.value(hledit,paste(1:ntv,collapse=","))
			e1$col.vars[tvind] = TRUE
			set.value(split.check[[tvind]],TRUE)
			e1$dset = e1$orig.dset
						
			set.value(fx, paste("      Current formula: >> ",tmp," <<   "))
			if(expected){
					e1$interaction.terms = paste(c(paste(e1$fterms[1:(e1$nv-1)],collapse="*"),e1$fterms[e1$nv]),collapse="+")
					set.value(modedit,e1$interaction.terms)
					e1$modchange = TRUE
				}
				
			rmbX(env = e1,init =TRUE)
			if(get.value(vcdcheck)){
				vcdX(e1)
			}
			})
	sel.check[[i]] = icheckbox("", checked= TRUE, window = g5sel,id = i, handler=function(h,...){
			buttonid = h$obj$id
			e1$var.sel[buttonid] = !e1$var.sel[buttonid]
			
			reset.cat.ord = e1$var.pos[buttonid] %in% c(e1$nv,length(e1$genfterms))
						
			curr.ord = order(e1$var.pos)
			e1$fterms = e1$genfterms[curr.ord][which(e1$var.sel[curr.ord])]
			e1$nv = e1$nv - (-1)^(get.value(h$obj))
			
			#if the last variable has changed we have to change hledit, col.vars and the dset
			#e1$nv has changed to the new value
			if(reset.cat.ord){
				tvind = which(names(e1$orig.dset) == e1$fterms[e1$nv])
				ntv = nlevels(e1$orig.dset[,tvind])
				e1$cat.ord = 1:ntv
				set.value(hledit,paste(1:ntv,collapse=","))
				e1$col.vars[tvind] = TRUE
				set.value(split.check[[tvind]],TRUE)
				e1$dset = e1$orig.dset
				e1$ntc = nlevels(e1$orig.dset[,which(names(e1$orig.dset) == e1$fterms[length(e1$fterms)])])
			}
			
			if(length(e1$fterms) > 0){
				tmp = paste("~",paste(e1$fterms,collapse="+"),sep="")
				e1$f = as.formula(tmp)
				set.value(fx, paste("Current formula: >> ",tmp," <<"))
				e1$updateplot = FALSE
				set.value(hledit,paste(1:nlevels(e1$orig.dset[,which(names(e1$orig.dset) == e1$fterms[e1$nv])]),collapse=","))
				
				set.value(tslider,Z)
				set.value(mslider,Z)
				e1$max.scale = 1
				e1$freq.trans = 1
				e1$updateplot = TRUE

				if(expected){
					e1$interaction.terms = paste(c(paste(e1$fterms[1:(e1$nv-1)],collapse="*"),e1$fterms[e1$nv]),collapse="+")
					set.value(modedit,e1$interaction.terms)
					e1$modchange = TRUE
				}
				rmbX(env = e1,init =TRUE)
				e1$modchange = FALSE
				if(get.value(vcdcheck)){
					vcdX(e1)
				}
			}
		})
	split.check[[i]] = icheckbox("", checked= e1$col.vars[i], window = g5split,id = i, handler=function(h,...){
			buttonid = h$obj$id
			tvind = which(names(e1$orig.dset) == e1$fterms[e1$nv])
		if(buttonid != tvind){ #the target variable is always horizontal and cannot be changed 
			e1$col.vars[buttonid] = !e1$col.vars[buttonid]
			


	#redraw if variable is in current selection
			if(e1$var.sel[buttonid]){
				# this could be more efficient:
				rmbX(env = e1,init = TRUE)
				if(get.value(vcdcheck)){
					vcdX(e1)
				}
			}
			
		}else{
		 	set.value(split.check[[buttonid]],TRUE)
		}
		})
}
	add(w,g)
	visible(w,TRUE)
}

# ------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------ #


space = function(nlvl,vsplit,gap.mult = 1.1,gap.prop = 0.1,last.col.zero = F, last.row.zero = F){
	# 	>>> function to compute the gaps widths for the rmb plot
	# 	>>> might be extended in further versions of the package
	
	# ----- parameter check ----- #
		if( !(  (gap.prop > 0)&(gap.prop < 1) ) ){
			stop(simpleError("Wrong gap.prop specification!"))
		}
	
	col.nlvl = nlvl[vsplit]
	row.nlvl = nlvl[!vsplit]								

	nc = prod(col.nlvl)
	nr = prod(row.nlvl)

	ncv = length(col.nlvl)
	nrv = length(row.nlvl)

	col.spaces = rep(gap.prop,nc-1)
	row.spaces = rep(gap.prop,nr-1)

	if( ncv > 1 ){
		cind = 1:(nc-1)
		cur = 1
		for( i in (ncv-1):1 ){ 
			cur = cur*col.nlvl[i+1]
			curind = which( (cind %% cur) == 0)
			col.spaces[ curind  ] = col.spaces[ curind  ] * gap.mult
		}
	}
	if( nrv > 1 ){
		rind = 1:(nr-1)
		cur=1
		for( i in (nrv-1):1 ){ 
			cur = cur*row.nlvl[i+1]
			curind = which( (rind %% cur) == 0)
			row.spaces[ curind  ] = row.spaces[ curind  ] * gap.mult
		}
	}
	if(last.col.zero){
		col.spaces[which(col.spaces == gap.prop)] = 0
		col.spaces = col.spaces/gap.prop
	}
	if(last.row.zero){
		row.spaces[which(row.spaces == gap.prop)] = 0
		row.spaces = row.spaces/gap.prop
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


vcdX = function(env,...){
	if(!is.null(env$expected)){
		exp.ls = lapply(strsplit(attr(terms(as.formula(paste("~",env$interaction.terms,collapse=""))),"term.labels"),":"),function(x) match(x,env$genfterms))
		exp.ls = c(exp.ls,1:env$nv)
	}else{
		exp.ls = NULL
	}
	
	sv = env$col.vars[ match(env$fterms, env$genfterms)  ]
	
	sv[env$nv] = FALSE
	
		if(env$lab){
				vcdlab = labeling_border(boxes = TRUE,labels= c(rep(TRUE,env$nv-1),FALSE))
				vcdmargins = NULL 
		}else{
			vcdlab = FALSE
			vcdmargins = 1
		}
		
		
	vcdtype = ifelse(env$use.expected.values,"expected","observed")
	if(!is.null(env$expected)){
	mosaic(formula = env$f,split_vertical=sv, data = env$dset,
			main = "", shade = TRUE, legend = TRUE, expected = exp.ls,labeling = vcdlab,margins = vcdmargins,type=vcdtype)	
	}else{
	mosaic(formula = env$f,split_vertical=sv, data = env$dset,
			main = "",  expected = exp.ls,
			labeling = vcdlab,margins = vcdmargins,type=vcdtype, gp = gpar(fill =env$colv))	
		
		}
	
	return(invisible(TRUE))
	
}

