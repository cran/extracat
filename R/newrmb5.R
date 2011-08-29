
# -------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------- #

draw = function(H,W,X,Y,alpha = 1,border = "black",bg = "white",vp=NULL, lwd = 1){
	grid.rect(x =  unit(X, "npc"), y = unit(Y, "npc"),
          width = unit(W, "npc"), height = unit(H, "npc"),
          just = c("left","bottom"),
          default.units = "npc", name = NULL,
          gp=gpar(col=border,fill = bg, alpha = alpha, lwd = lwd), draw = TRUE, vp = vp
	)
}


rmb = function(x,...){
	UseMethod("rmb")
}
rmb.table = function(x, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
 model.opt = list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(),...){
	
	vn <- names(dimnames(x))
	if(is.numeric(col.vars)){
			vn <- c(vn[-col.vars],vn[col.vars])
			col.vars <- sort(seq_along(	vn ) %in% col.vars)
	}
	
	formula <- as.formula(paste("~",paste(vn,collapse="+"),sep=""))
	
	data <- as.data.frame(ftable(x))
	rmb.formula(formula, data, col.vars = col.vars, spine = spine, circular = circular,
			 eqwidth = eqwidth, cat.ord = cat.ord, freq.trans = freq.trans, max.scale = max.scale,
			 use.na = use.na, expected = expected, residuals = residuals, model.opt = model.opt, 
			 gap.prop = gap.prop, gap.mult = gap.mult, col = col,col.opt = col.opt, label = label,
			 label.opt = label.opt)
}

rmb.ftable = function(x, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
 model.opt = list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(),...){
	
	rv <- names(attr(x,"row.vars"))
	cv <- names(attr(x,"col.vars"))
	col.vars <- NULL
	#if(is.null(col.vars)){
		col.vars <- c(rep(F,length(rv)),rep(T,length(cv)))	
	#}
	formula <- as.formula(paste("~",do.call(paste,c(as.list(c(rv,cv)),sep="+")),sep=""))
	data <- as.data.frame(x)
	rmb.formula(formula, data, col.vars = col.vars, spine = spine, circular = circular,
			 eqwidth = eqwidth, cat.ord = cat.ord, freq.trans = freq.trans, max.scale = max.scale,
			 use.na = use.na, expected = expected, residuals = residuals, model.opt = model.opt, 
			 gap.prop = gap.prop, gap.mult = gap.mult, col = col,col.opt = col.opt, label = label,
			 label.opt = label.opt)
}


rmb.glm = function(x, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
 model.opt = list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(),...){
	
	model.family <- x$family$family
	stopifnot(model.family %in% c("poisson","binomial"))
	
	if(model.family == "poisson"){
	
		fi <- which(names(x$data) == "Freq")
		
		tms <- terms(x$formula)
		vars <- attr(tms,"term.labels")[ which(attr(tms,"order") == 1)]
		formula <- paste( toString(x$formula[2]),"~", paste(vars,collapse="+"),sep="" ) 
		
		if("resid.type" %in% names(model.opt)){
			rtp <- model.opt$resid.type	
		}else{
			rtp <- "pearson"
		}
		data <- subtable(x$data, match(vars,names(x$data)[-fi]), keep.zero=TRUE)
	
		x <- update(x, data = data)
		residuals <- residuals(x,type=rtp)
		expected <- x$fitted

	}
	if(model.family == "binomial"){
		
		tms <- terms(x$formula)
		vars <- attr(tms,"term.labels")[ which(attr(tms,"order") == 1)]
		vars <- c(vars[-1],vars[1])
		formula <- paste( "Freq ~ ",toString(x$formula[2]), paste(vars,collapse="+"),sep="" ) 
		
		nv <- length(vars)+1
		
		unweighted <- all(x$prior.weights==1)
		data <- subtable(x$data, 1:nv, keep.zero=TRUE, allfactor = TRUE)
		
		if(unweighted){
			x <- update(x, data = data, weights = data$Freq)
		}else{
			x <- update(x, data = data)
		}
		if("resid.type" %in% names(model.opt)){
			rtp <- model.opt$resid.type	
		}else{
			rtp <- "pearson"
		}
		
		residuals <- residuals(x,type=rtp)
		expected <- x$fitted
	}
	
	rmb.formula(formula, data, col.vars = col.vars, spine = spine, circular = circular,
			 eqwidth = eqwidth, cat.ord = cat.ord, freq.trans = freq.trans, max.scale = max.scale,
			 use.na = use.na, expected = expected, residuals = residuals, model.opt = model.opt, 
			 gap.prop = gap.prop, gap.mult = gap.mult, col = col,col.opt = col.opt, label = label,
			 label.opt = label.opt)
}



#rmb.polr = function(x, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
#                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
# model.opt <- list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(),...){
#	
#	model.family <- x$family$family
#	stopifnot(model.family %in% c("poisson","binomial"))
#	
#
#		vars <- dimnames(attr(x$terms,"factors"))[[1]]
#		vars <- c(vars[-1],vars[1])
#		nv <- length(vars)
#		formula <- paste( "Freq~", paste(vars,collapse="+"),sep="" ) 
#		
#		data <- x$model
#		
#		unweighted <- ncol(data) == nv
#		
#		if(unweighted){
#			data <- subtable(data, 1:nv, keep.zero=FALSE, allfactor = TRUE)
#			x <- update(x, data = data, weights = data$Freq)
#		}else{
#			names(data)[which(!(names(data) %in% vars))] <- "Freq"	
#		}
#		
#		residuals <- residuals(x,type=resid.type)
#		expected <- x$fitted
#	
#	rmb.formula(formula, data, col.vars = col.vars, spine = spine, circular = circular,
#			 eqwidth = eqwidth, cat.ord = cat.ord, freq.trans = freq.trans, max.scale = max.scale,
#			 use.na = use.na, expected = expected, residuals = residuals, model.opt = model.opt, 
#			 gap.prop = gap.prop, gap.mult = gap.mult, col = col,col.opt = col.opt, label = label,
#			 label.opt = label.opt)
#}




rmb.formula = function(formula, data, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
 model.opt = list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(),...){
	
# use.expected.values = FALSE, mod.type = "poisson", resid.type = "pearson",
# resid.display = "both",max.rat = 2.0, resid.max = NULL, cut.rs = 5
# boxes = TRUE, lab.tv = FALSE, varnames = TRUE, abbrev = FALSE, lab.cex = 1.2, yaxis = TRUE

min.alpha <- 0.1

# ----- get model parameters from model.opt or use defaults ------ #
if( "use.expected.values" %in% names(model.opt) ){
	use.expected.values <- model.opt$use.expected.values
}else{
	use.expected.values <- FALSE
}
if( "mod.type" %in% names(model.opt) ){
	mod.type <- model.opt$mod.type
}else{
	mod.type <- "poisson"
}
if( "resid.type" %in% names(model.opt) ){
	resid.type <- model.opt$resid.type
}else{
	resid.type <- "pearson"
}
if( "resid.display" %in% names(model.opt) ){
	resid.display <- model.opt$resid.display
}else{
	resid.display <- "both"
}
if( "max.rat" %in% names(model.opt) ){
	max.rat <- model.opt$max.rat
}else{
	max.rat <- 2.0
}
if( "resid.max" %in% names(model.opt) ){
	resid.max <- model.opt$resid.max
}else{
	resid.max <- NULL
}
if( "cut.rs" %in% names(model.opt) ){
	cut.rs <- model.opt$cut.rs
}else{
	cut.rs <- 5
}

# ----- get label parameters from label.opt or use defaults ------ #
if( "yaxis" %in% names(label.opt) ){
	yaxis <- label.opt$yaxis
}else{
	yaxis <- TRUE
}
if( "boxes" %in% names(label.opt) ){
	boxes <- label.opt$boxes
}else{
	boxes <- TRUE
}
if( "lab.tv" %in% names(label.opt) ){
	lab.tv <- label.opt$lab.tv
}else{
	lab.tv <- FALSE
}
if( "varnames" %in% names(label.opt) ){
	varnames <- label.opt$varnames
}else{
	varnames <- TRUE
}
if( "abbrev" %in% names(label.opt) ){
	abbrev <- label.opt$abbrev
}else{
	abbrev <- FALSE
}
if( "lab.cex" %in% names(label.opt) ){
	lab.cex <- label.opt$lab.cex
}else{
	lab.cex <- 1.2
}




# ----- parameter check 1 ----- #
		
		dset <- data.frame(data)
		f <- as.formula(formula)
		cut.rv <- TRUE
		if(max.scale > 1 ){
			stop(simpleError("Wrong max.scale specification! (0,1]"))
		}
				
		if(!( (gap.prop < 1) & (gap.prop > 0) ) ){
			stop(simpleError("Wrong gap.prop specification! (0,1)"))
		}
		if( !((min.alpha <= 1) & (min.alpha >= 0)) ){
			stop(simpleError("Wrong alpha specification! [0,1]"))
		}
		
		stopifnot(is.logical(abbrev)|is.numeric(abbrev))
		stopifnot(is.logical(label)|is.numeric(label))
		
		if( !("base.alpha" %in% ls()) ){
			base.alpha <- 1	
		}
	
	
		
# ----- terms extraction ----- #
	fterms <- attr(terms(f),"term.labels")
	nv <- length(fterms)
	
	if(is.numeric(col.vars)){
			col.vars <- seq_along(	fterms ) %in% col.vars
		}
	
	if( !is.null(expected) ){
		both <- resid.display %in% sapply(1:4, function(z) abbreviate("both",z))
		res.val.only <- resid.display %in% c(sapply(1:6, function(z) abbreviate("values",z)),sapply(1:6, function(z) abbreviate("value",z)))
		disp.res <- both | res.val.only
		if( disp.res & all( sapply(expected, function(z) all(z == TRUE))) ){
			if( mod.type == "poisson" ){
				expected <- list(1:(nv-1),nv)	
			}else{
				expected <- as.list(1:(nv-1))	
			}
		}
		if( resid.display %in% sapply(1:6, function(z) abbreviate("values",z)) ){
			col <- "rgb"
			col.opt$alpha <- 0.3 	
		}
	}else{
		res.val.only <- FALSE	
	}
	

	
			
	nv0 <- nv
	if(is.null(col.vars)){
		col.vars <- rep(c(T,F),nv)[1:nv]	
	}
		
	col.vars[nv] <- TRUE
	tv <- fterms[nv]
	
	stopifnot(is.logical(col.vars) & length(col.vars) == nv)
	
	ind <-  match(fterms,names(dset))
# 	>>> a dummy-variable is used to handle the case of no row-variables
	if( sum(!col.vars) < 1 ){
		dset <- data.frame(probability = rep("prob",nrow(dset)),dset)
		nv<-nv+1
		col.vars <- c(F,col.vars)
		ind <- c(1,ind+1)
	}
	if( sum(col.vars) < 2 ){
		dset <- data.frame(distribution = rep(tv,nrow(dset)),dset)
		nv<-nv+1
		col.vars <- c(T,col.vars)
		ind <- c(1,ind+1)
		if("probability" %in% names(dset)){ expected <- NULL }
	}

# 	>>> preparing the labels for plotting including the abbreviation option	
	if(is.numeric(label)){
			lab <- TRUE
			tmp <- rep(FALSE,nv)
			tmp[label] <- TRUE
			label <- tmp
	}else{
			lab <- any(label)	
			if(length(label) == 1 & lab){
				label <- rep(TRUE,nv)	
			}
			if(!lab){
				label <- rep(FALSE,nv)	
			}
	}
	
	if(is.numeric(abbrev)){
			abbr <- TRUE
	}else{
			abbr <- abbrev	
	}
	
	
	
	if(!abbr){
		rclabs <- lapply(dset[,ind],function(x) levels(as.factor(x)))
	}else{
		rclabs <- lapply(dset[,ind],function(x) abbreviate(levels(as.factor(x)),abbrev))
	}
	lab.tv <- ifelse(spine,FALSE,lab.tv)

	
	label[nv] <- lab.tv
	if( sum(label) == 0 ){
		lab <- FALSE	
	}
	rlabs <- lapply(which( (!col.vars)*label > 0),function(x) rclabs[[x]])
	clabs <- lapply(which(col.vars*label > 0),function(x) rclabs[[x]])
	
#	orig.labs <- lapply(dset[,ind],function(x) levels(as.factor(x)))
	if( length(terms(f)) > 2 ){
		names(dset)[which( suppressWarnings(names(dset) == terms(f)[[2]]))] <- "Freq"
	}
	dset <- subtable(dset,ind,keep.zero=FALSE,allfactor=TRUE)
	
	# ----- predefined expected values ----- #

if( is.numeric(expected) ){
	stopifnot( length(expected) == length(dset$Freq) )
	mod.type <- "gen"
	pred <- expected
}else{
	if(suppressWarnings(max(unlist(expected))) > (nv - (mod.type=="logit")) ){
		stop(simpleError("Wrong expected specification!"))
	}	
}

	if(is.null(cat.ord)){
		ntc <- nlevels(dset[,nv])
		cat.ord <- 1:ntc
	}	
	if(length(cat.ord) > 1 & !spine){
		levels(dset[,nv])[which(! 1:nv %in% cat.ord )] <- NA
		dset[,nv] <- factor(dset[,nv],levels=levels(dset[,nv])[rank(cat.ord)])
		if(lab.tv){
			clabs[[length(clabs)]] <- levels(dset[,nv])	
		}
	}
	
	if(!use.na){
		dset <- na.omit(dset)
	}else{
		for( i in 1:(ncol(dset)-1) ){
			ind <- is.na(dset[,i])
			if(any(ind)){
				levels(dset[,i]) <- c(levels(dset[,i]),"N/A")	
				dset[which(ind),i] <- "N/A"
			}	
		}	
	}
	
	
# ----- descriptive parameters ----- #	
	ntc <- nlevels(dset[,nv])
	ntc0 <- ntc
	#nlvl <- sapply(dset[,(sapply(dset,class)=="factor")], nlevels)
	nlvl <- sapply(dset[,-ncol(dset)], nlevels)
	col.nlvl <- nlvl[which(col.vars)]
	row.nlvl <- nlvl[which(!col.vars)]								

	nc <- prod(col.nlvl)
	nr <- prod(row.nlvl)

	ncv <- length(col.nlvl)
	nrv <- length(row.nlvl)




# ----- computation of the underlying (relative) frequencies ----- #		
	tt1 <- ftable(tapply(dset$Freq,as.list(dset[,1:nv]),sum),col.vars=which(col.vars))
	tt2 <- spread(ftable(tapply(dset$Freq,as.list(dset[,1:(nv-1)]),sum),col.vars=which(col.vars[1:(nv-1)])),ncol=ntc)

	tt1[is.na(tt1)] <- 0 
	tt2[is.na(tt2)] <- 1 

	H0 <- tt1/tt2
	H <- H0

# ----- a few more auxiliary variables ----- #	

	row.gap.prop <- ifelse(nr > 1, gap.prop * min(1,nrv/ncv), 0) 
	col.gap.prop <- ifelse(nc > ntc, gap.prop * min(1,ncv/nrv), 0) 
	
	rind <- which( (!col.vars)*label > 0)
	cind <- which(col.vars*label > 0)
	nrl <- length(rind)*lab
	ncl <- length(cind)*lab

# ----- model computation in expected mode ----- #	
	if(!is.null(expected) ){
		
		int.dset <- lapply(dset,as.integer)
		dset <- dset[do.call("order",int.dset),]
				
		if(!(mod.type %in% c("poisson","polr", "gen"))){
			stop(simpleError("Wrong mod.type specification!"))
		}
		
				
				
		# ------ logit response residuals ----- #
		if(mod.type == "polr"){ 
			if(resid.type != "response"){
				print(simpleWarning("Argument resid.type ignored. Only response residuals are implemented for polr models."))
				resid.type<-"response"	
			}
			nameslist <- names(dset)[1:(nv-1)]		
			nameslist <- nameslist[which(  !(nameslist %in% c("probability","distribution"))  )]
			single.terms <- paste( nameslist ,collapse = "+")
			interaction.terms <- do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(nameslist[x]),sep="*"))),sep="+"))
			
			full.terms <- paste(nameslist,collapse = "*")
			mod.formula <- as.formula(paste(fterms[nv0],"~",single.terms,"+",interaction.terms,sep=""))
			
			modS.formula <- as.formula(paste(fterms[nv0],"~",full.terms,sep=""))	
			modS <- polr( formula = as.formula(modS.formula),data = dset, weights = dset$Freq, method ="logistic")
			mod <- polr( formula = as.formula(mod.formula),data = dset, weights = dset$Freq, method ="logistic")
			pred <- predict(mod,newdata = subtable(dset,c(1:(nv-1)),keep.zero=F),type="probs")
			fit.values <- as.vector(t(pred))
			tt1r <- ftable(tapply(fit.values,as.list(dset[,1:nv]),sum),col.vars=which(col.vars))
			resid.mat <- (H-tt1r)
			p.value<-anova(mod,modS)[2,7]
			
			if( use.expected.values ){
				H0 <- tt1r/tt2
				H00 <- H
				H <- H0
			}else{
				H00 <- tt1r/tt2
				# H <- H
			}
			
		}
		# ----- poisson model ----- #	
		if(mod.type == "poisson"){
			nameslist <- names(dset)[1:(nv)]		
			nameslist <- nameslist[which(  !(nameslist %in% c("probability","distribution"))  )]
			single.terms <- do.call("paste",c(as.list( nameslist ),sep="+"))
			interaction.terms <- do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(nameslist[x]),sep="*"))),sep="+"))
						
			mod.formula <- paste("Freq~",single.terms,"+",interaction.terms,sep="")
			
			mod <- glm( as.formula(mod.formula),data = dset,  family = "poisson")
		    mod.residuals <- residuals(mod,type=resid.type)
			#resid.mat <-  ftable(tapply(mod.residuals,as.list(dset[,1:nv]),sum),col.vars=which(col.vars))
			resid.mat <- ftable(xtabs(mod.residuals~.,data=dset[,1:nv]),col.vars=which(col.vars))
			Df  <- mod$df.residual
			Dev <- mod$deviance
			p.value <- 1-pchisq(Dev,Df)
			
			pred <- predict(mod,type="response")
			#pred.mat <-  ftable(tapply(pred,as.list(dset[,1:nv]),sum),col.vars=which(col.vars))
			pred.mat <-  ftable(xtabs(pred~.,data=dset[,1:nv]),col.vars=which(col.vars))

			if( use.expected.values ){
				H0 <- pred.mat/tt2
				H00 <- H
				H <- H0
			}else{
				H00 <-  pred.mat/tt2
				
			}
		}
		if(mod.type == "gen"){
			
			p.value <- 1
			
			# pred is already defined by expected = c(...)
			pred.mat <- ftable( xtabs(pred~., data = data[,1:nv]), col.vars = which(col.vars))
			
		if(is.null(residuals)){	
				mod.residuals <- dset$Freq - pred
				resid.type <- "response"
			}else{
				stopifnot(length(residuals) == length(pred))
				mod.residuals <- residuals
			}

			resid.mat <- ftable( xtabs(mod.residuals~., data = data[,1:nv]), col.vars = which(col.vars))

			if( use.expected.values ){
				H0 <- pred.mat/tt2
				H00 <- H
				H <- H0
			}else{
				H00 <-  pred.mat/tt2
			}
		}
		
		
	}else{
		disp.res <- FALSE
		H00 <- H0
	}

# ----- modify parameters for rmbplot mode ----- #	
	if( spine | circular ){
		nlvl[nv] <- 1
		nc <- nc/ntc
		ntc <- 1
	}
# preparing the transformation parameters
if( is.list(freq.trans) ){
	stopifnot(freq.trans[[1]] == "sqrt")
	stopifnot(freq.trans[[2]] > 0)
		tfreq <- "sqrt"
		tfn <- freq.trans[[2]]
	
}else{
	tfreq <- freq.trans
	check <- FALSE
	if(is.null(tfreq)){
		tfreq <-"sqrt"
		tfn <- 1
		check <- TRUE
	}else{
		if( tfreq == "sqrt" ){
			tfn <- 2
			check <- TRUE
			}
		if( tfreq == "log" ){
			tfn <- exp(1)
			check <- TRUE
		}
	}
	if(!check){
		stop(simpleError("Wrong freq.trans specification!"))	
	}
}
# ----- computation of rectangle widths, heights and x/y coordinates in matrix form ----- #	

	W <- matrix(1,ncol=nc,nrow=nr)
	S <- space(nlvl,col.vars,gap.prop=gap.prop,gap.mult=gap.mult,last.col.zero = T,last.row.zero = F)

	x <- ( seq(0,1-1/nc,1/nc)   )*(1-col.gap.prop)
	if(nc > nlevels(dset[,nv])){
		x <- x + c(0,cumsum(S[[1]]))*col.gap.prop/sum(S[[1]])
	}
	if( suppressWarnings(any(S[[2]])) ){
		y <- rev((seq(0,1-1/nr,1/nr))*(1-row.gap.prop) + c(0,cumsum(S[[2]]))*row.gap.prop/sum(S[[2]])) #rev test
	}else{
		y <- 0	
	}
		
	X <- spread(t(matrix( x )),nrow=nr)
	Y0 <- spread(matrix( y ),ncol=nc)
	Y <- Y0
	X2 <- X[,seq(1,nc,ntc)]
	Y2 <- Y[,seq(1,nc,ntc)]


	H2 <- matrix(1,ncol=nc/ntc,nrow=nr)
	tt3 <- ftable(tapply(dset$Freq,as.list(dset[,1:(nv-1)]),sum),col.vars=which(col.vars[1:(nv-1)]))
	
	tt3[is.na(tt3)] <- 0
	W2 <- tt3 / max(tt3)

	if(nc > nlevels(dset[,nv])){
		width.cor <- (ntc-1)*S[[1]][1]/sum(S[[1]])*col.gap.prop
	}else{
		width.cor <- 0
	}
# ----- ---------------- ----- #	
# ----- plotting section ----- #	
# ----- ---------------- ----- #

# ----- opening the basic viewport ----- #	
	
	s0 <- 0.02 # border
	s1 <- 0.04 # yaxis
	s2 <- 0.06 # labs
	s3 <- 0.1  # expected scale
	#dev.new()
	grid.rect(gp=gpar(fill="white",col="white"))
	#if(exists("rmb.id",.GlobalEnv)){
		#popViewport()
		#grid.rect(gp=gpar(fill="white"))
	
	#}else{
	#	.GlobalEnv$rmb.id <- rpois(1,lambda=10000)	
	#}
	vp0 <- viewport(x = 0.5, y = 0.5, w = 1- 2*s0 , h = 1 - 2*s0, just = "centre", name = "vp0")
	pushViewport(vp0)
	vp1 <- viewport(x = 1-yaxis*s1 - s3*(!is.null(expected) ) , y = 0, w = 1- nrl*s2 - yaxis*2*s1 - s3*(!is.null(expected)), h = 1-ncl*s2, just = c("right", "bottom"), name = "vp1")
	pushViewport(vp1)
	# alpha values on a log-scale
	#alpha <- spread(W2,ncol=ntc)*base.alpha
	
	#alpha <- spread(exp(2*(W2-1)),ncol=ntc)*base.alpha
	alpha <- exp(2*(W2-1))
	alpha[alpha < min.alpha] <- min.alpha
	if(!is.null(expected)){
		alpha[alpha < 1] <- 1	
	}

# ----- preparing the color matrix ----- #	
	

			if(is.null(expected) | res.val.only ){
				
				
				
				#has any color rule been defined?
				if(any(c("hsv","hcl","rgb","seq","sequential","sqn","sqt","div","diverging","diverge") %in% col)){
				num.col <- ntc0^(spine | !eqwidth | (eqwidth & circular))
				if( col %in% c("hsv","rgb") ){
					col.def <- formals(rainbow)
					if(!("s" %in% labels(col.opt))){
						col.opt$s <- eval(col.def$s)
					}
					if(!("v" %in% labels(col.opt))){
						col.opt$v <- eval(col.def$v)
					}
					if(!("start" %in% labels(col.opt))){
						col.opt$start <- eval(col.def$start)
					}
					if(!("end" %in% labels(col.opt))){
						col.opt$end <- max(num.col-1,1)/num.col
					}
					if(!("alpha" %in% labels(col.opt))){
						col.opt$alpha <- eval(col.def$alpha)
					}
					colv <- rainbow(num.col,s = col.opt$s, v = col.opt$v, start = col.opt$start, end = col.opt$end, alpha = col.opt$alpha)
					}
				if( col == "hcl" ){
					col.def <- formals(rainbow_hcl)
					if(!("c" %in% labels(col.opt))){
						col.opt$c <- eval(col.def$c)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l <- eval(col.def$l)
					}
					if(!("start" %in% labels(col.opt))){
						col.opt$start <- eval(col.def$start)
					}
					if(!("end" %in% labels(col.opt))){
						col.opt$end <- 360 * (num.col - 1)/num.col
					}
				colv <- rainbow_hcl(num.col,c = col.opt$c, l = col.opt$l, start = col.opt$start, end = col.opt$end)
				}
				if( col %in% c("seq","sqt","sqn","sequential") ){
					col.def <- formals(sequential_hcl)
					if(!("h" %in% labels(col.opt))){
						col.opt$h <- eval(col.def$h)
					}
					if(!("c" %in% labels(col.opt))){
						col.opt$c <- eval(col.def$c.)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l <- eval(col.def$l)
					}
					if(!("power" %in% labels(col.opt))){
						col.opt$power <- eval(col.def$power)
					}
					colv <- rev(sequential_hcl(num.col,h = col.opt$h, c. = col.opt$c, l = col.opt$l, power = col.opt$power))
				}
				if( col %in% c("div","diverging","diverge") ){
					col.def <- formals(diverge_hcl)
					if(!("h" %in% labels(col.opt))){
						col.opt$h <- eval(col.def$h)
					}
					if(!("c" %in% labels(col.opt))){
						col.opt$c <- eval(col.def$c)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l <- eval(col.def$l)
					}
					if(!("power" %in% labels(col.opt))){
						col.opt$power <- eval(col.def$power)
					}
					colv <- diverge_hcl(num.col,h = col.opt$h, c = col.opt$c, l = col.opt$l, power = col.opt$power)
					
				}
					#if((!spine) & eqwidth){
					#	colv <- rep(rainbow(1,alpha=col$alpha),ntc)
					#}
				}else{
					colv <- col[1:ntc0]	
				}					
				col.Mat <- spread( t(rep(colv , nc/ntc)),nrow=nr)
				
				if(eqwidth & !res.val.only){
					
					dm <- dim(col.Mat)
					alpha <- spread(alpha, ncol = num.col)
					col.Mat <- mapply( function(x,y){
						z <- col2rgb(x)/255
						rgb( red = z[[1]], green = z[[2]], blue = z[[3]], alpha = y )
					}, x = col.Mat, y = alpha)
					
					dim(col.Mat) <- dm
					col.Mat <- spread(col.Mat, ncol = ntc0/num.col)
				}			
				
				col.Mat2 <- col.Mat
				
			}else{
				resid.mat[is.nan(resid.mat)] <- 0
				resid.mat[is.na(resid.mat)] <- 0
				if(is.null(resid.max)){
					resid.max <- max(abs(resid.mat),na.rm=T)	
				}
				resid.mat[resid.mat < -resid.max] <- -resid.max
				resid.mat[resid.mat > resid.max]  <-  resid.max
				rmax <- ifelse(abs(resid.max) < 10^{-6},1,resid.max)
				
				if( !use.expected.values ){
					resid.mat <- -resid.mat	
				}
				
				alf <- 0.05
				
				col.Mat <- apply(resid.mat,2,function(x){
						
						 if(resid.type == "pearson"){
							qx <- abs(x)
							
							qx = 2*(1 - sapply( qx, pnorm ))
							#print(qx)
							qx[qx < alf] <- alf
							
							val <- log( qx )/ log( alf )
							#print(val)
							return( hsv( h = 2/3*(sign(x) > 0), s = val, v = 0.7 ) )
						
						}else{
							qx <- abs(x/rmax)
							#rgb(sign(x) < 0,0,sign(x) > 0,alpha = exp(2*(qx-1))*base.alpha ) #qx # cut.rs/2
							qx <- qx - qx %% 1/(cut.rs-1)
							alf <- (exp(qx - 1)-exp(-1))/(1-exp(-1))
							return( rgb(sign(x) < 0,0,sign(x) > 0,alpha = alf*base.alpha ) )#qx # cut.rs/2
						}
					
					})
				colv.grey <- hsv(h=1, s = 0, v = seq(0.9,0.6,-0.3/(ntc0-1)))#spine, eqwidth, circular?
				grey.Mat <- spread( t(rep(colv.grey , nc/ntc)),nrow=nr)
				
#				dev.new()
#				plot(1, xlim=c(0,nrow(col.Mat)), ylim = c(0,ncol(col.Mat)))
#				
#				for(i in 1:nrow(col.Mat)){
#					for(j in 1:ncol(col.Mat)){
#						rect(i-1,j-1,i,j,col = col.Mat[i,j])
#					}
#				}
#				dev.new()
				if( both ){
					col.Mat2 <- col.Mat
					col.Mat3 <- grey.Mat
					col.Mat4 <- grey.Mat
					col.Mat[resid.mat >= 0] <- hsv(0,0,0,alpha=0) # red remains
					col.Mat2[resid.mat <= 0] <- hsv(0,0,0,alpha=0) # blue remains
					col.Mat3[resid.mat < 0] <- hsv(0,0,0,alpha=0) # grey for nonreds remains
					col.Mat4[resid.mat > 0] <- hsv(0,0,0,alpha=0) # grey for nonblues remains
				}else{
					col.Mat2 <- col.Mat	
				}
				
			}
				
# ----- plotting itself ----- #
# 	>>> divided into eqwidth <- T | eqwidth <- F >>> tfreq <- Id/sqrt/log
# 	>>> in spine mode, the target categories will be plotted one after another (see the for-loop)
# 	>>> the 3 diff. draw-function add the rectangles for the relative frequencies, the weights and the background respectively

if( eqwidth ){
	W2trans <- W2
	W3 <- W
	X3 <- X
}else{
	if( tfreq == "log" ){
		W2trans <- log(tt3+1)/max(log(tt3+1))
	}
	if( tfreq == "sqrt" ){
		W2trans <- tt3^(1/tfn) / max(  tt3^(1/tfn) )
	}
	W3 <- spread(W2trans,ncol=ntc)
	X3 <- spread(X[,seq(1,nc,ntc),drop=FALSE],ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.gap.prop) 
}

if( circular ){
		
		P1 = 0*X
		X3 = X2 + 1/ncol(W)/2*(1-col.gap.prop)
		Y3 = Y2 + 1/nrow(W)/2*(1-row.gap.prop)
		
		W2sqrt = tt3^(1/2) / max(  tt3^(1/2) )
		W3 = spread(W2sqrt,ncol=ntc)
		
		
		
	
	if( disp.res ){ 
		H0 <- apply(H0,1:2, function(z) max(0.001,z))
		ratMat = sqrt(H00/H0)
		
		if( any( ratMat[H0 > 0.001] > max.rat ) ){
			warning(cat("residual with ratio", max(ratMat[H0 > 0.001],na.rm=TRUE), " > ",max.rat," occurred!"))	
		}
		ratMat <- apply(ratMat,1:2, function(z) min(max.rat, z))
		
		if(eqwidth){
			R1 = matrix(1,ncol = ncol(W3)*ntc0, nrow = nrow(W3))/ncol(W3)/2*(1-gap.prop) # 0?
			mv <- 1
		}else{
			R1 = spread(W3/ncol(W3)/2*(1-gap.prop),ncol=ntc0)
			mv <- 1/ncol(W3)/2*(1-gap.prop)
		}
			
			R2 <- R1 * ratMat
			rscf <- max(R2[ratMat <= max.rat], na.rm = TRUE) * (ncol(W3)*2/(1-gap.prop))
			if(rscf < 1){
				rscf <- 1	
			}
			
			R2[is.na(R2)] <- 1/ncol(W3)/2*(1-gap.prop)
			R2 <- R2/rscf
			R1 <- R1/rscf
			
			R2 <- apply(R2, 1:2, function(z) min(z, mv))
			R1 <- apply(R1, 1:2, function(z) min(z, mv))
			
		for( i in 1:length(cat.ord) ){
				
				xind = seq( cat.ord[i], nc*ntc0, ntc0 )
				
				P2 = P1 + H0[,xind]

			#exp
			mapply( function(p1,p2,rad,bg,x,y){
				cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*360 ))
					grid.polygon(x = x+c(0,rad*cos(cc*2*pi)),y = y+c(0, rad*sin(cc*2*pi)),
					gp=gpar(fill=bg, lwd=1.5)
					)}, p1 = P1, p2=P2, rad = R1[,xind], bg = col.Mat[,xind], x = X3, y = Y3 )
				
			# obs
			mapply( function(p1,p2,rad,bg,x,y){
					cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*360 ))
					grid.polygon(x = x+c(0,rad*cos(cc*2*pi)),y = y+c(0, rad*sin(cc*2*pi)),
						gp=gpar(fill=bg)
					)}, p1 = P1, p2=P2, rad = R2[,xind], bg = col.Mat2[,xind], x = X3, y = Y3 )
			
			if( both ){
			#exp	
			mapply( function(p1,p2,rad,bg,x,y){
					cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*360 ))
					grid.polygon(x = x+c(0,rad*cos(cc*2*pi)),y = y+c(0, rad*sin(cc*2*pi)),
						gp=gpar(fill=bg, lwd=1.5)
					)}, p1 = P1, p2=P2, rad = R1[,xind], bg = col.Mat3[,xind], x = X3, y = Y3 )
				
			# obs
			mapply( function(p1,p2,rad,bg,x,y){
					cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*360 ))
					grid.polygon(x = x+c(0,rad*cos(cc*2*pi)),y = y+c(0, rad*sin(cc*2*pi)),
						gp=gpar(fill=bg)
					)}, p1 = P1, p2=P2, rad = R2[,xind], bg = col.Mat4[,xind], x = X3, y = Y3 )
				
			}	
				
				P1 = P2
			}
		}else{
			H0 <- apply(H0,1:2, function(z) max(0.0001,z))
			for( i in 1:length(cat.ord) ){
				xind = seq( cat.ord[i], nc*ntc0, ntc0 )
				P2 = P1 + H0[,xind]
				if( eqwidth ){	
					R = H2/nrow(H2)*(1-gap.prop)/2
				}else{
					R = W3/ncol(W3)/2*(1-gap.prop)
				}
				mapply( function(p1,p2,rad,bg,x,y){
					cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*360 ))
					grid.polygon(x = x+c(0,rad*cos(cc*2*pi)),y = y+c(0, rad*sin(cc*2*pi)),
						gp=gpar(fill=bg, lwd=1.5))
					}, p1 = P1, p2=P2, rad = R, bg = col.Mat[,xind], x = X3, y = Y3 )
				P1 = P2
			}
		}
	
}else{ # not circular
	
	if( disp.res ){
		H0 <- apply(H0,1:2, function(z) max(0.001,z))
		ratMat <- sqrt(H00/H0)
		ratMat <- apply(ratMat,1:2, function(z) min(max.rat, z))
	}
	if( disp.res & spine ){	
		if(eqwidth){
			R1 = matrix(1,ncol = ncol(W3)*ntc0/ntc, nrow = nrow(W3))/ncol(W3)/2*(1-gap.prop) # 0?
		}else{
			R1 = spread(W3/ncol(W3)/2*(1-gap.prop),ncol=ntc0/ntc)
		}
			
			R2 <- R1 * ratMat
			rscf <- max(R2[!is.na(R2)]) * (ncol(W3)*2/(1-gap.prop))
			rscf <- max(rscf,1)
	}else{
		rscf <- 1	
	}	
			#R2[is.na(R2)] <- 1/ncol(W3)/2*(1-gap.prop)
			#R2 <- R2/rscf
			#R1 <- R1/rscf


	draw(H2/nrow(H2)*(1-row.gap.prop),( W2trans/ncol(W2trans)*(1-col.gap.prop) + (round(W2trans,digits=7) > 0)*width.cor )/rscf,X2,Y2,alpha = 0.25, bg = "grey")
	draw(H2/nrow(H2)*(1-row.gap.prop),H2/ncol(H2)*(1-col.gap.prop) + width.cor,X2,Y2,alpha = 0.1, bg = "grey")

	if( spine ){
		for( i in 1:length(cat.ord) ){
			if(i > 1){ 
				Y <- Y + H/nrow(H)*(1-row.gap.prop)
			}
			xind <- seq( cat.ord[i], nc*ntc0, ntc0 )
			H <- as.matrix(H0[,xind])
			if( nr < 2 ){
				H <- t(H)	
			}
			draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop)/rscf,X3,Y,alpha = 1, bg = col.Mat[,xind], border = NA)
			if( disp.res ){
				#ratMat = sqrt(H00/H0)
				max.rat = max(ratMat)
						
				W3res = W3/ncol(W3)*(1-col.gap.prop)*ratMat[,xind]
				# obs
				draw(H/nrow(H)*(1-row.gap.prop),W3res/rscf,X3,Y,alpha = 1, bg = col.Mat2[,xind], border = NA)
					
				if( both ){
						#exp	
						draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop)/rscf,X3,Y,alpha = 1, bg = col.Mat3[,xind], border = NA)
						# obs
						draw(H/nrow(H)*(1-row.gap.prop),W3res/rscf,X3,Y,alpha = 1, bg = col.Mat4[,xind], border = NA)
				}	
			}
			draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop)/rscf,X3,Y,alpha = 1, bg = NA, lwd = 1.5)
			
		}
		Y <- Y0
	
	}else{
		
		H <- apply(H,c(1,2),function(x) return(min(x/max.scale,1)))
		
		if( disp.res ){
				draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat, border = NA)
				
				ratMat = sqrt(H00/H0)
				max.rat = max(ratMat)		
				
				W3res = W3/ncol(W3)*(1-col.gap.prop)*ratMat#[,xind]
				# obs
				draw(H00/nrow(H00)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat2, border = NA)
					
				if( both ){
						#exp	
						draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat3, border = NA)
						# obs
						draw(H00/nrow(H00)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat4, border = NA)
				}
				draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = NA, lwd =1.5)	
		}else{
			draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat)	
		}
	}
}
		
	if(yaxis & !circular){
		if(spine){
			rscf <- 1	
		}
		sapply(y,function(z){
			grid.yaxis(at =seq(z,z+1/length(y)*(1-row.gap.prop),1/length(y)*(1-row.gap.prop)/5), label = round(seq(0,max.scale,max.scale/5)/rscf,digits = 3), main = TRUE,
			edits = NULL, name = NULL,
			gp = gpar(), draw = TRUE, vp = NULL)
		})
		sapply(y,function(z){
			grid.yaxis(at =seq(z,z+1/length(y)*(1-row.gap.prop),1/length(y)*(1-row.gap.prop)/5), label = round(seq(0,max.scale,max.scale/5)/rscf,digits = 3), main = FALSE,
			edits = NULL, name = NULL,
			gp = gpar(), draw = TRUE, vp = NULL)
		})
	}
	
	upViewport()
# ----- ---------------- ----- #	
# ----- labeling section ----- #	
# ----- ---------------- ----- #		
	

	
	
	if(lab){
			if(any(col.vars*label)){
	# ----- labeling the x-axis ----- #	

	vpX <- viewport(x = 1-yaxis*s1 - s3*(!is.null(expected)), y = 1-ncl*s2, w = 1-nrl*s2 - yaxis*2*s1 - s3*(!is.null(expected)), h = ncl*s2, just = c("right", "bottom"), name = "vpX")
	pushViewport(vpX)
	cur <- 1
	zc <- cumprod(col.nlvl)[ which(label[which(col.vars)]) ]
    for(i in 1:ncl){
		vpXX <- viewport(x = 1, y = 1 - i*1/ncl, w = 1, h = 1/ncl, just = c("right", "bottom"), name = "vpXX")
		pushViewport(vpXX)
			if( varnames ){
				grid.text(  names(dset)[cind[i]],x = 0.5 , y = 5/6, just = "centre",gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
			}
			#z <- (nlvl[cind])[i]*cur
			z <- zc[i]
			grid.text(  rep( clabs[[i]] ,cur), x = X[1,seq(1,nc,nc/z)] + (X[1,seq(nc/z,nc,nc/z)] - X[1,seq(1,nc,nc/z)] +1/nc*(1-col.gap.prop))/2, y = 1/3, just = "centre",gp=gpar(cex=lab.cex)) 
			if( boxes ){
				grid.rect(  y = 1/3 , x = X[1,seq(1,nc,nc/z)] , just = c("left","centre"), width = X[1,seq(nc/z,nc,nc/z)] - X[1,seq(1,nc,nc/z)] +1/nc*(1-col.gap.prop), height = 1/2, gp = gpar(fill="transparent"))
			}
			cur <- z
		popViewport()
	}
		
	upViewport()
	}
		if(any( (!col.vars)*label)){
# ----- labeling the y-axis ----- #			
	vpY <- viewport(x = nrl*s2, y = 0, w = nrl*s2, h = 1-ncl*s2, just = c("right", "bottom"), name = "vpY")
	pushViewport(vpY)
	cur <- 1
	zr <- cumprod(row.nlvl)[ which(label[which(!col.vars)]) ]
	for(i in 1:nrl){
		
		vpYY <- viewport(x = 0 + (i-1)*1/nrl, y = 0, w = 1/nrl, h = 1, just = c("left", "bottom"), name = "vpYY")
		pushViewport(vpYY)
			if( varnames & (nr > 1) ){
				grid.text(  names(dset)[rind[i]],x = 1/6 , y = 0.5, just = "centre",rot=90,gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
			}			
			#z <- (nlvl[rind])[i]*cur
			z <- zr[i]
			grid.text(  rep( rlabs[[i]] ,cur), y = Y[seq(1,nr,nr/z),1] + (Y[seq(nr/z,nr,nr/z),1] - Y[seq(1,nr,nr/z),1] +1/nr*(1-row.gap.prop))/2, x = 2/3, just = "centre",rot=90,gp=gpar(cex=lab.cex)) 
			if( boxes ){
				grid.rect(  x = 2/3 , y = Y[seq(1,nr,nr/z),1] , just = c("centre","bottom"), height = Y[seq(nr/z,nr,nr/z),1] - Y[seq(1,nr,nr/z),1] +1/nr*(1-row.gap.prop), width = 1/2, gp = gpar(fill="transparent"))
			}
					
			cur <- z
		popViewport()
	}
	upViewport()
	}
	}
	# ----- scale for expected mode ----- #
	if(!is.null(expected) & !res.val.only){

		vpS <- viewport(x = 1, y = 0 , w = s3, h = 0.8, just = c("right","bottom"), name = "vpS")
		pushViewport(vpS)
		
		#alpha values for the residuals on a log-scale
		#alfav <- exp(seq(0,-cut.rs/2+0.5,-0.5))
		if( resid.type == "pearson" ){
			sval <- c(log( (1 - seq(1-alf, alf, -(1-2*alf)/(cut.rs-1) ) ) )/ log( alf ))
			sval <- c(sval,0, rev(sval))
			
			alfy <- qnorm(1-alf/2)
			minv <- min(-alfy,min(resid.mat[resid.mat < 0]))
			maxv <- max(alfy,max(resid.mat[resid.mat > 0]))
			ys <- c( minv, seq(-alfy,alfy,2*alfy/(2*cut.rs-1)), maxv)
			colv <- hsv( c(rep(0,cut.rs),0,rep(2/3,cut.rs)),  sval,  0.8)
			
			grid.rect(  x = 0.1 , y = ((ys - minv)/(maxv - minv))[-length(ys)], just = c("left","bottom"), height = diff(ys)/(maxv-minv), width = 0.3 , gp = gpar(fill=colv))
			grid.text( label=round(ys,digits=2), y = (ys - minv)/(maxv - minv), x = 0.45, just = c("left","centre"),rot=0,gp=gpar(cex=lab.cex))
		
			
		}else{
			alfav <- (exp(seq(0,-1,-1/(cut.rs-1)))-exp(-1))/(1-exp(-1))
			colv <- rgb(rep(c(1,0),each=cut.rs),0,rep(c(0,1),each=cut.rs),alpha = c(alfav,rev(alfav))*base.alpha)
			grid.rect(  x = 0.1 , y = seq(0,1-1/2/cut.rs,1/2/cut.rs) , just = c("left","bottom"), height = 1/2/cut.rs, width = 0.3 , gp = gpar(fill=colv))
			grid.text( label=round(seq(-resid.max,resid.max,resid.max/min(10,cut.rs)),digits=2), y = seq(0,1,1/2/min(10,cut.rs)), x = 0.45, just = c("left","centre"),rot=0,gp=gpar(cex=lab.cex))
		}
		#grid.rect(  x = 0.1 , y = 0, just = c("left","bottom"), height = 1, width = 0.3 , gp = gpar(fill=hsv(0,0,0,alpha=0.1)))
			
		upViewport()	
		vpS2 <- viewport(x = 1, y = 0.8  , w = 0.1, h = 0.2 , just = c("right","bottom"), name = "vpS2")
		pushViewport(vpS2)
		
		grid.text( label=resid.type, y = 0.5, x = 0.2, just = "centre",rot=90,gp=gpar(cex=lab.cex*1.2))
		grid.text( label="residuals", y = 0.5, x = 0.4, just = "centre",rot=90,gp=gpar(cex=lab.cex*1.2))
		grid.text( label=paste("p.value =",format(p.value,digits=2)), y = 0.5, x = 0.7, just = "centre",rot=90,gp=gpar(cex=lab.cex*1.1))
	}
	
	upViewport()
	op <- options()
	options(show.error.messages = FALSE)
	try(upViewport())
	options(op)
	return(invisible(TRUE))
}

# ------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------ #


space <- function(nlvl,col.vars,gap.mult = 1.1,gap.prop = 0.1,last.col.zero = F, last.row.zero = F){
	# 	>>> function to compute the gaps widths for the rmb plot
	# 	>>> might be extended in further versions of the package
	
	# ----- parameter check ----- #
		if( !(  (gap.prop > 0)&(gap.prop < 1) ) ){
			stop(simpleError("Wrong gap.prop specification!"))
		}

	col.nlvl <- nlvl[col.vars]
	row.nlvl <- nlvl[!col.vars]								

	nc <- prod(col.nlvl)
	nr <- prod(row.nlvl)

	ncv <- length(col.nlvl)
	nrv <- length(row.nlvl)

	col.spaces <- rep(gap.prop,nc-1)
	row.spaces <- rep(gap.prop,nr-1)

	if( ncv > 1 ){
		cind <- 1:(nc-1)
		cur <- 1
		for( i in (ncv-1):1 ){ 
			cur <- cur*col.nlvl[i+1]
			curind <- which( (cind %% cur) == 0)
			col.spaces[ curind  ] <- col.spaces[ curind  ] * gap.mult
		}
	}
	if( nrv > 1 ){
		rind <- 1:(nr-1)
		cur<-1
		for( i in (nrv-1):1 ){ 
			cur <- cur*row.nlvl[i+1]
			curind <- which( (rind %% cur) == 0)
			row.spaces[ curind  ] <- row.spaces[ curind  ] * gap.mult
		}
	}
	if(last.col.zero){
		col.spaces[which(col.spaces == gap.prop)] <- 0
		col.spaces <- col.spaces/gap.prop
	}
	if(last.row.zero){
		row.spaces[which(row.spaces == gap.prop)] <- 0
		row.spaces <- row.spaces/gap.prop
	}
	
	return( list(col.spaces,row.spaces))
}

# ------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------ #


spread <- function(M,ncol = 1,nrow = 1){
	# 	>>> this function turns each matrix entry into a nrow x ncol - matrix
	M <- as.matrix(M)
	M2 <- apply(M,2,function(x) rep(x,each = nrow)  )
	if( (nrow(M) == 1)&(nrow == 1) ){
		M2 <- matrix(M2,ncol = ncol(M))
	#dim(M2) <- dim(M)
	}
	
	M3 <- t(apply(t(M2),2,function(x) rep(x,each = ncol)  ))
	if( (ncol(M) == 1)&(ncol == 1) ){
		M3 <- matrix(M3,nrow = nrow(M2))
		#dim(M3) <- dim(M2)
	}
	return(M3)
}

