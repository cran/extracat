
#####################################################
##                    FUNKTION                     ##
#####################################################
scpcp<-function(data, 
                gap = 0.2,
                sort.individual=TRUE,
                level.width=0.2,
                polygon = TRUE,
                base.colour = alpha("black",0.7), # Farbe für Polygone oder Linien falls kein doodle
                lab.opt = list(rot = 0, col = 1, bg = TRUE, abbr = FALSE, abbr.var = 12, hide.sel = TRUE),
                sel =  NULL, # Vektor oder alternativ Beschreibung Ã¡ la "Cont=='Low'&Infl=='Low'"
                sel.hide = TRUE,
                sel.palette = NULL,
                col.opt = list(),
                plot = TRUE,
	            return.coords = !plot)
{
 # alternative fÃ¼r shadowtext?
  
  
  s.old <- Sys.time()
  sv <- 0
  sv.names <- c("init")
  
#####
				#text.background = TRUE,
               # text.rotation = 0, # in Grad
	           # text.abbreviation = TRUE,
               # text.colour = "black", # Text-Hintergrund-Farbe!!!


if( "rot" %in% names(lab.opt) ){
	text.rotation <- lab.opt$rot
}else{
	text.rotation <- 0
}
if( "las" %in% names(lab.opt) ){
	las <- lab.opt$las
}else{
	las <- 1
}

if( "hide.sel" %in% names(lab.opt) ){
	hide.sel <- as.logical(lab.opt$hide.sel)
}else{
	hide.sel <- TRUE
}

if( "abbr.var" %in% names(lab.opt) ){
	abbr.var <- lab.opt$abbr.var
}else{
	abbr.var <- 12
}

if( "bg" %in% names(lab.opt) ){
	text.background <- lab.opt$bg
}else{
	text.background <- TRUE
}
if( "abbr" %in% names(lab.opt) ){
	text.abbreviation <- lab.opt$abbr
}else{
	text.abbreviation <- FALSE
}
if( "col" %in% names(lab.opt) ){
	text.colour <- lab.opt$col
}else{
	text.colour  <- "black"
}
if( "alpha" %in% names(col.opt) ){
	colour.alpha <- col.opt$alpha
	col.opt$alpha <- 1
}else{
	if(!polygon){
		colour.alpha <- min(0.7,700/nrow(data))
	}else{
		colour.alpha <- 1
	}
}

  ##################### FORMATIERUNG  #####################

  # Formatierung, falls die Daten eine Variable "Freq" enth?lt
  if("Freq" %in% names(data)){  
    data <- untableSet(data)
  } else {
  	if(inherits(data,"table")){
  		data <- subtable(data,1:length(dim(data)),allfactor=TRUE)
  		data <- untableSet(data)
  	}else{
  		for(i in 1:ncol(data)){
			data[[i]] <- as.factor(data[[i]])
		}
  	}
     #data <- data.frame(lapply(as.data.frame(data), as.factor)) # Sicherstellen, dass alle Variablen vom Typ Factor sind
  }

  ##################### DOODLE  #####################
  

  # doodle kann entweder durch einen externen Vektor (sel), oder durch direkte ?bergabe der zu highlightenden F?lle im Funktionsaufruf initiiert werden
  # zb bei sel="Cont=='Low'&Infl=='Low'" 
  
  if(!is.null(sel)){
    doodle<-TRUE
    
    if(is.character(sel)){
      #attach(data)
      doodle.v <- with(data, as.factor(eval(parse(text=sel))) ) # direkte Evaluierung der zu highlightenden F?lle
      #detach(data) # attach/detach notwendig, damit die ?bergebenen Variablennamen (zB. Cont) verstanden werden
      
    } else { 
      doodle.v<-as.factor(sel)
    }
    
    s.now <- Sys.time()
    sv <- c(sv,s.now-s.old)
    s.old <- s.now
    sv.names <- c(sv.names,"doodle")  
  } else {
    doodle <- FALSE
  }

  sel.hide <- doodle & sel.hide # hide nur erlauben, falls Ã¼berhaupt gedoodled wird


  ##################### FORMATIERUNG 2  #####################
  
  if(doodle){
    data <- cbind(doodle.v,data)
  }

  # should be the case after untableSet(x) anyway, but in general...
  # untableset() sortiert die Daten von erster zu letzter Variabe, wir brauchen es andersrum
  base.order <- do.call(order,c(data, decreasing = FALSE))
  data <- data[base.order,]
  
  s.now <- Sys.time()
  sv <- c(sv,s.now-s.old)
  s.old <- s.now
  sv.names <- c(sv.names,"format")
  
  ##### NUE?TZLICHES  
  
  N <- nrow(data)
  m <- ncol(data)
  nm <- names(data)
  labels <- lapply(data, levels)
  
  if(text.abbreviation != FALSE){
  	labels.abbr <- lapply(labels, abbreviate, minlength=text.abbreviation)
  }else{
  	labels.abbr <- labels
  }
  


  ##################### BERECHNUNGEN  #####################
  
  # get the sequences ...
  
  sequences <- lapply(data, table)
  

  ##### LINIEN KOORDINATEN

  seq.list <- lapply(sequences, function(z){ # Berechnung aller Koordinaten f?r alle Linien (nicht nicht zusammenh?ngend)
    p <- c(0,cumsum(z/N))*(1-gap) # Proportionen der Kategorien in nicht-gap-Bereich
    k <- length(z) # Anzahl Kategorien
    gap.proportion <- gap/(k-1) # H?he jedes gap-Bereichs
    seqs <- list()
    for(i in 1:k){
      seqs[[i]] <- seq(p[i],p[i+1], (p[i+1]-p[i])/(z[i]-1)) + (i-1)*gap.proportion 
    }
    return(seqs)
  })
  

  ##### LINIEN IDs

  id.list <- mapply( function(y, z){ # Durchnummerierung
    lapply(y, function(w){
      which(z == w)
    })
  },y = labels, z = data, SIMPLIFY=F)
 

  s.now <- Sys.time()
  sv <- c(sv,s.now-s.old)
  s.old <- s.now
  sv.names <- c(sv.names,"seq/id")
  

  ##### LINIEN

  lines.unsorted <- mapply(function(y,z){ # Kombinieren von seq.list und id.list um komplette Linien zu erzeugen
    ret <- rep(0,N) 
    for(i in seq_along(y)){                
      ret[ z[[i]] ] <-   y[[i]]
    }
    return(ret)  
  }, y = seq.list, z = id.list )
  
  s.now <- Sys.time()
  sv <- c(sv,s.now-s.old)
  s.old <- s.now
  sv.names <- c(sv.names,"lines")
  

  ##################### SORTIERUNG  ##################### 
  
  if(sort.individual & !doodle){  # Sortiert jede Sequenz nach den R?ngen der jeweils linken Variable
    lines <- lines.unsorted
    e1 <- environment()
    for(i in 2:m){
      sapply(id.list[[i]], function(z){
        e1$lines[z,i]  <- lines[z,i][rank(lines[z,i-1])] 
        return(invisible(TRUE))
      })
    }  
  } else if(sort.individual & doodle){ #Sortiert nach der jeweils linken Variable UND dem doodle
    lines <- lines.unsorted
    e1 <- environment()
    for(i in 2:m){
      sapply(id.list[[i]], function(z){
        tapply(z,e1$data[z,1],function(y){
        e1$lines[y,i]  <- lines[y,i][rank(e1$lines[y,i-1])]
        })
      return(invisible(TRUE))
      })     
    }
  }  else {
    lines<-lines.unsorted
  }
    
  s.now <- Sys.time()
  sv <- c(sv,s.now-s.old)
  s.old <- s.now
  sv.names <- c(sv.names,"sort")
  
  ##################### PLOT-VORBEREITNG ###################### 
  
    
  # Abschneiden des ergaenzten Datensatzes auf den Usprungsdatensatz zur Darstellung

  if(sel.hide){ # falls gehighlighted oder gefÃ¤rbt wird (und versteckt!)
    doodle.v<-data[,1] # der sortierte Vektor soll weiterbenutzt werden
    data<-data[2:m]
    lines<-lines[,2:m]
    lines.unsorted <- lines.unsorted[,2:m]
    seq.list <- seq.list[2:m]
    id.list <- id.list[2:m]
    labels <- labels[2:m]
    labels.abbr <- labels.abbr[2:m]
    m <- ncol(data)
    nm <- nm[2:length(nm)]
  } else if(doodle){
    doodle.v <- data[,1] # der sortierte Vektor soll weiterbenutzt werden 
  } 

  

    

  ##### HILFSKOORDINATEN

  middles<- as.vector(rapply(seq.list, mean)) 
  
  s.now <- Sys.time()
  sv <- c(sv,s.now-s.old)
  s.old <- s.now
  sv.names <- c(sv.names,"mid") 
         
if(plot){
	
	  # Anpassung der Raender fuer bessere Optik
  if(is.character(sel)){
    par(mar=c(3,0,0,0)) # um nicht zu groÃŸe Margins zu bekommen, mit Platz fÃ¼r Vektor-Beschreibung
  } else {
    par(mar=c(2,0,0,0)) # nicht zu groÃŸe Margins 
  }  
  
  ##### Farbbestimmung 
  if(doodle){
  	if(is.null(sel.palette)){
  		ndv <- length(levels(doodle.v))
  		if(ndv < 9){
  			sel.palette <- 1:ndv
  		}else{
  			sel.palette <- "rgb"
  		}
		colour.alpha <- 0.7
	}
  	
  	
    ##if(length(levels(doodle.v))==2){ # bei binärem Vektor wird von highlighting ausgegangen
     ## doodle.colours <- alpha(1:2, colour.alpha)
    ##} else{
      doodle.colours <- getcolors(length(levels(doodle.v)), sel.palette ,col.opt = col.opt) 
      doodle.colours <- alpha(doodle.colours, colour.alpha)
   ## }
  } else {
    doodle.colours <- base.colour
  }
  
  	if( "border" %in% names(col.opt) ){
		rect.border <- rep(col.opt$border,length(doodle.colours))[1:length(doodle.colours)]
	}else{
		rect.border  <- doodle.colours
	}
  
  
  ##################### LINE-PLOT ###################### 
  
  if(!polygon){
    
    lines.doubled <- lines[,rep(1:m,each=2)] # Achsen verdoppeln um Platz f?r Beschriftung zu schaffen
    
    # Bestimmung der Koordinaten f?r die Achsen bei ?bergebener level.width
    xcoords <- c(1,1+level.width)
    for (i in 2:m){
      xcoords <- c(xcoords, c(i,i+level.width))
    }
    
   ##  Zeichung aller Linien gruppiert (und eingef?rbt) nach highlighting oder colouring
  #  plot(1, xlim=c(1,m+level.width),ylim=c(0,1), axes=F, xlab=NA, ylab=NA, type="l")     
   # for(j in 1:length(levels(doodle.v))){
   #   apply(lines.doubled[which(doodle.v==levels(doodle.v)[j]),],1,function(z) lines(xcoords,z, col=doodle.colours[j]))
   # }
 
 
         plot(1, xlim=c(1,m+level.width),ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA,panel.first ={
         	for(j in 1:length(levels(doodle.v))){
      			apply(lines.doubled[which(doodle.v==levels(doodle.v)[j]),],1,function(z){
      				lines(xcoords,z,col=doodle.colours[j])
      			}) 
    		}
         } , type="l") 
        # plot(1, xlim=c(1,m+level.width),ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA,panel.first ={
        # 	by(cbind(doodle.colours[as.integer(doodle.v)],lines.doubled),doodle.v,function(y){
        # 		apply(y,1,function(z){
      	#			lines(xcoords,z[-1],col=alpha(z[1],min(0.7,700/nrow(data))))
      	#		}) 
        # 	})
        # } , type="l") 
       
  } else {
    
    ##################### POLYGON-PLOT ######################
    
    if(polygon){     
      
      ##### POLYGONE
      
      # immer gleicher Aufbau:
      # Ausgangspunkt finden
      # Matrix mit allen Eckpunkten der Polygone erstellen
      # per apply fÃ¼r jedes Eckpunkte-Quadrupel ein Polygon zeichnen
      
      as.pol<- function(color = 1){
        cc<-as.factor(paste(data[,1])) # Polygone fÃ¼r die erste Achse
        M2 <- cbind(tapply(lines[,1],cc,min),tapply(lines[,1],cc,max),tapply(lines[,1],cc,max),tapply(lines[,1],cc,min)) 
        apply(M2,1,function(z){
          polygon(x= c(1,1,1+level.width,1+level.width), y = z, col = doodle.colours,border=rect.border) 
          return(invisible(TRUE))
        })

        for(i in 2:m){
          cc<-as.factor(paste(cc,data[,i])) # achsenverbindenden Polygone
          M <- cbind(
          		tapply(lines[,i-1],cc,min),
          		tapply(lines[,i-1],cc,max),
          		tapply(lines[,i],cc,max),
          		tapply(lines[,i],cc,min))
          		
          apply(M,1,function(z){
            polygon(x= c(i-1+level.width,i-1+level.width,i,i), y = z, col = doodle.colours,border=rect.border)  
            return(invisible(TRUE))
          })
          
          M1<-M[,c(3,4,4,3), drop = FALSE] # Polygone fÃ¼r alle weiteren Achsen
          apply(M1,1,function(z){
            polygon(x= c(i,i,i+level.width,i+level.width), y = z, col = doodle.colours ,border=rect.border)
            return(invisible(TRUE))
          })
        }       
      }

      doodling <- function(){
        
        ##### POLYGONE        
        # immer gleicher Aufbau:
        # Ausgangspunkt finden
        # Matrix mit allen Eckpunkten der Polygone erstellen
        # per apply fÃ¼r jedes Eckpunkte-Quadrupel ein Polygon zeichnen
        
        # fÃ¼r jedes Level des zu betrachtenden Vektors extra
        
        for(j in 1:length(levels(doodle.v))){
          cc<-as.factor(paste(data[which(doodle.v==levels(doodle.v)[j]),1]))
          M2 <- cbind(tapply(lines[which(doodle.v==levels(doodle.v)[j]),1],cc,min),
          tapply(lines[which(doodle.v==levels(doodle.v)[j]),1],cc,max),
          tapply(lines[which(doodle.v==levels(doodle.v)[j]),1],cc,max),
          tapply(lines[which(doodle.v==levels(doodle.v)[j]),1],cc,min))
          apply(M2,1,function(z){
            polygon(x= c(1,1,1+level.width,1+level.width), y = z, col = doodle.colours[j],border=rect.border[j])
            return(invisible(TRUE))
          })
          for(i in 2:m){
            cc<-as.factor(paste(cc,data[which(doodle.v==levels(doodle.v)[j]),i]))
            M <- cbind(tapply(lines[which(doodle.v==levels(doodle.v)[j]),i-1],cc,min),
            tapply(lines[which(doodle.v==levels(doodle.v)[j]),i-1],cc,max),
            tapply(lines[which(doodle.v==levels(doodle.v)[j]),i],cc,max),
            tapply(lines[which(doodle.v==levels(doodle.v)[j]),i],cc,min))
            
            apply(M,1,function(z){
              polygon(x= c(i-1+level.width,i-1+level.width,i,i), y = z, col = doodle.colours[j],border=rect.border[j])  
              return(invisible(TRUE))
            })
            
            M1<-M[,c(3,4,4,3), drop = FALSE]
            apply(M1,1,function(z){
              polygon(x= c(i,i,i+level.width,i+level.width), y = z, col = doodle.colours[j],border=rect.border[j])
              return(invisible(TRUE))
            })
          }
        }
      }
      
      
      ##### Eigentiches Zeichnen des Plots
      
      if(!doodle){
        plot(1, xlim=c(1,m+level.width),ylim=c(0,1), axes=F, xlab=NA, ylab=NA,panel.first = as.pol(), type="l")
      } else {
        plot(1, xlim=c(1,m+level.width),ylim=c(0,1), axes=F, xlab=NA, ylab=NA,panel.first = doodling(), type="l")
      }
    }
  }  
  
  s.now <- Sys.time()
  sv <- c(sv,s.now-s.old)
  s.old <- s.now
  sv.names <- c(sv.names,"plot")
  

  ##################### LEVELS & LABELS ##################### 
  
  middlesX <- list()
  for(i in 1:ncol(lines.unsorted)){
    middlesX[[i]] <-  tapply(lines.unsorted[,i],data[,i],mean)
  }
  
  if (text.abbreviation){
    ll <- labels.abbr
  } else {
    ll <- labels
  }

  #xx <- as.list( 0.5*level.width + seq_along(middlesX))
 # mapply(function(x,y,z){
 #     if(text.background){
 #       shadowtext(rep(x,length(y)), y, z, col="gray90", font=1, bg=text.colour, srt=text.rotation)
 #     } else {
 #       text(x, y, z, col=1, font=2, srt = text.rotation)
 #       text(x, y, z, col="gray90", font=1, srt = text.rotation)
 #     } 
 #     return(invisible(TRUE)) 
 #   }, x = xx, y = middlesX, z = l, SIMPLIFY=F)  
xx <- 0.5*level.width + seq_along(middlesX)
xx <- rep(xx,sapply(middlesX,length))
yy <- unlist(middlesX)
zz <- unlist(ll)
 

   if(text.background){
        bgtext(xx, yy, zz, col="gray90", font=1, bg=text.colour, srt=text.rotation)
     } else {
        text(xx, yy, zz, col=1, font=2, srt = text.rotation)
        text(xx, yy, zz, col="gray90", font=1, srt = text.rotation)
     } 
    
    s.now <- Sys.time()
    sv <- c(sv,s.now-s.old)
    s.old <- s.now
    sv.names <- c(sv.names,"writelvl")

  
  ##### PLOT-UNTERSCHRIFTEN
  if(abbr.var != FALSE){
  	nm <- sapply(nm, abbreviate, minlength=as.integer(abbr.var))
  }
  
  
  mtext(nm, side=1, line=0, at=(1:m)+0.5*level.width, font=2, las = las) # Variablennamen unter die jeweiligen Achsen
  if(is.character(sel) & !hide.sel){ # falls vorhanden, Beschreibung des Highlightings
    mtext(paste("Highlight:",sel), side=1, line=1, at=((m+1+level.width)/2), font=1, cex=0.8, col= "red")  
  }  

}#end if plot

  ##################### RETURN #####################  
  
  if(return.coords){
  	colnames(lines) <- paste("y",colnames(lines),sep=".")
  	ret <- as.data.frame(cbind(lines,data))
  	return(invisible(ret))
  }else{
  	return(invisible(TRUE))
  }
  

  ##################### LAUFZEITEN ##################### 
  
  sv<-c(sv,sum(sv))
  sv.names<- c(sv.names,"sum")
  names(sv) <- sv.names
  print(round(sv,3))

}

bgtext <- function(x, y, labels, col='white', bg='black',
	k = 8, r=0.1, ... ) {

	theta= seq(pi/4, 2*pi, length.out=k)
	#xy <- xy.coords(x,y)
	xr <- r*strwidth('A')
	yr <- r*strheight('A')

	text( rep(x,each = k+1) + c(cos(theta)*xr,0), rep(y,each = k+1) + c(sin(theta)*yr,0), 
		rep(labels,each=k+1), col= c(rep(bg,k),col), ... )

	#text(x, y, labels, col=col, ... ) 
	return(invisible(TRUE))
}

