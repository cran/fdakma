kma.similarity <-
function(
  x.f = NULL, y0.f = NULL, y1.f = NULL,
  x.g = NULL, y0.g = NULL, y1.g = NULL, 
  similarity.method, unif.grid = TRUE){
  
  similarity.method.available <- c('d0.pearson', 'd1.pearson', 
                                   'd0.L2', 'd1.L2', 'd0.L2.centered','d1.L2.centered')
  
  # Check if similarity.method has a feasible value
  if (!similarity.method %in% similarity.method.available)
  {
    stop('Value of "similarity.method" not valid. Possibles choices are: ',
         '"', similarity.method.available[1], '"', ' ',
         '"', similarity.method.available[2], '"', ' ',
         '"', similarity.method.available[3], '"', ' ',
         '"', similarity.method.available[4], '"', ' ',
         '"', similarity.method.available[5], '"', ' ',
         '"', similarity.method.available[6], '"'
    )
  }
  
  if (similarity.method=='d0.pearson' 
      || similarity.method=='d0.L2'
      || similarity.method=='d0.L2.centered')
  {work.with.deriv.sim = 0}
  if (similarity.method=='d1.pearson' 
      || similarity.method=='d1.L2'
      || similarity.method=='d1.L2.centered')
  {work.with.deriv.sim = 1}
  
  if (!unif.grid %in% c(TRUE,FALSE))
  { 
  	stop('Value of parameter "unif.grid" not valid. 
         It must be TRUE or FALSE') 
  }
  ##############################################
  
  var.temp <- FALSE
  var.temp2 <- FALSE
  
  MAX_NUMBER <- 10^(150)
  
  if ((length(x.f)==0 && length(y0.f)==0)
      || (length(x.g)==0 && length(y0.g)==0))
  {
    var.temp <- TRUE
    var.temp2 <- TRUE
    warning('Domains of two curves in kma.similarity function do not overlap. If using kma function of fdakma package, suggested reduction of "t.max" and/or "m.max" or set "fence" equal to TRUE')
     if (similarity.method=='d1.pearson' || similarity.method=='d0.pearson')
     {
     	OUT <- 0
     }
     if (similarity.method=='d0.L2' || similarity.method=='d1.L2' || similarity.method=='d0.L2.centered' || similarity.method=='d1.L2.centered')
     {
     	OUT <- MAX_NUMBER
     }
  }
  
  if (var.temp==FALSE)
  {
    ##############################################
    # Control of value of x0/1 and y0/1 of f and g
    if (work.with.deriv.sim==0)
    {
      if (length(x.f)==0 || length(y0.f)==0
          || length(x.g)==0 || length(y0.g)==0)
      {
        stop('You have to provide the abscissas (x.f and x.g) and the evaluations of functions (y0.f and y0.g) with the similarity.method chosen (',similarity.method,')')
      }
      if (length(y1.f)!=0
          || length(y1.g)!=0)
      {
        warning('You provided the evaluations of function first dervivatives (y1.f and/or y1.g) with the chosen similarity.method (',similarity.method,').
                These values has been ignored.')
      }
      
    }
    if (work.with.deriv.sim==1)
    {
      if (length(x.f)==0 || length(y1.f)==0
          || length(x.g)==0 || length(y1.g)==0)
      {
        stop('You have to provide the abscissa (x.f and x.g) and the evaluations of function first derivatives (y1.f and y1.g) with the chosen similarity.method (',similarity.method,')')
      }
      
      if (length(y0.f)!=0
          || length(y0.g)!=0)
      {
        warning('You provided the evaluations of functions (y0.f and/or y0.g) with the chosen similarity.method (',similarity.method,').
                These values has been ignored.')
      }
      
    }
    ##############################################
    # Assign the variables with which we will work
    if (work.with.deriv.sim==0)
    {
      y.f <- y0.f
      y.g <- y0.g
    }
    if (work.with.deriv.sim==1)
    {
      y.f <- y1.f
      y.g <- y1.g
    }
    
    
    
    #  Transform/copy functions f and g in arrays to work uniquely
    if( length(dim(y.f))==2 ){
      
      if (dim(y.f)[1]==1){
    	y.f <- t(y.f)
    	y.g <- t(y.g)
   	 }
      r.t <- dim(y.f)[2]
      n.camp.f <- dim(y.f)[1]
      n.camp.g <- dim(y.g)[1]
      
      function1 <- array(0,dim=c(1,n.camp.f,r.t))
      function1[, , 1:r.t] <- y.f
      function2 <- array(0,dim=c(1,n.camp.g,r.t))
      function2[, , 1:r.t] <- y.g
      
    }else{
      r.t <- 1
      if (class(y.f)=='numeric') 
      {
        n.camp.f <- length(y.f)
      }else{n.camp.f <- dim(y.f)[1]}  
      
      if (class(y.g)=='numeric') 
      {
        n.camp.g <- length(y.g)
      }else{n.camp.g <- dim(y.g)[1]}    
      
      function1 <- array(0,dim=c(1,n.camp.f,r.t))
      function1[,,1:r.t] <- y.f
      function2 <- array(0,dim=c(1,n.camp.g,r.t))
      function2[,,1:r.t] <- y.g
    }  
    
    # Find the common domain between the two functions
    
    if (unif.grid==TRUE)
    {
      # In this case we generate a uniform grid in the 
      # common domain fo the two functions
      estr.inf.x <- max( min(x.f), min(x.g))
      estr.sup.x <- min( max(x.f), max(x.g))  
      precisione.f <- length(intersect(which( x.f>=estr.inf.x ), which( x.f<=estr.sup.x )))
      precisione.g <- length(intersect(which( x.g>=estr.inf.x ), which( x.g<=estr.sup.x )))
      n.out.sim <- max(precisione.f,precisione.g)
      x.com.sim <- seq(estr.inf.x,estr.sup.x, length=n.out.sim)
      
      if (length(x.com.sim)<=1)
      { 
        var.temp <- TRUE
        var.temp2 <- TRUE
        warning('Domains of two curves in kma.similarity function do not overlap. If using kma function of fdakma package, suggested reduction of "t.max" and/or "m.max" or set "fence" equal to TRUE')
    	if (similarity.method=='d1.pearson' || similarity.method=='d0.pearson')
     	{
     		OUT <- 0
     	}
     	if (similarity.method=='d0.L2' || similarity.method=='d1.L2' || similarity.method=='d0.L2.centered' ||similarity.method=='d1.L2.centered')
     	{
     		OUT <- MAX_NUMBER
     	}
      }
      
      if (var.temp2==FALSE)
      {
        # Approximate the 2 functions on the grid x.com.sim
        y.f.appr <- array(NA,dim=c(1,length(x.com.sim),r.t))
        y.g.appr <- array(NA,dim=c(1,length(x.com.sim),r.t))
        for(l in 1:r.t)
        {
          y.f.appr[1,,l] <- approx(x.f,function1[1,,l],xout=x.com.sim)$y
          y.g.appr[1,,l] <- approx(x.g,function2[1,,l],xout=x.com.sim)$y
        }
      }
    }else{   
        # In this case the two functions are supposed to have 
        # both equal and uniform abscissas, fine enough to 
        # compute the similarity measure
        x.com.sim <- x.f  
        if (length(x.com.sim)<=1)
        { 
     	    warning('Domains of two curves in kma.similarity function do not overlap. If using kma function of fdakma package, suggested reduction of "t.max" and/or "m.max" or set "fence" equal to TRUE')
    	 	if (similarity.method=='d1.pearson' || similarity.method=='d0.pearson')
     		{
     			OUT <- 0
    		}
     		if (similarity.method=='d0.L2' || similarity.method=='d1.L2' || similarity.method=='d0.L2.centered' || similarity.method=='d1.L2.centered')
    		{
     			OUT <- MAX_NUMBER
     		}
        }
      
        if(var.temp2==FALSE)
        {
          y.f.appr <- array(NA,dim=c(1,length(x.com.sim),r.t))
          y.g.appr <- array(NA,dim=c(1,length(x.com.sim),r.t))
          for(l in 1:r.t)
          {
            y.f.appr[1,,l] <- approx(x.f,function1[1,,l],xout=x.com.sim)$y
            y.g.appr[1,,l] <- approx(x.g,function2[1,,l],xout=x.com.sim)$y
          }
        }
        
    }
    
    if (var.temp2==FALSE)
    {
      
      #############################################
      # Compute the similarity index
      if (similarity.method=='d1.pearson' 
          || similarity.method=='d0.pearson')
      {
        distance <- 0
          for (l in 1:r.t)
          {
            distance <- distance + abs(sum(y.f.appr[1,,l]*y.g.appr[1,,l]))/(sqrt(sum(y.f.appr[1,,l]*y.f.appr[1,,l]))*sqrt(sum(y.g.appr[1,,l]*y.g.appr[1,,l])))
          }
          OUT <- distance/r.t
      }
        
      if (similarity.method=='d0.L2' || similarity.method=='d1.L2'
          || similarity.method=='d0.L2.centered' || similarity.method=='d1.L2.centered')
      {
        # If we want to "depure" functions from the mean
        # (d0.L2.centered or d1.L2.centered)
        if (similarity.method=='d0.L2.centered' || 
              similarity.method=='d1.L2.centered')      
        {
          for (l in 1:r.t)
          { 
            y.f.appr[1,,l] <- y.f.appr[1,,l] - mean(y.f.appr[1,,l],na.rm=TRUE)
            y.g.appr[1,,l] <- y.g.appr[1,,l] - mean(y.g.appr[1,,l],na.rm=TRUE)
          }
        }
        
        # Divide by the domain length
        D <- abs( max(x.com.sim)-min(x.com.sim) ) 
        
        delta <- x.com.sim[2]-x.com.sim[1] 
        # because here we always have an uniform grid
        
        L2distance <- 0
        for(l in 1:r.t)
        {
          # Approximate the L2norm with the rectangles form the "right"
          # We don't consider the first value (the grid is fine enough)
          L2distance <- L2distance + 
            #sqrt( sum(((y.f.appr[1,,l]-y.g.appr[1,,l])^2)*delta,na.rm=TRUE) )/D
            sqrt(sum(((y.f.appr[1,2:length(y.f.appr[1,,l]),l]-
                         y.g.appr[1,2:length(y.f.appr[1,,l]),l])^2)*delta,na.rm=TRUE) )/D
          
          #L2distance <- L2distance + 
          #sqrt( sum(((y.f.appr[1,,l]-y.g.appr[1,,l])^2)*delta,na.rm=TRUE) )/D[l]
        }
        OUT <- L2distance/r.t
        
      
    }
    
    } 
      
  }# fine if (var.temp==FALSE)
  
  OUT
}
