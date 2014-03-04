kma <-
function(
  x, y0 = NULL, y1 = NULL, 
  n.clust = 1, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means',
  seeds = NULL,
  optim.method = 'L-BFGS-B',
  span = 0.15, t.max = 0.1, m.max = 0.1, n.out = NULL,
  tol = 0.01, fence = TRUE ,iter.max = 100, show.iter = 0, nstart=1, return.all=FALSE) 
{
  n.clust.original <- n.clust
  y0.original <- y0  
  y1.original <- y1  
  similarity.method.original <- similarity.method
  
  # Control if similarity.method has a feasible value 
  similarity.method.available <- c(
    'd1.pearson', 'd0.L2', 'd0.pearson', 'd1.L2', 
    'd0.L2.centered', 'd1.L2.centered')
  if (!similarity.method %in% similarity.method.available) 
  {
    stop('Value of "similarity" not valid. Possibles choices are: ',
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
  {work.with.deriv.kma = 0}
  if (similarity.method=='d1.pearson' 
      || similarity.method=='d1.L2'
      || similarity.method=='d1.L2.centered')
  {work.with.deriv.kma = 1}
  
  ###########################################################
  ### Coherence controls: similarity.method VS warping.method
  ###########################################################
  
  ###########################
  # CONTROLS ON d1.pearson: y1 obblig, y0 facolt(se dato lo disegno con i template integrati, sennò nada)
  if( (similarity.method=='d1.pearson' 
       || similarity.method=='d1.pearson.mean') 
      && length(dim(y1))==0 )
  {
    stop('Evaluations of original function first derivatives (y1) must be provided to compute the chosen measure')
  }
  
  ###########################
  # CONTROLS ON d0.pearson: y0 obblig, y1 NO "has been ignored"
  if((similarity.method=='d0.pearson' 
      || similarity.method=='d0.pearson.mean') 
     && length(dim(y0))==0)
  {
    stop('Evaluations of original functions (y0) must be provided to compute the chosen measure')
  }
  if((similarity.method=='d0.pearson' 
      || similarity.method=='d0.pearson.mean') 
     && length(dim(y1))!=0)
  {
    warning('You provided evaluations of original function first derivatives (y1) but you chose ',similarity.method,
            ' as similarity.method. So value of y1 has been ignored')
  }  
  
  ###########################
  # CONTROLS ON d0.L2: y0 obblig, y1 NO "has been ignored"
  if((similarity.method=='d0.L2' 
      || similarity.method=='d0.L2.centered') 
     && length(dim(y0))==0)
  {
    stop('Evaluations of original functions (y0) must be provided to compute the chosen measure')
  }  
  if((similarity.method=='d0.L2' 
      || similarity.method=='d0.L2.centered') 
     && length(dim(y1))!=0)
  {
    warning('You provided evaluations of original function first derivatives (y1) but you chose ',similarity.method,
' as similarity.method. So the value of y1 has been ignored')
  }  
  if (warping.method!='NOalignment')
  {
    if ((similarity.method=='d0.L2' || 
           similarity.method=='d0.L2.centered') 
        && warping.method!='shift')
    {
      stop('L2 norm is a coherent similarity measure only for shift transformation of the abscissas. 
           So if you want to have similarity.method="L2" you have to choose warping.method="shift".')
    }
  }
     
  ###########################
  # CONTROLS ON d1.L2: y1 obblig, y0 facolt(se dato lo disegno con i template integrati, sennò nada)
  if((similarity.method=='d1.L2' || similarity.method=='d1.L2.centered') && length(dim(y1))==0 )
  {
    stop('Evaluations of original function first derivatives (y1) must be provided to compute the chosen measure')
  }
  if (warping.method!='NOalignment')
  {
    if ((similarity.method=='d1.L2' || similarity.method=='d1.L2.centered') && warping.method!='shift')
    {
      stop('d1.L2 norm is a coherent similarity measure only for shift transformation of the abscissas. 
           So if you want to have similarity.method="d1.L2" you have to choose warping.method="shift".')
    }
  }
  ############################# END COHERENCE CONTROLS
  
  ### registration function ###
  ###-----------------------###
  
  best_warping <- function(coeff){
    
    st  <- coeff[1]*x.reg[i,] + coeff[2]  # registered abscissa
    index.temp <- 0
    
    template.t <- templates[[k]][1,]
    b <- !is.na(template.t)
      
    if (work.with.deriv.kma==1)
    {data.t <- approx(st, data1[i,,1], xout=x.out)$y}
      
    if (work.with.deriv.kma==0)
    {data.t <- approx(st, data0[i,,1], xout=x.out)$y}
      
    a <- !is.na(data.t)
    sel <- a&b
    data.t <- data.t[sel]
   template.t <- templates[[k]][,sel]
      
    if (r>1) 
        {
        	data_def <- data.t
        	for(l in 2:r){
                
                if (work.with.deriv.kma==1)
       			{data.t <- approx(st, data1[i,,l], xout=x.out)$y}
          
        		if (work.with.deriv.kma==0)
        		{data.t <- approx(st, data0[i,,l], xout=x.out)$y}
         
         		data.t <- data.t[sel]
         		data_def <- rbind(data_def, data.t)		
    	   	}  
    	   	data.t <- data_def        
        }

    
    
    x.out.temp <- x.out[sel]
        
    if (work.with.deriv.kma==1)
    {        
      distance <- kma.similarity( x.f = x.out.temp, y1.f = t(data.t), x.g = x.out.temp, y1.g = t(template.t),
          similarity.method = similarity.method, unif.grid = FALSE)
    }
    if (work.with.deriv.kma==0)
    {
      distance <- kma.similarity(  x.f = x.out.temp, y0.f = t(data.t), x.g = x.out.temp, y0.g = t(template.t),
          similarity.method = similarity.method, unif.grid = FALSE)
    }
      
    index.temp <- index.temp + distance
   
    
        if(similarity.method == 'd1.pearson'  || similarity.method == 'd0.pearson'|| similarity.method == 'd1.pearson.mean' || similarity.method == 'd0.pearson.mean'
    )
    {
    	o <- -index.temp
    }else{
    	o <- index.temp
    	}

    o
  }
  
  best_warping_only.dilation <- function(coeff){
    st  <- coeff*x.reg[i,]  # registered abscissa
     index.temp <- 0
    
    template.t <- templates[[k]][1,]
    b <- !is.na(template.t)
      
    if (work.with.deriv.kma==1)
    {data.t <- approx(st, data1[i,,1], xout=x.out)$y}
      
    if (work.with.deriv.kma==0)
    {data.t <- approx(st, data0[i,,1], xout=x.out)$y}
      
    a <- !is.na(data.t)
    sel <- a&b
    data.t <- data.t[sel]
    template.t <- templates[[k]][,sel]
      
    if (r>1) 
        {
        	data_def <- data.t
        	for(l in 2:r){
                
                if (work.with.deriv.kma==1)
       			{data.t <- approx(st, data1[i,,l], xout=x.out)$y}
          
        		if (work.with.deriv.kma==0)
        		{data.t <- approx(st, data0[i,,l], xout=x.out)$y}
         
         		data.t <- data.t[sel]
         		data_def <- rbind(data_def, data.t)		
    	   	}    
    	   	data.t <- data_def      
        }
    
    
    x.out.temp <- x.out[sel]
        
    if (work.with.deriv.kma==1)
    {        
      distance <- kma.similarity( x.f = x.out.temp, y1.f = t(data.t), x.g = x.out.temp, y1.g = t(template.t),
          similarity.method = similarity.method, unif.grid = FALSE)
    }
    if (work.with.deriv.kma==0)
    {
      distance <- kma.similarity(  x.f = x.out.temp, y0.f = t(data.t), x.g = x.out.temp, y0.g = t(template.t),
          similarity.method = similarity.method, unif.grid = FALSE)
    }
      
    index.temp <- index.temp + distance
   
    
    if(similarity.method == 'd1.pearson'  || similarity.method == 'd0.pearson'|| similarity.method == 'd1.pearson.mean' || similarity.method == 'd0.pearson.mean'
    )
    {
    	o <- -index.temp
    }else{
    	o <- index.temp
    	}

    o
  }
  
  best_warping_only.shift <- function(coeff){
    st  <- x.reg[i,] + coeff  # registered abscissa
      index.temp <- 0
    
    template.t <- templates[[k]][1,]
    b <- !is.na(template.t)
      
    if (work.with.deriv.kma==1)
    {data.t <- approx(st, data1[i,,1], xout=x.out)$y}
      
    if (work.with.deriv.kma==0)
    {data.t <- approx(st, data0[i,,1], xout=x.out)$y}
      
    a <- !is.na(data.t)
    sel <- a&b
    data.t <- data.t[sel]
    template.t <- templates[[k]][,sel]
    
    if (r>1) 
        {
        	data_def <- data.t
        	for(l in 2:r){
                
                if (work.with.deriv.kma==1)
       			{data.t <- approx(st, data1[i,,l], xout=x.out)$y}
          
        		if (work.with.deriv.kma==0)
        		{data.t <- approx(st, data0[i,,l], xout=x.out)$y}
         
         		data.t <- data.t[sel]
         		data_def <- rbind(data_def, data.t)		
    	   	}    
    	   	data.t <- data_def      
        }
    
    
    x.out.temp <- x.out[sel] 
    if (work.with.deriv.kma==1)
    {        
      distance <- kma.similarity( x.f = x.out.temp, y1.f = t(data.t), x.g = x.out.temp, y1.g = t(template.t),
          similarity.method = similarity.method, unif.grid = FALSE)
    }
    if (work.with.deriv.kma==0)
    {
      distance <- kma.similarity(  x.f = x.out.temp, y0.f = t(data.t), x.g = x.out.temp, y0.g = t(template.t),
          similarity.method = similarity.method, unif.grid = FALSE)
    }
    index.temp <- index.temp + distance
   
    
        if(similarity.method == 'd1.pearson'  || similarity.method == 'd0.pearson'|| similarity.method == 'd1.pearson.mean' || similarity.method == 'd0.pearson.mean'
    )
    {
    	o <- -index.temp
    }else{
    	o <- index.temp
    	}
    o
 }
  
  
  ### Parameters handling ###
  ###---------------------###
  
  #if (fence==TRUE)
  #{
    cheb.locale <- 1.5
    cheb.globale <- 1.5
  #}else{
    #cheb.locale <- 1000
    #cheb.globale <- 1000
  #}
  
  lim.while =1
  loop.max = 3
  o.shift = 0
  o.dilation = 0
  reg <- FALSE;
  
  ## clusters number
  if(n.clust != floor(n.clust) | n.clust != ceiling(n.clust)){
    
    warning('"n.clust" must be integer, the value has been approximated to the nearest integer to continue')
    n.clust <- round(n.clust)
  }
  
  if(n.clust < 0){
    stop('"n.clust" must be positive')
  }
  
  # Controls if seeds has feasible values
  if (length(dim(y1))!=0){n.obs <- dim(y1)[1]}
  if (length(dim(y0))!=0){n.obs <- dim(y0)[1]}
  
  # Controls if seeds has a feasible dimension
  
  if (is.null(seeds))
  {
  	 seeds <- sample(1:n.obs,n.clust)
  }
  if (length(dim(seeds))==0)
  {
  	seeds <- as.matrix(t(seeds))
  }
  if (length(seeds[1,])>n.clust)
  {
    stop('number of columns of "seeds" must be inferior or equal to n.clust')
  }
  
 
  if (length(which(seeds>0 & seeds<=n.obs ))!= length(seeds))
  {
     stop('At least a value of "seeds" is not valid (is negative or null or superior to the number of observations)')
  }
  
  
  # Controls if other variables have feasible values
  if (center.method!='k-means' && center.method!='k-medoids')
  {
    stop('"center.method" has not a feasible value ("k-means" or "k-medoids"')
  }
  if (warping.method!='NOalignment' && warping.method!='affine' && warping.method!='shift' && warping.method!='dilation')
  {
    stop('"warping.method" has not a feasible value')
  }
  if (show.iter!=0 && show.iter!=1)
  {
    stop('"show.iter" must be 0 or 1')
  }
  if (iter.max<=0)
  {
    stop('Value of "iter.max" must be positive')
  }
  if(iter.max != floor(iter.max) | iter.max != ceiling(iter.max))
  {  
    warning('"iter.max" must be integer, the value has been approximated to the nearest integer to continue')
    iter.max <- round(iter.max)
  }
  if (t.max<=0 || t.max>=1) 
  {
    stop('"t.max" must be such that 0<t.max<1')    
  }
  if (m.max<=0 || m.max>=1) 
  {
    stop('"m.max" must be such that 0<m.max<1')    
  }
  

  method.available <- c('L-BFGS-B', 'SANN')
  if (length(optim.method)!=0) 
  {      
    if (!optim.method %in% method.available) 
    {
      stop('Value of "optim.method" not valid. It must be one of the following methods (optim package methods): ',
           '"', method.available[1], '"', ' ',
           '"', method.available[2], '"')
    }
  }
  
  # Control if x has a feasible value
  if (length(x)==0) {stop('You did not provide the abscissa x')}
  
  if (warping.method=='affine')
  {
    reg <- TRUE;
    o.shift <- 1;
    o.dilation <- 1;
  }
  if (warping.method=='shift')
  {
    reg <- TRUE;
    o.shift <- 1;
    o.dilation <- 0;
  }
  if (warping.method=='dilation')
  {
    reg <- TRUE;
    o.shift <- 0;
    o.dilation <- 1;
  }
  
  ## data and abscissa 
  if (length(dim(y1))==2) { y1 <- as.matrix(y1) } # Commento svedesi
  if (length(dim(y0))==2) { y0 <- as.matrix(y0) }
  
  #y1 <- as.matrix(y1)
  x <- as.matrix(x)
  if (dim(x)[1]==1 && dim(x)[2]>1) { x <- x}
  if (dim(x)[1]>1 && dim(x)[2]==1) { x <- t(x)}
  
  if (length(dim(y0))!=0)
  {
    if( sum(!is.na(match(is(y0),c('matrix','array'))))==0 ){
      stop('unknown data type, y0 must be an array or a matrix') 
    }
  }
  
  if (length(dim(y1))!=0)
  {
    if( sum(!is.na(match(is(y1),c('matrix','array'))))==0 ){
      stop('unknown data type, y1 must be an array or a matrix') 
    }
  }
  
  ###########################################################
  #### If x is a vector, apply this abscissa to all curves      
  x.temp <- x
  
  if(dim(x)[1]==1) {
    x.temp <- x
    #if(similarity.method=='d1.pearson' || similarity.method=='d1.L2')
    #{
      for(i in 1:(n.obs-1)) {x <- rbind(x,x.temp)}
    #}
    #if(similarity.method=='d0.pearson' || similarity.method=='d0.L2')
    #{for(i in 1:(dim(y0)[1]-1)) {x <- rbind(x,x.temp)}}
  }
  ###########################################################
  
  if (work.with.deriv.kma==1)
  {
    if(dim(x)[1]!=dim(y1)[1] || dim(x)[2]!=dim(y1)[2]){
      stop('Abscissa and function first derivatives dimensions must agree')
    }
    if (length(dim(y0))!=0)
    {
      if(dim(x)[1]!=dim(y0)[1] || dim(x)[2]!=dim(y0)[2]){
        stop('Abscissa and original function dimensions must agree')
      }
    }
  }
    
  if (work.with.deriv.kma==0)
  {
    if(dim(x)[1]!=dim(y0)[1] || dim(x)[2]!=dim(y0)[2]){
      stop('Abscissa and original functions dimensions must agree')
    }
  }

  
  #n.obs <- dim(x)[1]
  n.camp <- dim(x)[2]
  if (length(n.out)==0) { n.out = round(1.1*n.camp) }
  
  
  
  # Here we create data0 and/or data1, the object
  # we are going to work on
  
  if( length(dim(y1))!=0 )
  {
    if( length(dim(y1))==3 ){
      r <- r.deriv <- dim(y1)[3]
      data1 <- y1
    }else{
      r <- r.deriv <- 1    
      data1 <- array(0,dim=c(n.obs,n.camp,r))
      data1[,,1:r] <- y1
      # if r = 1, data is NOT a matrix but a three dimensional array
      data1 <- data1[,,1:r,drop=FALSE]
    }
  }
  
  if( length(dim(y0))!=0 )
  {
    if( length(dim(y0))==3 ){
      r <- r.orig <- dim(y0)[3]
      data0 <- y0
    }else{
      r <- r.orig <- 1    
      data0 <- array(0,dim=c(n.obs,n.camp,r))
      data0[,,1:r] <- y0
      # if r = 1, data is NOT a matrix but a three dimensional array
      data0 <- data0[,,1:r,drop=FALSE]
    }
  }  
  
  if( length(dim(y0))!=0 && length(dim(y1))!=0 )
  {
    if (r.deriv!=r.orig)
    {stop('original functions dimensions must agree with first derivatives')}
  }
  

  similarity.method.available_mean <- c( 'd0.L2.centered', 
                                        'd1.L2.centered')
   
  ### coherence checks (2) ###
  ###----------------------###
  
  ####################
  # Controls about y0
  ####################
    
  # Already checked
  #if (length(dim(y0))==2) { y0 <- as.matrix(y0) }
  
  if(length(dim(y0))!=0){
    if(dim(y0)[1]!=n.obs || dim(y0)[2]!=n.camp){
      stop('original functions dimensions must agree with abscissa')
    }
  }
  if(length(dim(y1))!=0){
    if(dim(y1)[1]!=n.obs || dim(y1)[2]!=n.camp){
      stop('function first derivatives dimensions must agree with abscissa')
    }
  }
  

  ## additional check: NAs in the data array
  if (length(dim(y1))!=0)
  {
    NA.data.y1 <- 0
    for(l in 1:r){NA.data.y1 <- NA.data.y1 + sum(rowSums(apply(data1[,,l], 2, is.na))==n.camp)}
    if(NA.data.y1 > 0){stop('at least one curve in the dataset is completely NA')}
  }
  if (length(dim(y0))!=0)
  {
    NA.data.y0 <- 0
    for(l in 1:r){NA.data.y0 <- NA.data.y0 + sum(rowSums(apply(data0[,,l], 2, is.na))==n.camp)}
    if(NA.data.y0 > 0){stop('at least one curve in the dataset is completely NA')}
  }
    
  
  ### method initialization ###
  ###-----------------------###
  
  OUT <- NULL
  
	for(s in 1:nstart)
	{

  
  
  # registered abscissa
  x.reg <- x
  # abscissa for the template curves
  x.out  <- seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE), length=n.out)
  
  # similarity index
  index <- rep(0,n.obs)
  index.old <- rep(-1,n.obs)
  index.list <- NULL    # list of similarity indexes along iterations
  
  # warping functions
  m.list <- NULL   # list of shifts along iterations
  t.list <- NULL   # list of dilations along iterations
  m.list.temp <- NULL
  t.list.temp <- NULL
  limitewhile <- lim.while
  count.loop <- loop.max
  only.shift <- o.shift
  only.dilation <- o.dilation
  
  # warping functions range of variation
  # (shifts)
  min.temp <- diff(range(x,na.rm=TRUE))
  for(i in 1:n.obs){ min.temp <- min(min.temp, diff(range(x[i,],na.rm=TRUE)) ) }
  
  # the first component is dilation, the second is shift
  lower.warp <- c(1-m.max,-t.max*min.temp)
  upper.warp <- c(1+m.max, t.max*min.temp)
  flagglobale <- 0 # Control Variable: not important if you don't unterstand
  coeff.cheb.locale <- cheb.locale # local "whisker"
  coeff.cheb.globale <- cheb.globale # global "whisker"
  
  # vector containing group assignments (labels)
  labels <- c(rep(1,n.obs-1),2)
  labels.old <- rep(0,n.obs)
  
  # iterations index
  iter <- 0
  
  # templates (chosen among curves)
  templates <- NULL
  
  # Here centers are chosen for the first time
  # (randomly if seeds=NULL, otherwise as specified 
  # by user in "seeds")
  if (dim(seeds)[1]<s)
  {
   		 # Random choice of seeds(=centers indexes) of first iteration
   		 u <- sample(1:n.obs,n.clust)
  } else{
  		 u <- seeds[s,]	
  }
  
  #   Here we approximate centers in the abscissas x.out 
  #   (hence on n.out points between min(x) and max(x))
  for(j in 1:n.clust){
    temp <- NULL
    for(i in 1:r){
      
      if (work.with.deriv.kma==1)
      {temp1 <- approx(x[u[j],],data1[u[j],,i], xout=x.out)$y}
      
      if (work.with.deriv.kma==0)
      {temp1 <- approx(x[u[j],],data0[u[j],,i], xout=x.out)$y}
      
      temp <- rbind(temp,temp1)
    }
    templates <- c(templates, list(temp))
  }
  
  #################################################################
  if(center.method=='k-means')
  {
    similarity.orig.viacenterscluster <- 0
    similarity.orig.vialoess <- 1
    similarity.orig.viamedoid <- 0
  }
  if(center.method=='k-medoids')
  {
    similarity.orig.viacenterscluster <- 0
    similarity.orig.vialoess <- 0
    similarity.orig.viamedoid <- 1
    
  }
  
  #############################################
  #############################################
  #############################################
  #############################################
  
  #############################################
  # Compute the similarity between each curves and 
  # the medoid of original data
  
  x.center.orig <- seq(max(x, na.rm=TRUE),min(x, na.rm=TRUE),
                       length=n.out)
  #x.center.orig <- x.out
  
  y0.center.orig <- NULL
  y1.center.orig <- NULL
  
  if (similarity.orig.viamedoid==1)
  {
    # Where all functions are defined 
    # (the common domaine)
    x.com  <- seq(max(x[,1], na.rm=TRUE),min(x[,dim(x)[2]], na.rm=TRUE),
                  length=n.out)
    
    
    # Compute registered functions
    data0.reg<-array(0,dim=c(n.obs,n.out,r))
    data1.reg<-array(0,dim=c(n.obs,n.out,r))
    for(l in 1:r){
      for(k in 1:n.obs){
        if (work.with.deriv.kma==1)
        {data1.reg[k,,l]<-approx(x[k,],data1[k,,l],xout=x.com)$y}
        if (work.with.deriv.kma==0)
        {data0.reg[k,,l]<-approx(x[k,],data0[k,,l],xout=x.com)$y}
        
      }
    }
    
    #x.center.orig <- x.out
    #y0.center.orig <- NULL
    #y1.center.orig <- NULL
    #y.center.orig <- NULL
    
    similarity.orig <- NULL
    #labels.iter1 <- rep(0,n.obs) 
    data.selec <- 1:n.obs
      
    y0 <- NULL
      
    distanze.matrix <- matrix(0,length(data.selec),length(data.selec))
    distanze<-rep(0,length(data.selec))
    
    
    # First we find the medoid of original data, saving
    # all the similarity indexex in a matrix called
    # distanze.matrix
      for(i in 1:length(data.selec)){
        for(j in 1:length(data.selec)){
          
          b <- !is.na(data1.reg[i,,1])
          a <- !is.na(data1.reg[j,,1])
       	  sel <- a&b          
          if (work.with.deriv.kma==1)
          #{ecco <- kma.similarity(data1.reg[i,,l],data1.reg[j,,l],similarity.method)}
          { 	
          	b <- !is.na(data1.reg[i,,1])
          	a <- !is.na(data1.reg[j,,1])
       	  	sel <- a&b          
            ecco <- kma.similarity(
            x.f = x.com[sel], y1.f = data1.reg[i,sel,],
            x.g = x.com[sel], y1.g = data1.reg[j,sel,],
            similarity.method = similarity.method, unif.grid = FALSE)
          }
          
          if (work.with.deriv.kma==0)
          #{ecco <- kma.similarity(data0.reg[i,,l],data0.reg[j,,l],similarity.method)}
          {
          	b <- !is.na(data0.reg[i,,1])
          	a <- !is.na(data0.reg[j,,1])
       	  	sel <- a&b 
            ecco <- kma.similarity(
              x.f = x.com[sel], y0.f = data0.reg[i,sel,],
              x.g = x.com[sel], y0.g = data0.reg[j,sel,],
              similarity.method = similarity.method, unif.grid = FALSE)
          }
          distanze[i]<-distanze[i]+ecco
          distanze.matrix[i,j] <- distanze.matrix[i,j] + ecco
        }
      }  
    
    
    if(similarity.method == 'd1.pearson' 
       || similarity.method == 'd0.pearson'
         ) 
    { m<-which.max(distanze) }
    
    if(similarity.method == 'd1.L2' 
       || similarity.method == 'd0.L2'
       || similarity.method == 'd1.L2.centered'
       || similarity.method == 'd0.L2.centered'
    ) 
    { m<-which.min(distanze) }
    
    remind.center.orig <- data.selec[m]
    

    for(l in 1:r)
    {
      # Here we approx the medoid of original similarity 
      # on x.center.orig
      if (work.with.deriv.kma==1)
      {temp<-approx(x[data.selec[m],],data1[data.selec[m],,l],xout=x.center.orig)$y}
      if (work.with.deriv.kma==0)
      {temp<-approx(x[data.selec[m],],data0[data.selec[m],,l],xout=x.center.orig)$y}
      
      y0 <- rbind(y0, temp)
      
    }
    
    
    # Here we end the research and approximation of the medoid of original data
    
    # Now we compute (in reality we already did that, 
    # now we copy the values in similarity.orig) the 
    # similarity between each function
    # and the medoid we just found
    if (work.with.deriv.kma==1)
    {
      y1.center.orig <- c(y1.center.orig, list(y0))
    }
    if (work.with.deriv.kma==0)
    {
      y0.center.orig <- c(y0.center.orig, list(y0))
    }
    #y.center.orig <- c(y.center.orig, list(y0))
    
    #y.center.orig <- y0
    similarity.orig <- distanze.matrix[m,]
    
  }
  
  #############################################
  # Compute the similarity between each curves and 
  # the Frechet mean obtained via loess procedure.
  if (similarity.orig.vialoess==1)
  {
    # Here we compute the loess for orginal data and compute 
    # similarities indexes
    templates.iter1.vialoess <- NULL
    
    data.selec <- 1:n.obs
    y0 <- NULL
        
    for(l in 1:r){
      
      if (work.with.deriv.kma==1)
      {temp.iter1.vialoess <- loess(as.numeric(data1[data.selec,,l]) ~ as.numeric(x[data.selec,]), span=span)}
      
      if (work.with.deriv.kma==0)
      {temp.iter1.vialoess <- loess(as.numeric(data0[data.selec,,l]) ~ as.numeric(x[data.selec,]), span=span)}
      
      temp.iter1.vialoess <- predict(temp.iter1.vialoess, x.center.orig)
      y0 <- rbind(y0, temp.iter1.vialoess)
      
    }
    
    templates.iter1.vialoess <- c(templates.iter1.vialoess, list(y0))
    
    # Calculus of quantities from first iteration
    # x.center.orig <- x.out
    if (work.with.deriv.kma==1)
    {
      y1.center.orig <-  templates.iter1.vialoess
    }
    if (work.with.deriv.kma==0)
    {
      y0.center.orig <-  templates.iter1.vialoess
    }
    #y.center.orig <-  templates.iter1.vialoess
    
    similarity.orig <- NULL  
    
    data.selec <- 1:n.obs
    
    for (i in 1:n.obs)
    {
      #index.temp.iter1.vec <- rep(0,n.clust)
      index.temp.iter1.vec <- rep(0,1)
      #for(k in 1:n.clust){
 
        index.temp.iter1 <- 0
        
        #template.t <- templates[[k]][l,]
        template.t <- templates.iter1.vialoess[[1]][1,]
        b <- !is.na(template.t)
        if (work.with.deriv.kma==1)
        {data.t <- approx(x[i,], data1[i,,1], xout=x.center.orig)$y}
          
        if (work.with.deriv.kma==0)
        {data.t <- approx(x[i,], data0[i,,1], xout=x.center.orig)$y}
          
        a <- !is.na(data.t)
        sel <- a&b
        data.t <- data.t[sel]
        
        x.center.orig.temp <- x.center.orig[sel]
        
        template.t <- templates.iter1.vialoess[[1]][,sel]
         
        if (r>1) 
        {
        	data_def <- data.t
        	for(l in 2:r){
                
                if (work.with.deriv.kma==1)
       			{data.t <- approx(x[i,], data1[i,,l], xout=x.center.orig)$y}
          
        		if (work.with.deriv.kma==0)
        		{data.t <- approx(x[i,], data0[i,,l], xout=x.center.orig)$y}
         
         		data.t <- data.t[sel]
         		data_def <- rbind(data_def, data.t)
         	}      
         	data.t <- data_def    
        }
      
         x.center.orig.temp <- x.center.orig[sel]
         
        if (work.with.deriv.kma==1)
          {
            distance <- kma.similarity(x.f = x.center.orig.temp, y1.f = t(data.t),  x.g = x.center.orig.temp, y1.g = t(template.t), similarity.method = similarity.method, unif.grid = FALSE)
          }
          
          if (work.with.deriv.kma==0)
          {
            distance <- kma.similarity(  x.f = x.center.orig.temp, y0.f = t(data.t), x.g = x.center.orig.temp, y0.g = t(template.t),  similarity.method = similarity.method, unif.grid = FALSE)
          }
          
          index.temp.iter1 <- index.temp.iter1 + distance
          similarity.orig  <- c(similarity.orig,index.temp.iter1)
        }
     
    }
    
  # END Calculus of quantities from first iteration
 
  
  #################################################################
  
  m.final.temp <- rep(1, n.obs)
  t.final.temp <- rep(0, n.obs)
  
  ### k-mean alignment method ###
  ###-------------------------###
  
  still.in <- TRUE
    
   if (similarity.method == 'd1.L2'  || similarity.method == 'd0.L2' || similarity.method == 'd1.L2.centered' || similarity.method == 'd0.L2.centered')
  {
  	index.old <- rep(1, n.obs)
  }  
    
  while(( (((sum((index - index.old) < tol) < n.obs) && (similarity.method == 'd1.pearson' || similarity.method == 'd0.pearson')) || ((sum((-index + index.old) < tol) < n.obs) && (similarity.method == 'd1.L2'  || similarity.method == 'd0.L2' || similarity.method == 'd1.L2.centered' || similarity.method == 'd0.L2.centered')) ) && still.in) && iter < iter.max)
  {   
    
    # labels
    labels.old <- labels
    # similarity measure
    index.old <- index
    
     if ((similarity.method == 'd1.L2'  || similarity.method == 'd0.L2' || similarity.method == 'd1.L2.centered' || similarity.method == 'd0.L2.centered' || similarity.method=='d0ed1.L2') && iter==0)
 	{
  		index.old <- rep(10^5, n.obs)
  	}
  
    templates.old <- templates
    
    
    index <- NULL
    
    # iterations index
    iter <- iter + 1
  
    ### assignment and alignment step ###
    ###-------------------------------###
    teigs <- NULL # shift
    meigs <- NULL # dilations
    
    if (show.iter==1)
    {
      output.text0 <- paste("Starting step: ", s, collapse = "\t")
  	  output.text02 <- paste("Seeds: ", u)
      output.text1 <- paste('Num.cluster: ', n.clust, collapse='\t')
      output.text2 <- paste('Alignment: ', warping.method, collapse='\t')
      output.text3 <- paste('Iteration: ', iter, collapse='\t')
      #print(   '*********************************************');
      print(output.text0)
   	  print(output.text02)
      print( output.text1 );
      print( output.text2 );
      print( output.text3 );
      print('*********************************************');
    
    }
    
    for(i in 1:n.obs){
      
      warping_temp <- matrix(-100,2,length(templates))
      index_temp <- NULL
      
      for(k in 1:length(templates)){        
        
        if(reg)
        {
          # Optimization step
          
          # shift
          if (only.shift==1 && only.dilation==0)
          {
            if (optim.method=='L-BFGS-B')
            {
              result <- optim(c(0), best_warping_only.shift, method=optim.method, lower = lower.warp[2], upper = upper.warp[2])
              if (result$convergence!=0)
              {
              	result <- optim(c(0), best_warping_only.shift, method=optim.method,control=list(parscale=1e-15), lower = lower.warp[2], upper = upper.warp[2])
              }
            }
            if (optim.method=='SANN')
            {
              f.sann <- function(x){
                #step1 <- runif(1, lower.warp[1], upper.warp[1])
                step2 <- runif(1, lower.warp[2], upper.warp[2])
                x <- step2
                return(x)
              }
              result <- optim(c(0), best_warping_only.shift, gr=f.sann, method=optim.method)
            }
            
            warping_temp[2,k] <- result$par
            warping_temp[1,k] <- 1
          }
          
          # dilation
          if (only.dilation==1 && only.shift==0) 
          {
            if (optim.method=='L-BFGS-B')
            {
              result <- optim(c(1), best_warping_only.dilation, method=optim.method, lower = lower.warp[1], upper = upper.warp[1])
                if (result$convergence!=0)
               {
               	result <- optim(c(1), best_warping_only.dilation, method=optim.method,control=list(parscale=1e-15), lower = lower.warp[2], upper = upper.warp[2])
              }

            } 
            if (optim.method=='SANN')
            {
              f.sann <- function(x){
                step1 <- runif(1, lower.warp[1], upper.warp[1])
                #step2 <- runif(1, lower.warp[2], upper.warp[2])
                x <- step1
                return(x)
              }
              result <- optim(c(1), best_warping_only.dilation, gr=f.sann, method=optim.method)
            }
            
            warping_temp[1,k] <- result$par
            warping_temp[2,k] <- 0
          }
          
          # affine
          if (only.dilation==only.shift) 
          {
            if (optim.method=='L-BFGS-B')
            { 
              result <- optim(c(1,0), best_warping, method=optim.method, lower = lower.warp, upper = upper.warp)
               if (result$conergence!=0)
               {
               	result <- optim(c(1,0), best_warping, method=optim.method,control=list(parscale=c(1e-15, 1e-15)), lower = lower.warp[2], upper = upper.warp[2])
              }

            }
            if (optim.method=='SANN')
            {
              f.sann <- function(x){
                step1 <- runif(1, lower.warp[1], upper.warp[1])
                step2 <- runif(1, lower.warp[2], upper.warp[2])
                x <- c(step1,step2)
                return(x)
              }
              result <- optim(c(1,0), best_warping, gr=f.sann, method=optim.method)
            }
            
            warping_temp[,k] <- result$par
          }
          
          
          if(similarity.method == 'd1.pearson' 
             || similarity.method == 'd0.pearson'
          ) 
          {index_temp <- c(index_temp, -result$value)}
          
          if(similarity.method == 'd1.L2' 
             || similarity.method == 'd0.L2'
             || similarity.method == 'd1.L2.centered'
             || similarity.method == 'd0.L2.centered'
          ) 
          {index_temp <- c(index_temp, +result$value)}

        }else{
          
          warping_temp[,k] <- c(1,0)
          temp <- 0
          template.t <- templates[[k]][1,]
          b <- !is.na(template.t)
          
          if (work.with.deriv.kma==1)
          {data.t <- approx(x.reg[i,], data1[i,,1], xout=x.out)$y}
           
          if (work.with.deriv.kma==0)
          {data.t <- approx(x.reg[i,], data0[i,,1], xout=x.out)$y}
            
          a <- !is.na(data.t)    
          sel <- a&b 
          data.t <- data.t[sel]       
          template.t <- templates[[k]][,sel]
          x.out.temp <- x.out[sel]
          
          if (r>1) 
        	{

			data_def <- data.t
          	for(l in 2:r){
            
            	if (work.with.deriv.kma==1)
          		{data.t <- approx(x.reg[i,], data1[i,,l], xout=x.out)$y}
          	 
          		if (work.with.deriv.kma==0)
          		{data.t <- approx(x.reg[i,], data0[i,,l], xout=x.out)$y}
            
           		data.t <- data.t[sel]
           		data_def <- rbind(data_def, data.t) 
           	}
           	data.t <- data_def
          	}
          	
          	x.out.temp <- x.out[sel]
          	
            #temp <- temp + kma.similarity(data.t,template.t,similarity.method)
            
            if (work.with.deriv.kma==1)
            {
              temp <- temp + kma.similarity( x.f = x.out.temp, y1.f = t(data.t),  x.g = x.out.temp, y1.g = t(template.t),
                similarity.method = similarity.method, unif.grid = FALSE)   
            }
            
            if (work.with.deriv.kma==0)
            {
              temp <- temp + kma.similarity(  x.f = x.out.temp, y0.f = t(data.t),  x.g = x.out.temp, y0.g = t(template.t),
                similarity.method = similarity.method, unif.grid = FALSE)   
            }
            
          index_temp <- c(index_temp, temp)
          
        }
        
      }# end cycle on the clusters
            
      if(similarity.method == 'd1.pearson' 
         || similarity.method == 'd0.pearson') 
      {
        # similarity index
        index  <- c(index,  max(index_temp))
        # cluster label
        labels[i] <- which.max(index_temp)
      }
      if(similarity.method == 'd1.L2' || similarity.method == 'd0.L2' || similarity.method == 'd1.L2.centered' || similarity.method == 'd0.L2.centered'
      ) 
      {
        # similarity index
        index  <- c(index,  min(index_temp))
        # cluster label
        labels[i] <- which.min(index_temp)
      }
      
      # optimal warping function shift and dilation
      meigs <- c(meigs, warping_temp[1,labels[i]])
      teigs <- c(teigs, warping_temp[2,labels[i]])
      
    }# end cycle on the curves
    
    if (((sum(index)-sum(index.old))>0 && (similarity.method == 'd1.pearson' || similarity.method == 'd0.pearson') ) || ((sum(index)-sum(index.old))<0 && (similarity.method == 'd1.L2'  || similarity.method == 'd0.L2' || similarity.method == 'd1.L2.centered' || similarity.method == 'd0.L2.centered')) || (center.method=='k-means' && iter==2)|| iter==1)
    {
    
    if (fence==TRUE)
    {
      #####################       FENCE       ######################################
      if (iter==1) {limitewhile <- 1} else {limitewhile <- lim.while}
      
      if(reg){
        
        contatorewhile <- 1
        
        # Normalize shifts and dilations
        for(k in unique(labels)){
          
          data.selec <- which(labels==k)
          meigs[data.selec] <- meigs[data.selec]/mean(meigs[data.selec])
          teigs[data.selec] <- (teigs[data.selec] - mean(teigs[data.selec]))/mean(meigs[data.selec])
          
        }
        
        B1.meigs <- B2.meigs <- B1.teigs <- B2.teigs <- NULL
        
        ########################   CALCOLO BOUNDS ITERAZIONE     ########################
        Q3.meigs <- summary(meigs)[5] # third quartile
        Q1.meigs <- summary(meigs)[2] # first quartile
        Q3.teigs <- summary(teigs)[5] # third quartile
        Q1.teigs <- summary(teigs)[2] # first quartile
        B1.meigs <- Q1.meigs - coeff.cheb.locale*(Q3.meigs -Q1.meigs)
        B2.meigs <- Q3.meigs + coeff.cheb.locale*(Q3.meigs -Q1.meigs)
        B1.teigs <- Q1.teigs - coeff.cheb.locale*(Q3.teigs -Q1.teigs)
        B2.teigs <- Q3.teigs + coeff.cheb.locale*(Q3.teigs -Q1.teigs)
        
        #       B1.meigs <- mean(meigs) - coeff.cheb.locale*sd(meigs)
        #       B2.meigs <- mean(meigs) + coeff.cheb.locale*sd(meigs)
        #       B1.teigs <- mean(teigs) - coeff.cheb.locale*sd(teigs)
        #       B2.teigs <- mean(teigs) + coeff.cheb.locale*sd(teigs)
        
        if (only.shift==1 && only.dilation==0) { ind.out.it <- which( (teigs > B2.teigs) | (teigs < B1.teigs) ) }
        if (only.dilation==1 && only.shift==0) { ind.out.it <- which( (meigs > B2.meigs) | (meigs < B1.meigs) ) }
        if ( only.shift == only.dilation )
        {
          ind.out.it <- which( (meigs > B2.meigs) | (meigs < B1.meigs) | (teigs > B2.teigs) | (teigs < B1.teigs) )
        }
        
        ########################    CALCOLO BOUNDS GLOBALI       ########################
        ind.out.temp <- NULL     
        contatoreloop <- 0
        ind.out.OLD <- n.obs + 1
        
        while (contatorewhile <= limitewhile) # limitewhile=3
        {
          m.list.temp <- m.list
          t.list.temp <- t.list
          m.list.temp <- c(m.list.temp,list(meigs))
          t.list.temp <- c(t.list.temp,list(teigs))
          
          if (iter > 1)
          {
            # calcolo shift e dilat globali
            m.final.temp <- rep(1,n.obs)
            t.final.temp <- rep(0,n.obs)
            
            for(h in 1:(iter)){
              m.final.temp <- m.final.temp*m.list.temp[[h]]
              t.final.temp <- t.final.temp*m.list.temp[[h]] + t.list.temp[[h]]
            }
            ## NORMALIZE FINAL.TEMP
            for(k in 1:n.clust){
              data.selec <- which(labels==k)
              m.final.temp[data.selec] <- m.final.temp[data.selec]/mean(m.final.temp[data.selec])
              t.final.temp[data.selec] <- (t.final.temp[data.selec] - mean(t.final.temp[data.selec]))/mean(m.final.temp[data.selec])
            }
            
            # Calcolo Q1 e Q3 globali
            Q3.m.final.temp <- summary(m.final.temp)[5] # third quartile
            Q1.m.final.temp <- summary(m.final.temp)[2] # first quartile
            Q3.t.final.temp <- summary(t.final.temp)[5] # third quartile
            Q1.t.final.temp <- summary(t.final.temp)[2] # first quartile
            
            B1.m.final.temp <- Q1.m.final.temp - coeff.cheb.globale*(Q3.m.final.temp - Q1.m.final.temp)
            B2.m.final.temp <- Q3.m.final.temp + coeff.cheb.globale*(Q3.m.final.temp - Q1.m.final.temp)
            B1.t.final.temp <- Q1.t.final.temp - coeff.cheb.globale*(Q3.t.final.temp - Q1.t.final.temp)
            B2.t.final.temp <- Q3.t.final.temp + coeff.cheb.globale*(Q3.t.final.temp - Q1.t.final.temp)
            
            #         B1.m.final.temp <- mean(m.final.temp) - coeff.cheb.globale*sd(m.final.temp)
            #         B2.m.final.temp <- mean(m.final.temp) + coeff.cheb.globale*sd(m.final.temp)
            #         B1.t.final.temp <- mean(t.final.temp) - coeff.cheb.globale*sd(t.final.temp)
            #         B2.t.final.temp <- mean(t.final.temp) + coeff.cheb.globale*sd(t.final.temp)
            
            if (only.shift == 1 && only.dilation==0) { ind.out.temp <- which( (t.final.temp > B2.t.final.temp) | (t.final.temp < B1.t.final.temp) ) }
            if (only.dilation == 1 && only.shift==0) { ind.out.temp <- which( (m.final.temp > B2.m.final.temp) | (m.final.temp < B1.m.final.temp) ) }
            if ( only.shift == only.dilation ) 
            {
              ind.out.temp <- which( (m.final.temp > B2.m.final.temp) | (m.final.temp < B1.m.final.temp) | (t.final.temp > B2.t.final.temp) | (t.final.temp < B1.t.final.temp) )
            }
          }
          else {
            B1.m.final.temp <- B1.meigs
            B2.m.final.temp <- B2.meigs
            B1.t.final.temp <- B1.teigs
            B2.t.final.temp <- B2.teigs
            ind.out.temp <- ind.out.it
          }
          # FINE CALCOLO BOUNDS GLOBALI
          
          if (contatorewhile > 1) { ind.out <- sort(ind.out.temp); }
          else { ind.out <- sort ( unique (c( ind.out.it, ind.out.temp) ) ) }
          if (length(ind.out)==0) break;
          if ( isTRUE(all.equal(ind.out, ind.out.OLD)) ) {contatoreloop <- contatoreloop + 1}
          if (contatoreloop == count.loop) break;
          ind.out.OLD <- ind.out
          contatorewhile <- contatorewhile +1
          
          lwn <- matrix(NA,length(ind.out),2)
          uwn <- matrix(NA,length(ind.out),2)
          
          if (iter==1) 
          {
            lwn[,1] <- rep(B1.meigs,length(ind.out))
            uwn[,1] <- rep(B2.meigs,length(ind.out))
            lwn[,2] <- rep(B1.teigs,length(ind.out))
            uwn[,2] <- rep(B2.teigs,length(ind.out))
          }
          if (iter > 1)
          {
            i.seq <- 0
            for (gh in ind.out)
            {
              i.seq = i.seq + 1
              
              lwn[i.seq,1] <- B1.m.final.temp / (m.final.temp[gh] / meigs[gh] )
              uwn[i.seq,1] <- B2.m.final.temp / (m.final.temp[gh] / meigs[gh] )
              lwn[i.seq,2] <- B1.t.final.temp - (t.final.temp[gh] - teigs[gh] )
              uwn[i.seq,2] <- B2.t.final.temp - (t.final.temp[gh] - teigs[gh] )
            }
          }
          
          i.seq <- 0
          
          for (i  in ind.out)
          {
            i.seq = i.seq + 1
            
            warping_temp <- matrix(-100,2,length(templates))
            index_temp <- NULL
            
            #  CONTROLLO INTERSEZIONI E FISSA I LOWER.WARP.NEW E UPPER.WARP.NEW
            flag1<-NULL
            if (only.shift==1 && only.dilation==0)
            {flag1 <- (uwn[i.seq,2] < B1.teigs) | (lwn[i.seq,2] > B2.teigs)        }
            if (only.shift==0 && only.dilation==1)
            {flag1 <- (uwn[i.seq,1] < B1.meigs) | (lwn[i.seq,1] > B2.meigs)        }
            if (only.shift==only.dilation)
            {flag1 <- (uwn[i.seq,1] < B1.meigs) | (lwn[i.seq,1] > B2.meigs) | (uwn[i.seq,2] < B1.teigs) | (lwn[i.seq,2] > B2.teigs)        }
            
            if (flag1) # se non vi e' intersezione
            { 
              flagglobale <- flagglobale+1
              lower.warp.new <- c( lwn[i.seq,1] - 10^(-4) , lwn[i.seq,2] - 10^(-4) )
              upper.warp.new <- c( uwn[i.seq,1] + 10^(-4) , uwn[i.seq,2] + 10^(-4) )
            }
            
            if (flag1==FALSE) # se vi e' intersezione
            {
              lower.warp.new <- c( ( max(B1.meigs, lwn[i.seq,1]) - 10^(-4) ), (max(B1.teigs, lwn[i.seq,2]) - 10^(-4)) )
              upper.warp.new <- c( ( min(B2.meigs, uwn[i.seq,1]) + 10^(-4) ), (min(B2.teigs, uwn[i.seq,2]) + 10^(-4)) )
            }
            
            #  FINE CONTROLLO INTERSEZIONI E FISSA I LOWER.WARP.NEW E UPPER.WARP.NEW
            for(k in 1:length(templates))
            {
              init1<-mean(upper.warp.new[1],lower.warp.new[1]); 
              init2<-mean(upper.warp.new[2],lower.warp.new[2]);
              
              if ( (lower.warp.new[1]<=1) && (upper.warp.new[1]>=1)) { init1<-1; }
              if ( (lower.warp.new[2]<=0) && (upper.warp.new[2]>=0)) { init2<-0; }
              
              # Optimization step
              
              # shift
              if (only.shift==1 && only.dilation==0)
              {
                if (optim.method=='L-BFGS-B')
                {
                  result <- optim(c(init2), best_warping_only.shift, method=optim.method, lower = lower.warp.new[2], upper = upper.warp.new[2])
                   if (result$convergence!=0)
                   {
                   	result <- optim(c(init2), best_warping_only.shift, method=optim.method,control=list(parscale=1e-15), lower = lower.warp[2], upper = upper.warp[2])
                   }
                }
                if (optim.method=='SANN')
                {
                  f.sann <- function(x){
                    #step1 <- runif(1, lower.warp.new[1], upper.warp.new[1])
                    step2 <- runif(1, lower.warp.new[2], upper.warp.new[2])
                    x <- c(step2)
                    return(x)
                  }
                  result <- optim(c(init2), best_warping_only.shift, gr=f.sann, method=optim.method)
                }
                warping_temp[2,k] <- result$par
                warping_temp[1,k] <- 1
              }
              
              # dilation
              if (only.dilation==1 && only.shift==0) 
              {
                if (optim.method=='L-BFGS-B')
                {
                  result <- optim(c(init1), best_warping_only.dilation, method=optim.method, lower = lower.warp.new[1], upper = upper.warp.new[1])
                  if (result$convergence!=0)
                   {
                   	result <- optim(c(init1), best_warping_only.dilation, method=optim.method,control=list(parscale=1e-15), lower = lower.warp[2], upper = upper.warp[2])
                 
                   }
                }
                if (optim.method=='SANN')
                {
                  f.sann <- function(x){
                    step1 <- runif(1, lower.warp.new[1], upper.warp.new[1])
                    #step2 <- runif(1, lower.warp.new[2], upper.warp.new[2])
                    x <- c(step1)
                    return(x)
                  }
                  result <- optim(c(init1), best_warping_only.dilation, gr=f.sann, method=optim.method)
                }
                warping_temp[1,k] <- result$par
                warping_temp[2,k] <- 0
              }
              
              # affine
              if (only.dilation==only.shift) 
              {
                if (optim.method=='L-BFGS-B')
                {
                  result <- optim(c(init1,init2), best_warping, method=optim.method, lower = lower.warp.new, upper = upper.warp.new)
                   if (result$conergence!=0)
                   {
                   	result <- optim(c(init1, init2), best_warping, method=optim.method,control=list(parscale=c(1e-15, 1e-15)), lower = lower.warp[2], upper = upper.warp[2])
                   }
                }
                if (optim.method=='SANN')
                {
                  f.sann <- function(x){
                    step1 <- runif(1, lower.warp.new[1], upper.warp.new[1])
                    step2 <- runif(1, lower.warp.new[2], upper.warp.new[2])
                    x <- c(step1, step2)
                    return(x)
                  }
                  result <- optim(c(init1,init2), best_warping, gr=f.sann, method=optim.method)
                }
                
                warping_temp[,k] <- result$par
              }
              
              if(similarity.method == 'd1.pearson' 
                 || similarity.method == 'd0.pearson'
              ) 
              {index_temp <- c(index_temp, -result$value)}
              
              if(similarity.method == 'd1.L2' 
                 || similarity.method == 'd0.L2'
                 || similarity.method == 'd1.L2.centered'
                 || similarity.method == 'd0.L2.centered'
              ) 
              {index_temp <- c(index_temp, +result$value)} 
              
            }
            
            if(similarity.method == 'd1.pearson' 
               || similarity.method == 'd0.pearson'
            ) 
            {
              # similarity index
              index[i]  <- max(index_temp)
              # cluster label
              labels[i] <- which.max(index_temp)
            }
            
            if(similarity.method == 'd1.L2' 
               || similarity.method == 'd0.L2'
               || similarity.method == 'd1.L2.centered'
               || similarity.method == 'd0.L2.centered'
            ) 
            {
              # similarity index
              index[i]  <- min(index_temp)
              # cluster label
              labels[i] <- which.min(index_temp)
            }
            
            # optimal warping function shift and dilation
            meigs[i] <- warping_temp[1,labels[i]]
            teigs[i] <- warping_temp[2,labels[i]]
          }
          
          # Normalize shifts and dilations
          for(k in unique(labels)){
            
            data.selec <- which(labels==k)
            meigs[data.selec] <- meigs[data.selec]/mean(meigs[data.selec])
            teigs[data.selec] <- (teigs[data.selec] - mean(teigs[data.selec]))/mean(meigs[data.selec])
            
          }
          
        } # fine for (contatorefor in 1:3) FINE DEL BREAK
        
      }# fine if(reg)
      
      ########    END FENCE     ###########################
    }
    
    
    ### apply the selected registration to the whole data set ###
    ###-------------------------------------------------------###
    
    for(i in 1:n.obs){ x.reg[i,] <- meigs[i]*x.reg[i,] + teigs[i] }
    
    ### list of shifts and dilations along iterations
    m.list <- c(m.list,list(meigs))
    t.list <- c(t.list,list(teigs))
    m.list.temp <- NULL
    t.list.temp <- NULL
    m.list.temp <- m.list
    t.list.temp <- t.list
    
    ### template identification step ###
    ###------------------------------###
        
    # templates common abscissa
    x.out  <- seq(min(x.reg, na.rm=TRUE),max(x.reg, na.rm=TRUE), length=n.out)
    
    # qui templates
    
    if (center.method=='k-means')
    {
      templates <- NULL
      
      for(k in sort(unique(labels))){
        
        data.selec <- which(labels==k)
        y0 <- NULL
        
        for(l in 1:r){
          
          if (work.with.deriv.kma==1)
          {temp <- loess(as.numeric(data1[data.selec,,l]) ~ as.numeric(x.reg[data.selec,]), span=span)}
          if (work.with.deriv.kma==0)
          {temp <- loess(as.numeric(data0[data.selec,,l]) ~ as.numeric(x.reg[data.selec,]), span=span)}

          temp <- predict(temp, x.out)
          y0 <- rbind(y0, temp)
          
        }
        
        templates <- c(templates, list(y0))
      }
      
    }
    
    if (center.method=='k-medoids')
    {
      # Where all functions are defined 
      # (the common domaine)
      x.com  <- seq(max(x.reg[,1], na.rm=T),min(x.reg[,dim(x.reg)[2]], na.rm=T), length=n.out)
            
      # Compute registered functions
      data0.reg<-array(0,dim=c(n.obs,n.out,r))
      data1.reg<-array(0,dim=c(n.obs,n.out,r))
      for(l in 1:r){
        for(k in 1:n.obs){
          
          if(similarity.method == 'd1.L2' 
             || similarity.method == 'd1.pearson'
             || similarity.method == 'd1.L2.centered'
             || similarity.method == 'd1.pearson.mean'
          ) 
          {data1.reg[k,,l]<-approx(x.reg[k,],data1[k,,l],xout=x.com)$y}
          
          if(similarity.method == 'd0.L2' 
             || similarity.method == 'd0.pearson'
             || similarity.method == 'd0.L2.centered'
             || similarity.method == 'd0.pearson.mean'
          ) 
          {data0.reg[k,,l]<-approx(x.reg[k,],data0[k,,l],xout=x.com)$y}
        }
      }
      
      templates <- NULL
      
      for(k in sort(unique(labels))){
        
        data.selec <- which(labels==k)
        y0 <- NULL
        
          
          distanze<-rep(0,length(data.selec))
          for(i in 1:length(data.selec)){
            for(j in 1:length(data.selec)){
              
              if (work.with.deriv.kma==1)
              {	b <- !is.na(data1.reg[i,,1])
          		a <- !is.na(data1.reg[j,,1])
       	  		sel <- a&b 
              	distanze[i]<-distanze[i]+kma.similarity(x.f = x.com[sel], y1.f = data1.reg[i,sel,],  x.g = x.com[sel], y1.g = data1.reg[j,sel,],  similarity.method = similarity.method, unif.grid = FALSE)
              	}
              if (work.with.deriv.kma==0)
              {
              	b <- !is.na(data0.reg[i,,1])
          		a <- !is.na(data0.reg[j,,1])
       	  		sel <- a&b 
              	distanze[i]<-distanze[i]+kma.similarity( x.f = x.com[sel], y0.f = data0.reg[i,sel,], x.g = x.com[sel], y0.g = data0.reg[j,sel,],  similarity.method = similarity.method, unif.grid = FALSE)
              	}
            }
          }
          
          if(similarity.method == 'd1.L2'  || similarity.method == 'd0.L2'  || similarity.method == 'd1.L2.centered' || similarity.method == 'd0.L2.centered') 
          {m<-which.min(distanze)}
          if(similarity.method == 'd1.pearson' || similarity.method == 'd0.pearson' ) 
          {m<-which.max(distanze)}         
          
         for (l in 1:r)
         {
         	if (work.with.deriv.kma==1)
         	 {temp<-approx(x.reg[data.selec[m],],data1[data.selec[m],,l],xout=x.out)$y}
          	if (work.with.deriv.kma==0)
         	 {temp<-approx(x.reg[data.selec[m],],data0[data.selec[m],,l],xout=x.out)$y}
          
          	y0 <- rbind(y0, temp)
         } 
        templates <- c(templates, list(y0))
      }
    }
    
    } else{
  	
  	index <- index.old
  	labels <- labels.old
  	templates <- templates.old
  	iter <- iter-1
  	still.in <- FALSE
  	warning("The algorithm stops because the last iteration is associated to  a decrease of the total similarity; the results provided are associated to the second to last iteration.")
  
  }

  } # end k-mean alignment cycle
    
  if(iter==iter.max) warning('reached maximum number of iterations, method stops')
  
  
  ### calculate global shift and dilation ###
  ###-------------------------------------###
  
  m.final <- rep(1, n.obs)
  t.final <- rep(0, n.obs)
  
  # composition
  for(i in 1:iter){
    m.final <- m.final*m.list[[i]]
    t.final <- t.final*m.list[[i]] + t.list[[i]]
  }
  
  # normalization of global shifts and dilations
  for(k in 1:n.clust){
    data.selec <- which(labels==k)
    m.final[data.selec] <- m.final[data.selec]/mean(m.final[data.selec])
    t.final[data.selec] <- (t.final[data.selec] - mean(t.final[data.selec]))/mean(m.final[data.selec])
  }
  
  ### final warping of the abscissa ###
  ###-------------------------------###
  x.reg <- x
  
  for(i in 1:n.obs){ x.reg[i,] <- m.final[i]*x[i,] + t.final[i] }
  
  ### final templates ###
  ###-----------------###  
  x.out  <- seq(min(x.reg, na.rm=TRUE),max(x.reg, na.rm=TRUE), length=n.out)
  
  #   here templates final k-means
  if (center.method=='k-means')
  {
    templates <- NULL
    
    for(k in sort(unique(labels))){
      
      data.selec <- which(labels==k)
      y0 <- NULL
      
      for(l in 1:r){
        
        if (work.with.deriv.kma==1)
        {temp <- loess(as.numeric(data1[data.selec,,l]) ~ as.numeric(x.reg[data.selec,]), span=span)}
        
        if (work.with.deriv.kma==0)
        {temp <- loess(as.numeric(data0[data.selec,,l]) ~ as.numeric(x.reg[data.selec,]), span=span)}
        
        temp <- predict(temp, x.out)
        y0 <- rbind(y0, temp)
        
      }
      
      templates <- c(templates, list(y0))
    }
  }
  
  #   here templates final k-medoids
  if (center.method=='k-medoids')
  { 
    remind.centers <- rep(0,length(unique(labels)))
    
    # Where all functions are defined (the common domaine)
    x.com  <- seq(max(x.reg[,1], na.rm=T),min(x.reg[,dim(x.reg)[2]], na.rm=T), length=n.out)
    
    # Compute registered functions
    data0.reg<-array(0,dim=c(n.obs,n.out,r))
    data1.reg<-array(0,dim=c(n.obs,n.out,r))
    for(l in 1:r){
      for(k in 1:n.obs){
        
        if (work.with.deriv.kma==1)
        {data1.reg[k,,l]<-approx(x.reg[k,],data1[k,,l],xout=x.com)$y}
        if (work.with.deriv.kma==0)
        {data0.reg[k,,l]<-approx(x.reg[k,],data0[k,,l],xout=x.com)$y}
      }
    }
    
    templates <- NULL
    
    for(k in sort(unique(labels))){
      
      data.selec <- which(labels==k)
      y0 <- NULL
      distanze<-rep(0,n.obs)
      
      for(i in data.selec){
        for(j in data.selec){
            
            if (work.with.deriv.kma==1)
            {
            	b <- !is.na(data1.reg[i,,1])
          		a <- !is.na(data1.reg[j,,1])
       	  		sel <- a&b 
            	distanze[i]<-distanze[i]+kma.similarity( x.f = x.com[sel], y1.f = data1.reg[i,sel,], x.g = x.com[sel], y1.g = data1.reg[j,sel,],  similarity.method = similarity.method, unif.grid = FALSE)
            	}
            
            if (work.with.deriv.kma==0)
			{
				b <- !is.na(data0.reg[i,,1])
          		a <- !is.na(data0.reg[j,,1])
       	  		sel <- a&b 
            	distanze[i]<-distanze[i]+kma.similarity(  x.f = x.com[sel], y0.f = data0.reg[i,sel,], x.g = x.com[sel], y0.g = data0.reg[j,sel,],    similarity.method = similarity.method, unif.grid = FALSE)
            }

        }
      }        
      
      if(similarity.method == 'd1.pearson'  || similarity.method == 'd0.pearson'
      ) 
      {m<-which.max(distanze[data.selec])}
      
      if(similarity.method == 'd1.L2' || similarity.method == 'd0.L2' || similarity.method == 'd1.L2.centered'   || similarity.method == 'd0.L2.centered'
      )
      {m<-which.min(distanze[data.selec])}
        
      if (work.with.deriv.kma==1)
      {
        for (l in 1:r)
        { 
          temp<-approx(x.reg[data.selec[m],],data1[data.selec[m],,l],xout=x.out)$y          
          y0 <- rbind(y0, temp)          
        }
      }
        
      if (work.with.deriv.kma==0)
      {
        for (l in 1:r)
        {
          temp<-approx(x.reg[data.selec[m],],data0[data.selec[m],,l],xout=x.out)$y
          y0 <- rbind(y0, temp)
        }
      }
        
      remind.centers[k] <- data.selec[m]
      templates <- c(templates, list(y0))
      
    } # end of for(k in sort(unique(labels)))
          
      
  }
  
  
  # Controls if n.clust cluster have been founded.
  # Otherwise gives a warning
  n.clust.final.topass <- unique(labels)
  n.clust.final = length(unique(labels))
  if (n.clust.final<n.clust)
  {
    dif <- n.clust - n.clust.final
    warning(dif, ' empty cluster(s) founded. 
n.clust: ', n.clust,'
n.clust.final: ', n.clust.final)
  }  

  
  
  # here we will work on the aligned abscissa grid
  
  
  ###############################################
  ####  Compute the centers of aligned data  ####
  ####    (from derivatives centers)         ####
  ###############################################
  ##########################
  # compute the  function derivatives centers
  # INTEGR K, optimal constant later if y0!=NULL
  
  # Variables in the code:
  # - x.out : abscissas of centers of function first derivatives
  # - templates : evaluations of centers of function first derivatives
  # - xt : abscissas of centers of original functions
  # - integr : evaluations of centers of original functions
  
  #x.center.orig <- seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE), length=n.out)
  
  if (work.with.deriv.kma==1)
  {
    if (center.method=='k-means')
    {
      yt <- templates
      xt <- NULL
      integr=NULL
      r.t <- r
      
      if( length(y0.original)!=0 )
      {
        xt <- x.out
        integr <- yt
        
        for (k in sort(unique(labels)))
        {
          ind <- which(!is.na(yt[[k]]))
          ind.NA <- which(is.na(yt[[k]]))
          
          for (i in 1:length(integr[[k]]))
          {
            
            if(is.na(yt[[k]][i])) { integr[[k]][i] <- NA }
            if(!is.na(yt[[k]][i]))
            {
              if (i==1) { integr[[k]][1] <- yt[[k]][1] }
              if (i>1) 
              {
                if (is.na(yt[[k]][i-1]))
                {
                  if(i==ind[1]) {integr[[k]][i]<-yt[[k]][i]
                  }else{
                    last.value <- max(which(ind<i))
                    integr[[k]][i] <- integr[[k]][last.value]
                  }
                }else{
                  # approx integral with rectangle 
                  # (i), (i-1)
                  integr[[k]][i] <- integr[[k]][i-1]+yt[[k]][i]*(xt[i]-xt[i-1])
                }
              }
            }
          }
        }
      
        # Now we have the integrals of centers of function first derivatives,
        # but without knowing the integral constant.
        # So let's find the optimal constant of the integral
        
        # To decide how to compute the quadratic error
        quadratic.error.viakma.similarity <- 0
        quadratic.error.viaL2norm <- 0
        quadratic.error.viaL2norm.w <- 1
        
        ###################### Auxiliary function:
        # Compute the L2 norm between two curves (even 
        # multidimensional curves) without considering 
        # NA's values
        L2norm.w <- function(f1,f2,weight)
        {
          n.obs.t <- 1
          if( length(dim(f1))==3 ){
            
            r.t <- dim(f1)[3]
            n.camp <- dim(f1)[2]
            function1 <- f1
            function2 <- f2
            
          }else{
            r.t <- 1
            if (class(f1)=='numeric') 
            {
              n.camp <- length(f1)
            }else{n.camp <- dim(f1)[2]}        
            
            function1 <- array(0,dim=c(n.obs.t,n.camp,r.t+1))
            function1[,,1:r.t] <- f1
            function1 <- function1[,,1:r.t,drop=FALSE]
            function2 <- array(0,dim=c(n.obs.t,n.camp,r.t+1))
            function2[,,1:r.t] <- f2
            function2 <- function2[,,1:r.t,drop=FALSE]
          }
          
          
          L2distance <- 0
          for(l in 1:r.t)
          {
            L2distance <- L2distance + sum((function1[,,l]-function2[,,l])^2*weigth,na.rm=TRUE)
          }
          L2distance
        }
        
        L2norm <- function(f1,f2)
        {
          n.obs.t <- 1
          if( length(dim(f1))==3 ){
            
            r.t <- dim(f1)[3]
            n.camp <- dim(f1)[2]
            function1 <- f1
            function2 <- f2
            
          }else{
            r.t <- 1
            if (class(f1)=='numeric') 
            {
              n.camp <- length(f1)
            }else{n.camp <- dim(f1)[2]}        
            
            function1 <- array(0,dim=c(n.obs.t,n.camp,r.t+1))
            function1[,,1:r.t] <- f1
            function1 <- function1[,,1:r.t,drop=FALSE]
            function2 <- array(0,dim=c(n.obs.t,n.camp,r.t+1))
            function2[,,1:r.t] <- f2
            function2 <- function2[,,1:r.t,drop=FALSE]
          }
          
          
          L2distance <- 0
          for(l in 1:r.t)
          {
            L2distance <- L2distance + sum((function1[,,l]-function2[,,l])^2,na.rm=TRUE)
          }
          L2distance
        }
        
        ######################
        labels.final <- labels  
        
        # Approx y0 (original data) on the same abscissa than the integral
        # i.e. xt
        y0.appr <- array( 0,dim=c( dim(data0)[1], length(xt), r.t ) )
        for (h in 1:n.obs)
        {
          for (l in 1:r.t)
          {
            y0.appr[h,,l] <- approx(x.reg[h,],data0[h,,l], xout=xt)$y    
          }
        }
        
        # Create weigths to "weigth" the quadratic error later
        weigth <- array( 0,dim=c( 1, length(xt), r.t ))
        
        for (h in 1:n.obs)
        {
          for (l in 1:r.t)
          {
            for (g in 1:length(xt))
            {
              
              if (!is.na(y0.appr[h,g,l]))
              {
                weigth[1,g,l] <- weigth[1,g,l] +1
              }
              
            }
          }   
        }
        # normalize the weigths
        #weigth <- weigth / sum(weigth)
        
        # If we want to compute the wadratic error via kma.similarity
        # we have to weigth first the functions in y0.appr
        if (quadratic.error.viakma.similarity==1)
        {
          for (h in 1:n.obs)
          {
            for (l in 1:r.t)
            {
              for (g in 1:length(xt))
              {
                
                if (!is.na(y0.appr[h,g,l]))
                {
                  y0.appr[h,g,l] <- y0.appr[h,g,l]*weigth[1,g,l]
                }
                
              }
            }   
          }
        }
        
        # Define function which compute the quadratic error
        # It will be minimized later in order to find the right
        # integral constant
        quad.error <- function(const)
        {
          
          integr.group <- integr[[group]] + const
          
          y0.appr.group <- array( 0,dim=c( length(which(labels.final==group)), length(xt), r.t ) )
          for (l in r.t)
          {
            y0.appr.group[,,l] <- y0.appr[which(labels.final==group),,l]
          }
          
          distance.sum <- 0
          for (dd in 1:length(which(labels.final==group)))
          {
            ##### Without weights (L2norm)
            if (quadratic.error.viaL2norm==1)
            {distance.sum <- distance.sum + L2norm (integr.group,y0.appr.group[dd,,])}
            
            ##### With weigths (L2norm.w) OK!!!
            if (quadratic.error.viaL2norm.w==1)
            {distance.sum <- distance.sum + L2norm.w (integr.group,y0.appr.group[dd,,],weigth)}
            
            ##### With weigths (kma.similarity)
            if (quadratic.error.viakma.similarity==1)
            {
              if (work.with.deriv.kma==1)
              {distance.sum <- distance.sum + kma.similarity(
                x.f = xt, y1.f = integr.group,
                x.g = xt, y1.g = y0.appr.group[dd,,],
                similarity.method = similarity.method, unif.grid = FALSE)}
              
              if (work.with.deriv.kma==0)
              {distance.sum <- distance.sum + kma.similarity(
                x.f = xt, y0.f = integr.group,
                x.g = xt, y0.g = y0.appr.group[dd,,],
                similarity.method = similarity.method, unif.grid = FALSE)}
              
              #distance.sum <- distance.sum + kma.similarity (
              #integr.group,y0.appr.group[dd,,],similarity.method)
            }
          }
          distance.sum
        }
        
        # Using optim function we can find the constant minimizing
        # the weigthed quadratic error (with respect 
        # to the relative cluster) and update the centers
        constant.opt <- rep(0, n.clust.final)
        for (group in 1:n.clust.final)
        {
          result <- optim(0, quad.error, method='BFGS')
          constant.opt[group] <- result$par
          # Update the centers of original data
          integr[[group]] <- integr[[group]] + constant.opt[group]
        }
                
        ###############################################
        ####     End centers of aligned data       ####
        ###############################################
      }
      
    } # end if (center.method=='k-means')
    
    if (center.method=='k-medoids')
    {
      # In remind.centers we have centers indexes saved before
      # we don't need to integrate. We have the indexes of 
      # function first derivatives centers saved in 
      # remind.centers. So we simply take the original 
      # functions as centers (those corresponding to 
      # the indexes in remind.centers).
      
      # centers of function first derivates 
      x.out <- x.out
      templates <- templates
      # centers of original functions 
      xt <- NULL
      integr <- NULL
      
      if( length(y0.original)!=0 )
      {
        # centers of original functions approximated 
        # on x.out to be coherent with other centers
        integr <- templates        
        for (k in sort(unique (labels)))
        {          
          #integr[[k]] <- data0[remind.centers[k],,]
          integr[[k]] <- approx(x.reg[remind.centers[k],],data0[remind.centers[k],,],xout=x.out)$y
        }
        xt <- x.out
      } 
    } # end if (center.method=='k-medoids') in case:'d1.pearson'
    
  } # end in case:'d1.pearson' or 'd1.L2'
  
  if (work.with.deriv.kma==0)
  {
    # we already worked on original functions, but we have all saved 
    # in the variables related to the function first derivatives
    # so, simply, we "swap" the values
    integr <- templates
    templates <- NULL
    xt <- x.out
    #x.out <- NULL
  } # end in case:'d0.L2' or 'd0.pearson'
  
  
  # here we will work on the original abscissa grid
  
  # If y0 is given as input, let's compute the center 
  # of orginal data.
  # If center.method='k.medoids' it will be the 
  # corrispondent original function, whose index has been 
  # saved in remind.center.orig
  # If center.method='k.means' we have to compute the 
  # integral of y1.center.orig on the same abscissa
  
  if (length(y0.original)!=0)
  {
    if (work.with.deriv.kma==1)
    {
      if (center.method=='k-medoids')
      {
        # In remind.center.orig there is the index of the 
        # medoid of original function first derivatives
        temp <- NULL
        y0 <- NULL
        
        for(l in 1:r)
        {
          # Here we approx the medoid of original 
          # similarity on x.out                    
          temp<-approx(x[remind.center.orig,],data0[remind.center.orig,,l],xout=x.center.orig)$y
          y0 <- rbind(y0, temp)
        }
        y0.center.orig <- c(y0.center.orig, list(y0))
      }
      
      if (center.method=='k-means')
      {
        # We have to integrate the loess of function first 
        # derivatives to obtain the loess of original data
        
        # Here we have the loess of function first derivatives:
        # y1.center.orig defined on x.out
        integral_list <- function(list_in.x ,list_in.y)
        {
          if( length(list_in.y)!=0 ) {list_out.y <- list_in.y}
          
          if( length(list_in.y)!=0 )
          {
            list_out.x <- list_in.x
            list_out.y <- list_in.y
            
            for (k in 1:length(list_in.y))
            {
              ind <- which(!is.na(list_in.y[[k]]))
              ind.NA <- which(is.na(list_in.y[[k]]))
              
              for (i in 1:length(list_in.y[[k]]))
              {
                
                if(is.na(list_in.y[[k]][i])) { list_out.y[[k]][i] <- NA }
                if(!is.na(list_in.y[[k]][i]))
                {
                  if (i==1) { list_out.y[[k]][1] <- list_in.y[[k]][1] }
                  if (i>1) 
                  {
                    if (is.na(list_in.y[[k]][i-1]))
                    {
                      if(i==ind[1]) {list_out.y[[k]][i]<-list_in.y[[k]][i]
                      }else{
                        last.value <- max(which(ind<i))
                        list_out.y[[k]][i] <- list_out.y[[k]][last.value]
                      }
                    }else{
                      # approx integral with rectangle 
                      # (i), (i-1)
                      list_out.y[[k]][i] <- list_out.y[[k]][i-1]+list_in.y[[k]][i]*(list_in.x[i]-list_in.x[i-1])
                    }
                  }
                }
              }
            }
          }
          
          list_out.y
          
        } # end integral_list
        
        y0.center.orig <- integral_list(x.center.orig, y1.center.orig)
        
        
        # Here we find the optimal integral constant 
        
        # Approx y0 (original data) on the same abscissa than the integral
        # i.e. xt
        y0.appr <- array( 0,dim=c( n.obs, length(x.center.orig), r ) )
        for (h in 1:n.obs)
        {
          for (l in 1:r)
          {
            y0.appr[h,,l] <- approx(x[h,],data0[h,,l], xout=x.center.orig)$y    
          }
        }
        
        # Create weigths to "weigth" the quadratic error later
        weigth <- array( 0,dim=c( 1, length(x.center.orig), r ))
        
        for (h in 1:n.obs)
        {
          for (l in 1:r)
          {
            for (g in 1:length(x.center.orig))
            {
              
              if (!is.na(y0.appr[h,g,l]))
              {
                weigth[1,g,l] <- weigth[1,g,l] +1
              }
              
            }
          }   
        }
        
        
        ##########################################
        quad.error.y0.center.orig <- function(const)
        {
          
          y0.center.orig[[1]] <- y0.center.orig[[1]] + const
          
          distance.sum <- 0
          for (dd in 1:n.obs)
          {
            ##### With weigths (L2norm.w) OK!!!
            distance.sum <- distance.sum + L2norm.w (y0.center.orig[[1]],y0.appr[dd,,],weigth)
          }
          distance.sum
        }
        
        # Using optim function we can find the constant minimizing
        # the weigthed quadratic error (with respect 
        # to the relative cluster) and update the centers
        constant.opt <- rep(0, 1)
        for (group in 1:1)
        {
          result <- optim(0, quad.error.y0.center.orig, method='BFGS')
          constant.opt[group] <- result$par
          # Update the centers of original data
          y0.center.orig[[group]] <- y0.center.orig[[group]] + constant.opt[group]
        }
        
        ##########################################
        
        
      }
      
    }
  }
  
  
  
  # Convert centers form list to matrix
  # if data are 2-D
  if (r==1)
  {
    
    list_to_matrix <- function (list.in)
    {
      if (length(list.in)==0) {out <- list.in}
      if (length(list.in)!=0)
      {
        list.in.x <- length(list.in)
        list.in.y <- length(list.in[[1]])
        matrix.temp <- matrix(NA,list.in.x,list.in.y)
        for (i in 1:list.in.x)
        {
          for (j in 1:list.in.y)
          {
            matrix.temp[i,j] <- list.in[[i]][j]
          }
        }
        out <- matrix.temp
      }
      
      out
    }
    
    #y.center.orig <- list_to_matrix(y.center.orig)
    y0.center.orig <- list_to_matrix(y0.center.orig)
    y1.center.orig <- list_to_matrix(y1.center.orig)
    
    integr <- list_to_matrix(integr)
    templates <- list_to_matrix(templates)
    
  }
  
  
  total_result <- list(
        
    iterations = iter,
    x = x, # cambio x.orig=orig.absc
    y0 = y0.original, # aggiunto, prima non c'era
    y1 = y1.original, # cambio y1.orig=orig.data
    n.clust = n.clust.original,
    warping.method = warping.method, # cambio warping.method=Registration  
    similarity.method = similarity.method.original,
    center.method = center.method,
    
    x.center.orig = x.center.orig, # cambio x.center.orig=x.templates.iter1
    y0.center.orig = y0.center.orig, # cambio y.center.orig=y.templates.iter1
    #x.center.orig = x.center.orig,
    y1.center.orig = y1.center.orig,
    similarity.orig = similarity.orig,
    
    x.final = x.reg, # cambio x.reg=reg.absc
    n.clust.final = n.clust.final, # cambio n.clust.final=n.clust
    
    #x0.centers.final = xt,
    #x1.centers.final = x.out, # 
    x.centers.final = x.out, # 
    y0.centers.final = integr,
    y1.centers.final = templates, # 
    
    labels = labels, # cambio labels.final=labels
    similarity.final = index, # cambio similarity.final=similarity
    
    dilation.list = m.list,
    shift.list = t.list,
    dilation = m.final,
    shift = t.final)
   
  OUT <- c(OUT, list(total_result))  
}
if (return.all==TRUE)
{
	return(OUT)
}else
{
	mean_sim <- rep(0, nstart)
	for (i in 1:nstart)
	{
		mean_sim[i] <- mean(OUT[[i]]$similarity.final)
		if (similarity.method == 'd1.pearson' || similarity.method == 'd0.pearson') 
		{
			j <- which.max(mean_sim)
		}
		if (similarity.method == 'd1.L2'  || similarity.method == 'd0.L2' || similarity.method == 'd1.L2.centered' || similarity.method == 'd0.L2.centered')
		{
			j <- which.min(mean_sim)
		}
	}
	return(OUT[[j]])
}

}
