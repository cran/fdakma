kma.show.results <-
function (Result)
{

  if (length(Result)==0)
  {
    stop('First parameter (Result) must not be NULL')
  }

  spessore.centers <- 3
  #############################################################
  ######################     0       ##########################
  ## Handling data
  
  ############### Output ok kma function
  iterations <- Result$iterations 
  x <- as.matrix(Result$x)
  y0 <- Result$y0
  y1 <- Result$y1 
  n.clust <- Result$n.clust
  warping.method <- Result$warping.method 
  similarity.method <- Result$similarity.method 
  center.method <- Result$center.method
  x.center.orig <- Result$x.center.orig
  y0.center.orig <- Result$y0.center.orig
  y1.center.orig <- Result$y1.center.orig
  
  similarity.orig <- Result$similarity.orig
  x.final <- Result$x.final
  n.clust.final <- Result$n.clust.final
  x.centers.final <- Result$x.centers.final 
  y1.centers.final <- Result$y1.centers.final #[[2]] 
  x.centers.final <- Result$x.centers.final 
  y0.centers.final <- Result$y0.centers.final #[[2]] 
  labels.final <- Result$labels #labels <- Result$labels 
  similarity.final <- Result$similarity.final 
  dilation.list <- Result$dilation.list
  shift.list <- Result$shift.list 
  dilation <- Result$dilation 
  shift <- Result$shift 
  ##############
  
  if( length(dim(y1))!=0 )
  {
    n.camp <- dim(y1)[2]
    n.obs <- dim(y1)[1] 
    
    if( length(dim(y1))==3 ){
      r <- dim(y1)[3]
    }else{
      r <- 1    
    }
  }
  
  if( length(dim(y0))!=0 )
  {
    n.camp <- dim(y0)[2]
    n.obs <- dim(y0)[1] 
    
    if( length(dim(y0))==3 ){
      r <- dim(y0)[3]
    }else{
      r <- 1    
    }
  }
  
  

  # If x is a vector, apply this abscissa to all curves  
#   if(dim(x)[1]==1 && dim(y1)[1]!=1) {
#     x.temp <- x
#     for(i in 1:(dim(y1)[1]-1)) 
#       x <- rbind(x,x.temp)
#   }  
  
  sim.final <- mean(similarity.final)
  
  # Colours choice
  labels.unique <- sort(unique(labels.final))
  myrainbow <- c('red','blue','green3','orange','grey','yellow')
  myrainbow.dark<-c('darkred','darkblue','darkgreen','darkorange','black','brown')
  myrainbow <- c( myrainbow, rainbow(length(labels.unique) ) )
  myrainbow <- myrainbow[1:n.clust]
  myrainbow.dark <- c( myrainbow.dark, rainbow(length(labels.unique) ) )
  myrainbow.dark <- myrainbow.dark[1:n.clust]
    
  # Colors curves
  colours.random <- rainbow(n.obs)
  colours.bygroup <- rep(0,n.obs)
  
  for (k in labels.unique)
  {
    colours.bygroup[ which(labels.final==k) ] <- myrainbow[k] 
    colours.bygroup[ which(labels.final==k) ] <- myrainbow[k] 
  }
  
  colours.bygroup.dark <- rep(0,n.obs)
  for (k in labels.unique)
  {
    colours.bygroup.dark[ which(labels.final==k) ] <- myrainbow.dark[k] 
    colours.bygroup.dark[ which(labels.final==k) ] <- myrainbow.dark[k] 
  }
  
  # Colors warping functions
  colours.warping <- rep(0,n.obs)
  for (k in labels.unique)
  {
    colours.warping[ which(labels.final==k) ] <- myrainbow[k] 
    colours.warping[ which(labels.final==k) ] <- myrainbow[k] 
  }
  
  colours.templates.iter1 <- myrainbow
  colours.templates.last <- myrainbow
  
  
  #############################################################
  ######################     1 Data       #####################
  ## Plot original DATA + templates of first iteration
  ## (template given by loess/medoid applied to all original derivative)
  
  if(length(y0)!=0)
  {
    dev.new()
    
    matplot(t(x),t(y0),type='l', col=colours.random, xlab="x", ylab='y')
    title(main="Original Data")
    
    for (k in 1:dim(y0.center.orig)[1])    
    {
      lines(x.center.orig,
            #y.center.orig[[k]],
            y0.center.orig[k,],
            lwd=4, 
            col= colours.templates.iter1[k]) 
    }
    
  }
  
  #############################################################
  ######################     2 Data       #####################
  ## Plot aligned DATA coloured by clusters + templates
  
  if(length(y0)!=0)
  {
    dev.new()
    
    tex <- paste(c('k = ', n.clust), collapse='')
    matplot(t(x.final),t(y0),type='l', col=colours.bygroup.dark, xlab="x", ylab=tex)
    
    title2 <- c("Registration: ", warping.method)
    title2 <- paste(title2, collapse='')
    title22 <- c('Aligned Data')
    title2def <- c(title2, title22)
    title(main=title2def)  
    
    # Draw the centers of data aligned 
    for (k in 1:n.clust.final)
    {
      
      lines(x.centers.final,
            #y0.centers.final[[k]],
            y0.centers.final[k,],
            lwd=spessore.centers, 
            # qui spessore centers data
            col= colours.templates.iter1[k]) 
      
    }
    
    # Legend by group
    text <- rep(0,length(labels.unique))
    for (i in 1:length(labels.unique))
    {
      mamma <- c('Cluster ', i)
      mamma <- paste(mamma, collapse="")
      text[i] <- mamma
    }
    lty <- rep(1,length(text))
    
    legend('topleft', 
           legend=text, 
           col=colours.templates.last, 
           lty=lty, 
           cex = 0.6)
    
  }
  
  
  #############################################################
  ######################     1 Derivative      ################
  ## Plot original derivative + templates of first iteration
  ## (template given by loess applied to all original derivative)
        
  if(length(y1)!=0 
     && (similarity.method=='d1.pearson' 
         || similarity.method=='d1.L2'
         || similarity.method=='d1.L2.centered')
     )
  {
    dev.new()
    
    matplot(t(x),t(y1),type='l', col=colours.random, xlab="x", ylab="y")
    #title(main="Original Data", sub=Result$warping.method)
    title(main="Original Derivative")
    
    #for (k in 1:nrow(y.center.orig[[1]]))
    for (k in 1:dim(y1.center.orig)[1])    
    {
      lines(x.center.orig,
            #y.center.orig[[k]],
            y1.center.orig[k,],
            lwd=4, 
            col= colours.templates.iter1[k]) 
    }
  }

  
  #############################################################
  ######################     2 Derivative      ################
  ## Plot aligned derivative coloured by clusters + templates
  
  if(length(y1)!=0 
     && (similarity.method=='d1.pearson' 
         || similarity.method=='d1.L2'
         || similarity.method=='d1.L2.centered')
    )
  {
    dev.new()
    
    matplot(t(x.final),t(y1),type='l', col=colours.bygroup.dark, xlab="x", ylab="y")
    
    title2 <- c("Registration: ", Result$warping.method)
    title2 <- paste(title2, collapse='')
    title22 <- c('Aligned Derivative')
    title2def <- c(title2, title22)
    title(main=title2def)  
    
    #title(main="Aligned Data", sub=Result$warping.method ) 
    
    for ( k in 1:n.clust.final)
    {
      lines(Result$x.centers.final,
            #t(Result$y1.centers.final[[k]]),
            #Result$y1.centers.final[[k]],
            Result$y1.centers.final[k,],
            lwd=spessore.centers, 
            # qui spessore centers derivatives
            col= colours.templates.last[k]) 
      
    }
    
    # Legend by group
    text <- rep(0,length(labels.unique))
    for (i in 1:length(labels.unique))
    {
      mamma <- c('Cluster ', i)
      mamma <- paste(mamma, collapse="")
      text[i] <- mamma
    }
    lty <- rep(1,length(text))
    
    legend('topleft', 
           legend=text, 
           col=colours.templates.last, 
           lty=lty, 
           cex = 0.6)
    
  }
  
  
  #############################################################
  ######################     3       ##########################
  ## Plot warping functions
  
  dev.new()
  
  plot(t(x[1,]),t(x.final[1,]),xlim=c(min(x,na.rm = TRUE),max(x,na.rm = TRUE)), 
       ylim=c(min(x,na.rm = TRUE),max(x,na.rm = TRUE)), type='l', lwd=1, col=colours.warping[1], 
       xlab="x", ylab="y",asp=1)
  for (i in 2:n.obs)
  {
#     lines(t(x[i,]),t(x.final[i,]),type='l', lwd=1, col=colours.bygroup[i],
#           xlab="x", ylab="y",asp=1)
    lines(t(x[i,]),t(x.final[i,]),type='l', lwd=1, col=colours.warping[i],
          xlab="x", ylab="y",asp=1)
  }
  
  # Legend warping functions
  #title2 <- c("Warping functions", " - ", "Registration: ", Result$warping.method)
  title3 <- c("Registration: ", Result$warping.method)
  title3 <- paste(title3, collapse='')
  title33 <- c('Warping Functions')
  title3def <- c(title3, title33)
  title(main=title3def)  
  
  abline(v=min(x))
  abline(v=max(x))
  
  text <- rep(0,length(labels.unique))
  for (i in 1:length(labels.unique))
  {
    mamma <- c('Cluster ', i)
    mamma <- paste(mamma, collapse="")
    text[i] <- mamma
  }
  lty <- rep(1,length(text))
  
  legend(
        x=min(x),
        #y=max(x)+0.20,
        y=max(x),
        #'topleft',
         legend=text, 
         col=colours.templates.last, 
         lty=lty, 
         cex = 0.6)
  
  #############################################################
  ######################     4       ##########################
  ## Boxplots: - similarities of original data with respect 
  ##             to first templates
  
  
  sim.orig <- Result$similarity.orig
  sim.final <- Result$similarity.final
  
  etichette <- rep(0,2)
  etichette[1] <- 'orig. data'
  etichette[2] <- paste('k =',n.clust)          
  
  dev.new()
  
  if (similarity.method=='d1.pearson' || similarity.method=='d0.pearson'
      || similarity.method=='d1.pearson.mean' || similarity.method=='d0.pearson.mean')
  {
    boxplot(sim.orig, sim.final, 
          notch=F, 
          boxwex = 0.3, 
          col=c('grey','orange') 
          ,ylim=c( min(sim.orig,sim.final), 1)
          )
  }
  if (similarity.method=='d0.L2' || similarity.method=='d1.L2'
      || similarity.method=='d0.L2.centered' || similarity.method=='d1.L2.centered')
  {
    boxplot(sim.orig, sim.final, 
            notch=FALSE, 
            boxwex = 0.3, 
            col=c('grey','orange') 
            ,ylim=c( min(sim.orig,sim.final), max(sim.orig,sim.final))
    )
  }
  
  
  title4 <- c("Registration: ", Result$warping.method)
  title4 <- paste(title4, collapse='')
  title44 <- c('Boxplot Similarity Indexes')
  title4def <- c(title4, title44)
  title(main=title4def)  
  #title('Similarity Indexes', sub=Result$warping.method) 
  axis(1,at=1:2,labels=etichette,las=0)
      
}
