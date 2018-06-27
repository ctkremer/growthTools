#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Growth rate estimation routines:   ####

# Developed by CTK for NSF Dimensions project

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#' Helper function for smoothed lagged/saturating abundance equations
#' 
sqfunc<-function(x,b,s){
  (1/2)*sqrt(b*(4*s+b*x^2))
}


#' Equations for modeling abundace time series
#' 
#' Intended to allow the extraction of exponential growth rates from time series 
#' while accounting for the presence of initial lags in growth, saturating abundances,
#' or both in the same time series. These equations provide smoothed piecewise linear
#' functions, where lagged or saturated portions of the time series maintain constant 
#' abundance, and elsewhere abundance increases linearly.
#' 
#' This approach for lag, saturation, and lag+saturation are based on:
#' https://stats.stackexchange.com/questions/149627/piecewise-regression-with-constraints
#' which invokes a smooth approximation to a piecewise linear function,
#' where parameter s determines the smoothness around break-points. Generally, as s->0 
#' this smooth model approximates more closely the piecewise linear one. The s term 
#' could be fit explicitly, but for now it is fixed at a small number (1E-10).
#' 
#' Note: Currently, only the linear model without lag or saturation can produce negative
#' growth rate estimates. The lagged/saturating models will be extended to allow this 
#' possibility in future versions of this package.
#' 
#' @param x Time variable
#' @param a Initial abundance at time = 0
#' @param logb ln(slope) of the increasing linear portion of the time series
#' @param B1 Time point where abundance starts to increase (leaves lag phase)
#' @param B2 Time point where abundance stops increasing (saturates)
#' @param s Smoothing parameter; as this term -> 0, these continuous functions approach true piecewise equations
#' 
#' @return Abundance at time x as a function of model parameters
#' 
#' @examples 
#' 
#' curve(lag(x,5,0,4,s=1E-10),0,10,col='green',ylim=c(0,11),ylab='Abundance')
#' curve(sat(x,0.9,0,8,s=1E-10),0,10,col='red',add=T)
#' curve(lagsat(x,5.1,0,4,8,s=1E-10),0,10,col='blue',add=T)
#' 
#' @export
lagsat<-function(x,a,logb,B1,B2,s=1E-10){
  a + (1/2)*exp(logb)*(B2-B1) + sqfunc(B1-x,exp(logb),s) - sqfunc(B2-x,exp(logb),s)
}

#' @describeIn lagsat Lagged increasing linear function
#' @export
lag<-function(x,a,logb,B1,s=1E-10){
  sqfunc(B1-x,exp(logb),s)-(exp(logb)/2)*(B1-x)+a
}

#' @describeIn lagsat Saturating linear function
#' @export
sat<-function(x,a,logb,B2,s=1E-10){
  a + (1/2)*exp(logb)*(B2) + sqfunc(-x,exp(logb),s) - sqfunc(B2-x,exp(logb),s)
}

#' Extract exponential growth rate assuming exponential growth
#' 
#' This function fits a linear model to ln(abundance) data.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a linear model regressing ln(abundance) on time
#' 
#' @examples 
#' 
#' @export
get.gr<-function(x,y,plotQ=F,fpath=NA,id=''){
  lm1<-lm(y~x)
  
  if(plotQ){
    if(!is.na(fpath)){
      pdf(fpath)
      plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      abline(lm1,col='red')
      dev.off()
    }else{
      plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      abline(lm1,col='red')
    }
  }
  
  return(lm1)
}

#' Extract exponential growth rate assuming lagged exponential growth
#' 
#' This function fits a smoothed piecewise linear model to ln(abundance) data, with 
#' the assumption that abundances are nearly constant for several time points, before 
#' exponential growth kicks in.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a nonlinear least-squares regression model
#' 
#' @examples 
#' 
#' @export
#' @import minpack.lm
#' @import zoo
get.gr.lag<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  #slopes <- rollapply(data, 3, localslope, by.column=F)
  
  fit.lag<-nlsLM(y ~ lag(x,a,logb,B1,s=1E-10),
                 start = c(B1=mean(x)-(mean(x)-min(x))/2, a=min(y), logb=0),data = data,
                 control = nls.control(maxiter=1000, warnOnly=TRUE))
  cfs<-data.frame(t(coef(fit.lag)))
  
  #start=c(B1=4,B2=8,a=1,logb=0),
  #start = c(B1=mean(x), a=min(y)+1E-10, logb=log(max(slopes))),
  #start = c(B1=min(x)+1E-10, B2=max(x)-1E-10, a=min(y)+1E-10, logb=log(max(slopes))),
  
  if(plotQ){
    if(!is.na(fpath)){
      pdf(fpath)
      plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      curve(lag(x,cfs$a,cfs$logb,cfs$B1,s=1E-10),min(x),max(x),add=T,col='red')
      dev.off()
    }else{
      plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      curve(lag(x,cfs$a,cfs$logb,cfs$B1,s=1E-10),min(x),max(x),add=T,col='red')
    }
  }
  return(fit.lag)
}

#' Extract exponential growth rate assuming exponential growth that saturates
#' 
#' This function fits a smoothed piecewise linear model to ln(abundance) data, with 
#' the assumption that abundances increase linearly at first, but then saturate and
#' remain constant.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a nonlinear least-squares regression model
#' 
#' @examples 
#' 
#' @export
#' @import minpack.lm
#' @import zoo
get.gr.sat<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  slopes <- rollapply(data.frame(x=x,y=y), 3, localslope, by.column=F)
  a.guess<-coef(lm(y~x))[[1]]
  
  fit.sat<-nlsLM(y ~ sat(x,a,logb,B1,s=1E-10),
                 #start = c(B1=mean(x)+(max(x)-mean(x))/2, a=1.5*a.guess, logb=-0.1),
                 #start=c(B1=9,a=-0.2933,logb=0),
                 start=c(B1=mean(x)+(max(x)-mean(x))/2,a=a.guess,logb=round(log(max(slopes)),5)),
                 data = data,
                 control = nls.control(maxiter=1000, warnOnly=TRUE))
  cfs<-data.frame(t(coef(fit.sat)))
  
  if(plotQ){
    if(!is.na(fpath)){
      pdf(fpath)
      plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      curve(sat(x,cfs$a,cfs$logb,cfs$B1,s=1E-10),min(x),max(x),add=T,col='red')
      dev.off()
    }else{
      plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      curve(sat(x,cfs$a,cfs$logb,cfs$B1,s=1E-10),min(x),max(x),add=T,col='red')
    }
  }
  return(fit.sat)
}

#' Extract exponential growth rate assuming lagged exponential growth that saturates
#' 
#' This function fits a smoothed piecewise linear model to ln(abundance) data, with 
#' the assumption that abundances are nearly constant for several time points, before 
#' exponential growth kicks in; subsequently, growth saturates and abundances become 
#' constant again.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plotQ logical; should the fit be plotted?
#' @param fpath character; path specifying where plot should be saved, if generated
#' @param id Label corresponding to the population/strain/species of interest
#' 
#' @return This function returns a nonlinear least-squares regression model
#' 
#' @examples 
#' 
#' @export
#' @import minpack.lm
get.gr.lagsat<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  #slopes <- rollapply(data.frame(x=x,y=y), 3, localslope, by.column=F)
  #log(max(slopes))
  #a.guess<-coef(lm(y~x))[[1]]
  # round(log(max(slopes)),5)
  
  fit.lagsat<-nlsLM(y ~ lagsat(x,a,logb,B1,B2,s=1E-10),
                    start = c(B1=mean(x)-(mean(x)-min(x))/2,B2=mean(x)+(max(x)-mean(x))/2, a=min(y)+0.1, logb=0),
                    data = data,
                    control = nls.control(maxiter=1000, warnOnly=TRUE))
  cfs<-data.frame(t(coef(fit.lagsat)))
  
  if(plotQ){
    if(!is.na(fpath)){
      pdf(fpath)
      plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      curve(lagsat(x,cfs$a,cfs$logb,cfs$B1,cfs$B2,s=1E-10),min(x),max(x),add=T,col='red')
      dev.off()
    }else{
      plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
      curve(lagsat(x,cfs$a,cfs$logb,cfs$B1,cfs$B2,s=1E-10),min(x),max(x),add=T,col='red')
    }
  }
  return(fit.lagsat)
}

#' Local Slope function
#' 
#' Helper function to calculate and extract the slope of a basic linear regression 
#' relating y to x; the resulting value is used to obtain a reasonable starting guess
#' for the slopes of the piecewise linear functions in \code{lag}, \code{sat}, and
#'  \code{lagsat}
#' 
#' @param d A data frame containing two columns, x and y
#' 
#' @return Slope of the linear regression
#' 
#' @export
localslope<-function (d) {
  m <- lm(y~x, as.data.frame(d))
  return(coef(m)[2])
}

#' Detect model failure
#' 
#' Helper function used to determine whether an attempt at an mle2 fit failed
#' 
#' @param x An object
#' 
#' @return TRUE/FALSE depending on whether or not 'try-error' is listed as (one) of the classes of object x.
detect<-function(x){
  !(c('try-error') %in% class(x))
}


#' Extract exponential growth rate from a time series of ln(population abundance)
#' 
#' This meta-function takes a time series of abundance, and attempts to extract an 
#' estimate of exponential growth rate, using one or more of a suite of possible methods.
#' These methods allow for the possibility that exponential growth may lag or saturate,
#' or both, over the course of the time series. All selected methods are used to fit 
#' models to the time series. Subsequently, model comparison (based on AIC) is used to 
#' determine which model best fits the focal data.
#' 
#' @param x Time steps
#' @param y ln(abundance)
#' @param plot.best.Q logical; should the best fitting model be plotted?
#' @param fpath character; if best model is to be plotted, provide the file path for saving the plot
#' @param methods Must be a character vector containing one or more of \code{'linear'}, \code{'lag'}, \code{'sat'}, or \code{'lagsat'}
#' @param id Label corresponding to the population/strain/species of interest; used to determine the title and file name of saved plot, if any.
#' 
#' @return A data frame containing the identity of the best model, the content of the best model, the estimated slopes of the increasing linear portion of the regressions (ie, exponential growth rate), and the full list of all models fit.
#' 
#' 
#' @examples 
#' 
#' @export
#' @import bbmle
get.growth.rate<-function(x,y,id,plot.best.Q=F,fpath=NA,methods=c('linear','lag','sat','lagsat')){
  
  if(sum(methods %in% c('linear','lag','sat','lagsat'))==0){
    print('Error! None of the specified methods matched a currently implemented approach')
  }
  
  # Initialize empty data structures
  modlist<-list(gr=NA,gr.lag=NA,gr.sat=NA,gr.lagsat=NA)
  class(modlist$gr)<-class(modlist$gr.lag)<-class(modlist$gr.sat)<-class(modlist$gr.lagsat)<-'try-error'
  gr<-gr.lag<-gr.sat<-gr.lagsat<-NA
  slope.gr<-slope.gr.lag<-slope.gr.sat<-slope.gr.lagsat<-NA
  
  # Fill data structures for each given method, if requested by user:
  if('linear' %in% methods){
    gr<-get.gr(x,y) 
    slope.gr<-coef(gr)[[2]]
    modlist$gr<-gr
  }
  if('lag' %in% methods){
    gr.lag<-try(get.gr.lag(x,y))  
    slope.gr.lag<-ifelse(prod(class(gr.lag)!='try-error'),
                         exp(coef(gr.lag)[3]),NA)
    modlist$gr.lag<-gr.lag
  }
  if('sat' %in% methods){
    gr.sat<-try(get.gr.sat(x,y))
    slope.gr.sat<-ifelse(prod(class(gr.sat)!='try-error'),
                         exp(coef(gr.sat)[3]),NA)
    modlist$gr.sat<-gr.sat
  }
  if('lagsat' %in% methods){
    gr.lagsat<-try(get.gr.lagsat(x,y))
    slope.gr.lagsat<-ifelse(prod(class(gr.lagsat)!='try-error'),
                            exp(coef(gr.lagsat)[4]),NA)
    modlist$gr.lagsat<-gr.lagsat
  }
  
  # determine which fits occured and were successful
  successful.fits<-sapply(modlist,detect)
  
  if(sum(successful.fits)==0){
    print('Error! All results for requested methods failed!')
    break();
  }
  
  # assemble model names, contents, and slopes, but only for successful fits
  mod.names<-c('gr','gr.lag','gr.sat','gr.lagsat')[successful.fits]
  mod.list<-list(gr,gr.lag,gr.sat,gr.lagsat)[successful.fits]
  slope.ests<-c(slope.gr,slope.gr.lag,slope.gr.sat,slope.gr.lagsat)[successful.fits]
  
  # compare successful models
  aictab<-AICtab(mod.list,mnames = mod.names)
  best.mod.id<-which(mod.names==attr(aictab,"row.names")[1])
  
  # format output:
  result<-list(best.slope=slope.ests[[best.mod.id]],
               best.model=as.character(mod.names[[best.mod.id]]),
               best.model.contents=list(mod.list[[best.mod.id]]),
               slopes=slope.ests,
               models=list(gr=gr,gr.lag=gr.lag,gr.sat=gr.sat,gr.lagsat=gr.lagsat))
  #print(result)
  
  if(plot.best.Q){
    if(!is.na(fpath)){
      fpath<-paste(fpath,id[1],'.pdf',sep='')
    }
    
    # want to show the best model in the requested model set... given methods options.
    # no model can end up in the model set if not requested, so this should be o.k. as written
    gigo<-switch(result$best.model,
                 gr=get.gr(x,y,plotQ=T,fpath=fpath,id=id[1]),
                 gr.lag=get.gr.lag(x,y,plotQ=T,fpath=fpath,id=id[1]),
                 gr.sat=get.gr.sat(x,y,plotQ=T,fpath=fpath,id=id[1]),
                 gr.lagsat=get.gr.lagsat(x,y,plotQ=T,fpath=fpath,id=id[1]))
  }
  
  return(result)
}


# These two functions help to trim data sets to exclude temperature treatments that are well beyond a species' thermal niche, by identifying the minimum and maximum temperatures where positive growth was observed, and determining the next highest (lowest) temperature treatment.
funky.low<-function(x,y){
  dat<-data.frame(x,y) 
  dat2<-dat %>% group_by(x) %>% summarise(max.y=max(y))
  Ts<-sort(unique(dat2$x))
  
  minT<-min(dat2$x[dat2$max.y>0])
  res<-Ts[max(c(1,which(Ts==minT)-1))]
  
  return(res)
}

funky.high<-function(x,y){
  dat<-data.frame(x,y) 
  dat2<-dat %>% group_by(x) %>% summarise(max.y=max(y))
  Ts<-sort(unique(dat2$x))
  
  maxT<-max(dat2$x[dat2$max.y>0])
  res<-Ts[min(c(length(Ts),which(Ts==maxT)+1))]
  
  return(res)
}