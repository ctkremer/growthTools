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
#' @param b slope of the increasing linear portion of the time series, must be >=0
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
lagsat<-function(x,a,b,B1,B2,s=1E-10){
  a + (1/2)*b*(B2-B1) + sqfunc(B1-x,b,s) - sqfunc(B2-x,b,s)
}

#' @describeIn lagsat Lagged increasing linear function
#' @export
lag<-function(x,a,b,B1,s=1E-10){
  sqfunc(B1-x,b,s)-(b/2)*(B1-x)+a
}

#' @describeIn lagsat Saturating linear function
#' @export
sat<-function(x,a,b,B2,s=1E-10){
  a + (1/2)*b*(B2) + sqfunc(-x,b,s) - sqfunc(B2-x,b,s)
}

#' @describeIn lagsat Floored decreasing linear function
#' @export
flr<-function(x,a,b,B2,s=1E-10){
  b <- -1*b
  a - (1/2)*(b)*(B2) - sqfunc(-x,b,s) + sqfunc(B2-x,b,s)
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
  
  fit.lag<-try(nlsLM(y ~ lag(x,a,b,B1,s=1E-10),
                 start = c(B1=mean(x)-(mean(x)-min(x))/2, a=min(y), b=1),data = data,
                 lower = c(B1=-Inf,a=-Inf,b=0.0001),
                 control = nls.control(maxiter=1000, warnOnly=TRUE)),silent=T)
  if(class(fit.lag)=='try-error'){
    fit.lag<-try(nlsLM(y ~ lag(x,a,b,B1,s=1E-10),
                       start = c(B1=10, a=min(y), b=1),data = data,
                       lower = c(B1=-Inf,a=-Inf,b=0.0001),
                       control = nls.control(maxiter=1000, warnOnly=TRUE)),silent=T)
  }
  if(class(fit.lag)=='try-error'){
    if(!grepl(attr(fit.lag,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.lag,"condition"))
    }
    #print('fit.lag failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.lag)))
    
    if(plotQ){
      if(!is.na(fpath)){
        pdf(fpath)
        plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        curve(lag(x,cfs$a,cfs$b,cfs$B1,s=1E-10),min(x),max(x),n = 400,add=T,col='blue')
        curve(lag(x,cfs$a,cfs$b,cfs$B1,s=1E-10),cfs$B1,max(x),n = 400,add=T,col='red')
        dev.off()
      }else{
        plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        curve(lag(x,cfs$a,cfs$b,cfs$B1,s=1E-10),min(x),max(x),n = 400,add=T,col='blue')
        curve(lag(x,cfs$a,cfs$b,cfs$B1,s=1E-10),cfs$B1,max(x),n = 400,add=T,col='red')
      }
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
  
  fit.sat<-try(nlsLM(y ~ sat(x,a,b,B2,s=1E-10),
                 start=c(B2=mean(x)+(max(x)-mean(x))/2,a=a.guess,b=round(max(c(slopes,0.0001)),5)),
                 data = data,
                 lower = c(B2=-Inf,a=-Inf,b=0.0001),
                 control = nls.control(maxiter=1000, warnOnly=TRUE)),silent=T)
  if(class(fit.sat)=='try-error'){
    fit.sat<-try(nlsLM(y ~ sat(x,a,b,B2,s=1E-10),
                       start=c(B2=10,a=a.guess,b=round(max(c(slopes,0.0001)),5)),
                       data = data,
                       lower = c(B2=-Inf,a=-Inf,b=0.0001),
                       control = nls.control(maxiter=1000, warnOnly=TRUE)),silent=T)
  }
  if(class(fit.sat)=='try-error'){
    if(!grepl(attr(fit.sat,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.sat,"condition"))
    }
    #print('fit.sat failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.sat)))
    
    if(plotQ){
      if(!is.na(fpath)){
        pdf(fpath)
        plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        curve(sat(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),max(x),add=T,col='blue')
        curve(sat(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),cfs$B2,add=T,col='red')
        dev.off()
      }else{
        plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        curve(sat(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),max(x),add=T,col='blue')
        curve(sat(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),cfs$B2,add=T,col='red')
      }
    }
  }  
  
  return(fit.sat)
}

#' Extract exponential growth rate assuming exponential death that hits a floor
#' 
#' This function fits a smoothed piecewise linear model to ln(abundance) data, with 
#' the assumption that abundances decrease linearly at first, but then hit a floor
#' and remain constant - consistent with a population that declines to the detection
#' limit of fluoresence.
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
get.gr.flr<-function(x,y,plotQ=F,fpath=NA,id=''){
  
  data<-data.frame(x=x,y=y)
  slopes <- rollapply(data.frame(x=x,y=y), 3, localslope, by.column=F)
  a.guess<-coef(lm(y~x))[[1]]
  
  fit.flr<-try(nlsLM(y ~ flr(x,a,b,B2,s=1E-10),
                 start=c(a=a.guess,b=min(c(-0.1,round(min(slopes),5))),B2=mean(x)),
                 data = data,
                 upper = c(a=Inf,b=-0.00000001,B2=Inf),
                 control = nls.control(maxiter=1000, warnOnly=TRUE)),silent=T)
  
  if(class(fit.flr)=='try-error'){
    if(!grepl(attr(fit.flr,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.flr,"condition"))
    }
    #print('fit.flr failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.flr)))
    
    if(plotQ){
      if(!is.na(fpath)){
        pdf(fpath)
        plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        curve(flr(x,cfs$a,cfs$b,cfs$B2,s=1E-10),cfs$B2,max(x),add=T,col='blue')
        curve(flr(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),cfs$B2,add=T,col='red')
        dev.off()
      }else{
        plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        curve(flr(x,cfs$a,cfs$b,cfs$B2,s=1E-10),cfs$B2,max(x),add=T,col='blue')
        curve(flr(x,cfs$a,cfs$b,cfs$B2,s=1E-10),min(x),cfs$B2,add=T,col='red')
      }
    }
  }
  return(fit.flr)
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
  
  fit.lagsat<-try(nlsLM(y ~ lagsat(x,a,b,B1,B2,s=1E-10),
                    start = c(B1=mean(x)-(mean(x)-min(x))/2,B2=mean(x)+(max(x)-mean(x))/2, a=min(y)+0.1, b=1),
                    data = data,
                    lower = c(B1=-Inf,B2=-Inf,a=-Inf,b=0.0001),
                    control = nls.control(maxiter=1000, warnOnly=TRUE)),silent=T)
  if(class(fit.lagsat)=='try-error'){
    if(!grepl(attr(fit.lagsat,"condition"),pattern='singular gradient matrix')){
      print(attr(fit.lagsat,"condition"))
    }
    #print('fit.lagsat failed after two tries')
  }else{
    cfs<-data.frame(t(coef(fit.lagsat)))
    
    if(plotQ){
      if(!is.na(fpath)){
        pdf(fpath)
        plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        curve(lagsat(x,cfs$a,cfs$b,cfs$B1,cfs$B2,s=1E-10),min(x),max(x),add=T,col='blue')
        curve(lagsat(x,cfs$a,cfs$b,cfs$B1,cfs$B2,s=1E-10),cfs$B1,cfs$B2,add=T,col='red')
        dev.off()
      }else{
        plot(y~x,xlab='Time (days)',ylab='ln(fluorescence)',main=id)
        curve(lagsat(x,cfs$a,cfs$b,cfs$B1,cfs$B2,s=1E-10),min(x),max(x),add=T,col='blue')
        curve(lagsat(x,cfs$a,cfs$b,cfs$B1,cfs$B2,s=1E-10),cfs$B1,cfs$B2,add=T,col='red')
      }
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
#' @param methods Must be a character vector containing one or more of \code{'linear'}, \code{'lag'}, \code{'sat'}, \code{'flr'}, or \code{'lagsat'}
#' @param id Label corresponding to the population/strain/species of interest; used to determine the title and file name of saved plot, if any.
#' @param model.selection control parameter to specify which IC metric to use in model selection; default is AICc, which corrects for small sample sizes and converges asymptotically on AIC.
#' @param internal.r2.cutoff control parameter specifying the R2 criteria that may be applied to drop fits where the number of observations in the exponential portion is equal to 3. The default value of zero permits all fits of 3 obs to be considered.
#' @param zero.time if TRUE, shift time axis so that each time series starts at time = 0
#' 
#' @return A data frame containing the identity of the best model, the content of the best model, the estimated slopes of the increasing linear portion of the regressions (ie, exponential growth rate), the standard errors associated with these slopes, the IC table used to determine the best model, and the full list of all models fit. See vignette for details.
#' 
#' 
#' @examples 
#' 
#' @export
#' @import bbmle
get.growth.rate<-function(x,y,id,plot.best.Q=F,fpath=NA,methods=c('linear','lag','sat','flr','lagsat'),model.selection=c('AICc'),internal.r2.cutoff=0,verbose=FALSE,zero.time=T){
  
  # thin vectors if abundance measure is NA
  x<-x[!is.na(y)]
  y<-y[!is.na(y)]
  
  if(zero.time){
    x <- x-min(x,na.rm = T)
  }
  
  if(sum(methods %in% c('linear','lag','sat','flr','lagsat'))==0){
    print('Error! None of the specified methods matched a currently implemented approach')
  }
  
  if(verbose){print(paste('data set id = ',id))}
  
  # Initialize empty data structures
  modlist<-list(gr=NA,gr.lag=NA,gr.sat=NA,gr.flr=NA,gr.lagsat=NA)
  class(modlist$gr)<-class(modlist$gr.lag)<-class(modlist$gr.sat)<-class(modlist$gr.flr)<-class(modlist$gr.lagsat)<-'try-error'
  gr<-gr.lag<-gr.sat<-gr.flr<-gr.lagsat<-NA
  slope.gr<-slope.gr.lag<-slope.gr.sat<-slope.gr.flr<-slope.gr.lagsat<-NA
  slope.n.gr<-slope.n.gr.lag<-slope.n.gr.sat<-slope.n.gr.flr<-slope.n.gr.lagsat<-NA
  slope.r2.gr<-slope.r2.gr.lag<-slope.r2.gr.sat<-slope.r2.gr.flr<-slope.r2.gr.lagsat<-NA
  se.gr<-se.gr.lag<-se.gr.sat<-se.gr.flr<-se.gr.lagsat<-NA
  
  if(length(unique(x))==2){
    print('Caution: only two unique time points, high risk of over-fitting. Methods other than "linear" are likely to fail')
  }
  
  # if there are more than two unique time points with data:
  if(length(unique(x))>=2){
  
    # Fill data structures for each given method, if requested by user:
    if('linear' %in% methods){
      gr<-get.gr(x,y) 
      slope.gr<-coef(gr)[[2]]
      slope.n.gr<-length(x)
      slope.r2.gr<-get.R2(gr,y)
      se.gr<-sqrt(diag(vcov(gr)))['x']
      modlist$gr<-gr
    }
    if('lag' %in% methods){
      gr.lag<-try(get.gr.lag(x,y))  
      if(prod(class(gr.lag)!='try-error')){
        b1.cutoff <- coef(gr.lag)[1]-0.1  # where does exponential phase begin?
        pds.lag <- predict(gr.lag)[x>=b1.cutoff] # predicted values above this cutoff
        obs.lag <- y[x>=b1.cutoff] # observed values above this cutoff
        
        slope.gr.lag <- unname(coef(gr.lag)['b'])
        slope.n.gr.lag <- length(x[x>=b1.cutoff])  # how many observations above cutoff
        slope.r2.gr.lag <- 1-sum((pds.lag-obs.lag)^2)/sum((obs.lag-mean(obs.lag))^2)  # r2
        se.gr.lag <- sqrt(diag(vcov(gr.lag)))['b']
        
        # if exponential portion is based on fewer than 3 observations, re-classify this fit as
        # resulting in an error. This removes it from consideration as a 'best model', allowing
        # a different model to succeed.
        if(slope.n.gr.lag<3  | (slope.n.gr.lag==3 & slope.r2.gr.lag<internal.r2.cutoff)){
          class(gr.lag)<-'try-error'          
        }
      }else{
        slope.gr.lag<-NA
        slope.n.gr.lag<-NA
        slope.r2.gr.lag<-NA
        se.gr.lag<-NA
      }
      modlist$gr.lag<-gr.lag
    }
    if('sat' %in% methods){
      gr.sat<-try(get.gr.sat(x,y))
      if(prod(class(gr.sat)!='try-error')){
        b2.cutoff <- coef(gr.sat)[1]+0.1  # where does exponential phase end?
        pds.sat <- predict(gr.sat)[x<=b2.cutoff] # predicted values below this cutoff
        obs.sat <- y[x<=b2.cutoff] # observed values below this cutoff
        
        slope.gr.sat <- unname(coef(gr.sat)['b'])
        slope.n.gr.sat <- length(x[x<=b2.cutoff])  # how many observations below cutoff
        slope.r2.gr.sat <- 1-sum((pds.sat-obs.sat)^2)/sum((obs.sat-mean(obs.sat))^2)  # r2
        se.gr.sat <- sqrt(diag(vcov(gr.sat)))['b']
        
        # if exponential portion is based on fewer than 3 observations, re-classify this fit as
        # resulting in an error. This removes it from consideration as a 'best model', allowing
        # a different model to succeed.
        if(slope.n.gr.sat<3 | (slope.n.gr.sat==3 & slope.r2.gr.sat<internal.r2.cutoff)){
          class(gr.sat)<-'try-error'          
        }
      }else{
        slope.gr.sat <- NA
        slope.n.gr.sat <- NA
        slope.r2.gr.sat <- NA
        se.gr.sat <- NA
      }
      modlist$gr.sat<-gr.sat
    }
    
    if('flr' %in% methods){
      gr.flr<-try(get.gr.flr(x,y))
      if(prod(class(gr.flr)!='try-error')){
        b2.cutoff <- coef(gr.flr)['B2']+0.1  # where does exponential phase end?
        pds.flr <- predict(gr.flr)[x<=b2.cutoff] # predicted values below this cutoff
        obs.flr <- y[x<=b2.cutoff] # observed values below this cutoff
        
        slope.gr.flr <- unname(coef(gr.flr)['b'])
        slope.n.gr.flr <- length(x[x<=b2.cutoff])  # how many observations below cutoff
        slope.r2.gr.flr <- 1-sum((pds.flr-obs.flr)^2)/sum((obs.flr-mean(obs.flr))^2)  # r2
        se.gr.flr <- sqrt(diag(vcov(gr.flr)))['b']
        
        # if exponential portion is based on fewer than 3 observations, re-classify this fit as
        # resulting in an error. This removes it from consideration as a 'best model', allowing
        # a different model to succeed.
        if(slope.n.gr.flr<3 | (slope.n.gr.flr==3 & slope.r2.gr.flr<internal.r2.cutoff)){
          class(gr.flr)<-'try-error'          
        }
      }else{
        slope.gr.flr <- NA
        slope.n.gr.flr <- NA
        slope.r2.gr.flr <- NA
        se.gr.flr <- NA
      }
      modlist$gr.flr<-gr.flr
    }
    
    if('lagsat' %in% methods){
      gr.lagsat<-try(get.gr.lagsat(x,y))
      if(prod(class(gr.lagsat)!='try-error')){
        b1.cutoff <- coef(gr.lagsat)[1]-0.1  # where does exponential phase begin?
        b2.cutoff <- coef(gr.lagsat)[2]+0.1  # where does exponential phase end?
        pds.lagsat <- predict(gr.lagsat)[x<=b2.cutoff & x >=b1.cutoff] # predictions btwn cutoffs
        obs.lagsat <- y[x<=b2.cutoff & x >=b1.cutoff] # observed values between cutoffs
        
        
        slope.gr.lagsat <- unname(coef(gr.lagsat)['b'])
        slope.n.gr.lagsat <- length(x[x<=b2.cutoff & x >=b1.cutoff])  # how many obs btwn cutoffs
        slope.r2.gr.lagsat <- 1-sum((pds.lagsat-obs.lagsat)^2)/sum((obs.lagsat-mean(obs.lagsat))^2)  # r2
        se.gr.lagsat <- sqrt(diag(vcov(gr.lagsat)))['b']
        
        # if exponential portion is based on fewer than 3 observations, re-classify this fit as
        # resulting in an error. This removes it from consideration as a 'best model', allowing
        # a different model to succeed.
        if(slope.n.gr.lagsat<3  | (slope.n.gr.lagsat==3 & slope.r2.gr.lagsat<internal.r2.cutoff)){
          class(gr.lagsat)<-'try-error'          
        }

      }else{
        slope.gr.lagsat <- NA
        slope.n.gr.lagsat <- NA
        slope.r2.gr.lagsat <- NA
        se.gr.lagsat <- NA
      }
      modlist$gr.lagsat<-gr.lagsat
    }
  
    # determine which fits occured and were successful
    successful.fits<-sapply(modlist,detect)
  
    if(sum(successful.fits)==0){
      print('Error! All results for requested methods failed!')
      break();
    }
  
    # assemble model names, contents, and slopes, but only for successful fits
    mod.names<-c('gr','gr.lag','gr.sat','gr.flr','gr.lagsat')[successful.fits]
    mod.list<-list(gr,gr.lag,gr.sat,gr.flr,gr.lagsat)[successful.fits]
    slope.ests<-unname(c(slope.gr,slope.gr.lag,slope.gr.sat,slope.gr.flr,slope.gr.lagsat)[successful.fits])
    se.ests<-unname(c(se.gr,se.gr.lag,se.gr.sat,se.gr.flr,se.gr.lagsat)[successful.fits])
    slope.n.vals<-c(slope.n.gr,slope.n.gr.lag,slope.n.gr.sat,slope.n.gr.flr,slope.n.gr.lagsat)[successful.fits]
    slope.r2.vals<-c(slope.r2.gr,slope.r2.gr.lag,slope.r2.gr.sat,slope.r2.gr.flr,slope.r2.gr.lagsat)[successful.fits]
      
    # compare successful models
    switch(model.selection,
           AIC={ictab<-AICtab(mod.list,mnames = mod.names)},
           AICc={ictab<-AICctab(mod.list,mnames = mod.names)},
           BIC={ictab<-BICtab(mod.list,mnames = mod.names)},
           print("Error! Invalid IC method selected in get.nbcurve()"))
  
    best.mod.id<-which(mod.names==attr(ictab,"row.names")[1])
    
    # impose QC based on slope.n and slope.r2 here? or outside of function...
    
    # match model names to outputs:
    names(slope.ests)<-mod.names
    names(se.ests)<-mod.names
    names(slope.n.vals)<-mod.names
    names(slope.r2.vals)<-mod.names
      
    # format output:
    result<-list(best.slope=slope.ests[[best.mod.id]],
                 best.se=se.ests[[best.mod.id]],
                 best.model=as.character(mod.names[[best.mod.id]]),
                 best.model.rsqr=get.R2(mod.list[[best.mod.id]],y),
                 best.model.slope.n=slope.n.vals[[best.mod.id]],
                 best.model.slope.r2=slope.r2.vals[[best.mod.id]],
                 best.model.contents=list(mod.list[[best.mod.id]]),
                 slopes=slope.ests,
                 ses=se.ests,
                 slope.ns=slope.n.vals,
                 slope.rs=slope.r2.vals,
                 ictab=ictab,
                 models=list(gr=gr,gr.lag=gr.lag,gr.sat=gr.sat,gr.flr=gr.flr,gr.lagsat=gr.lagsat))
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
                 gr.flr=get.gr.flr(x,y,plotQ=T,fpath=fpath,id=id[1]),
                 gr.lagsat=get.gr.lagsat(x,y,plotQ=T,fpath=fpath,id=id[1]))
    }
  }else{
    print("Warning: fewer than two unique time points provided")
    
    result<-list(best.slope=NA,
                 best.se=NA,
                 best.model="NA",
                 best.model.rsqr=NA,
                 best.model.slope.n=NA,
                 best.model.slope.r2=NA,
                 best.model.contents=list(NA),
                 slopes=NA,
                 ses=NA,
                 slope.ns=NA,
                 slope.rs=NA,
                 models=list(gr=NA,gr.lag=NA,gr.sat=NA,gr.flr=NA,gr.lagsat=NA))
  }
  
  return(result)
}


# These two functions help to trim data sets to exclude temperature treatments that are well beyond a species' thermal niche, by identifying the minimum and maximum temperatures where positive growth was observed, and determining the next highest (lowest) temperature treatment.
#funky.low<-function(x,y){
#  dat<-data.frame(x,y) 
#  dat2<-dat %>% group_by(x) %>% summarise(max.y=max(y))
#  Ts<-sort(unique(dat2$x))
  
#  minT<-min(dat2$x[dat2$max.y>0])
#  res<-Ts[max(c(1,which(Ts==minT)-1))]
  
#  return(res)
#}

#funky.high<-function(x,y){
#  dat<-data.frame(x,y) 
#  dat2<-dat %>% group_by(x) %>% summarise(max.y=max(y))
#  Ts<-sort(unique(dat2$x))
#  
#  maxT<-max(dat2$x[dat2$max.y>0])
#  res<-Ts[min(c(length(Ts),which(Ts==maxT)+1))]
#  
#  return(res)
#}
