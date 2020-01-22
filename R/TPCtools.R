
#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Conduct fits of thermal reaction norms to various parametric equations

# - currently focus is on Norberg fits
# - expand to include double exponential
# - allow user-selected parametric equations and fitting algorithms

#----------------------------------------------------------------#
#----------------------------------------------------------------#


#' Norberg function
#' 
#' See Norberg et al. 2001, Thomas et al. 2012
#' 
#' Note that with the original formulation, 'Optimum temperature' is the temperature at
#' which growth rate coincides with the exponential portion of this equation. Only if b
#' is exactly 0 will this match the temperature at which a species achieves its highest
#' growth rate. See Thomas et al. 2012 for more on this, and \code{nbcurve()} for an
#' alternative formulation.
#' 
#' @param x Temperature
#' @param copt Competitive optimum temperature
#' @param w Thermal niche width
#' @param a ln(Exponential intercept)
#' @param b Exponential scaling
#' 
#' @return Predicted exponential growth rate at temperature x
#' 
#' @export
nbcurve.legacy<-function(x,copt,w,a,b){
  res<-exp(a)*exp(b*x)*(1-((x-copt)/(w/2))^2)
  res
}

#' Re-parameterized Norberg functions
#' 
#' With this formulation, optimum temperature is explicitly the temperature at which
#' a species achieves its highest growth rate.
#' 
#' @param x Temperature
#' @param topt Optimum temperature
#' @param w Thermal niche width
#' @param a Affects the y-intercept
#' @param b Exponential scaling
#' 
#' @return Predicted exponential growth rate at temperature x
#' 
#' @export
nbcurve<-function(x,topt,w,a,b){
  num <- -2-2*b*topt+2*b*x+sqrt(4+(b*w)^2)
  res<-exp(a+b*x)*(1-(num/(b*w))^2)
  res
}

#' Get optimum temperature from nbcurve.legacy()
#' 
#' Given the competitive optimum temperature and other legacy Norberg parameters, 
#' calculate the true optimum temperature.
#' 
#' @param copt Competitive optimum temperature
#' @param w Thermal niche width
#' @param b Exponential scaling
#' 
#' @return Optimum temperature
#' 
#' @export
get.topt<-function(copt,w,b){
  (1/(2*b))*(-2+2*b*copt+sqrt(4+b^2*w^2))
}

#' Calculate thermal niche limits
#' 
#' Note: this operates on nbcurve parameters (ie, uses topt not copt)
#' 
#' @param topt Optimum temperature
#' @param w Thermal niche width
#' @param b Exponential scaling
#' @param type One of either 'tmin' or 'tmax'; default is 'tmin'
#' 
#' @return Either the estimate of Tmin or Tmax corresponding to the provided Norberg parameters.
#' 
#' @export
get.tlim<-function(topt,w,b,type='tmin'){
  if(!(type %in% c('tmin','tmax'))) print('Error in get.tlim; unknown type requested')
  
  if(type=='tmax'){
    res <- (2+2*b*topt+b*w-sqrt(4+b^2*w^2))/(2*b)
  }
  if(type=='tmin'){
    res <- (2+2*b*topt-b*w-sqrt(4+b^2*w^2))/(2*b)
  }
  return(res)
}


#' Double exponential model
#' 
#' Re-parameterized from Thomas et al. 2016 to depend explicitly on Topt
#'
#' @param temperature Temperature(s) to calculate growth rate at
#' @param topt Optimum temperature (where growth rate is highest)
#' @param b1 Birth rate at 0 Celsius
#' @param b2 Temperature scaling of birth rate, > 0
#' @param d0 Temperature-independent death rate
#' @param d2 Temperature scaling of death rate, > 0 (and d2 > b2)
#'
#' @return Growth rate at one or more temperatures
#' 
#' @export   
decurve<-function(temperature,topt,b1,b2,d0,d2){
  res <- b1*exp(b2*temperature) - (d0 + ((b1*b2)/d2)*exp((b2-d2)*topt)*exp(d2*temperature))
  res
}


#' Double exponential model, phi version
#' 
#' Alternate parameterization for Topt and phi (see mathematica notebook)
#'
#' @param temperature Temperature(s) to calculate growth rate at
#' @param topt Optimum temperature (where growth rate is highest)
#' @param phi Birth rate at 0 Celsius
#' @param b2 Temperature scaling of birth rate, > 0
#' @param d0 Temperature-independent death rate
#' @param d2 Temperature scaling of death rate, > 0 (and d2 > b2)
#'
#' @return Growth rate at one or more temperatures
#' 
#' @export   
decurve2<-function(temperature,topt,phi,b2,d0,d2){
  res <- (phi/(b2*(b2-d2)))*exp(b2*(temperature-topt))-(d0+(phi/(d2*(b2-d2)))*exp(d2*(temperature-topt)))
  res
}


#' Fit Norberg curve to growth rate vs. temperature data
#' 
#' Note: this function fits a re-parameterized version of the original Norberg curve (Norberg et al. 2001, Thomas et al. 2012), altered to depend directly on a parameter that provides the true optimum temperature (ie, the temperature at which growth rate is maximal). See \code{nbcurve()} for more.
#' 
#' @param temperature Temperature
#' @param mu Exponential growth rate
#' @param method Specify which fitting algorithm to use, 'mle2' or 'grid.mle2'
#' @param plotQ Should regression be visualized?
#' @param conf.bandQ Should we calculate a confidence band around the regression? logical.
#' @param fpath If visual requested, and valid file path provided here, plot will be saved as a .pdf file. Default is NA.
#' @param id Character string providing any information ID'ing the specifc curve being fit; used to label plots, if any are requested. Default is NA.
#' @param suppress.grid.mle2.warnings logical; should warnings arising from grid.mle2 invocation be suppressed (TRUE), or displayed (FALSE)? Default is TRUE.
#' @param ... Additional arguments passed to grid.mle2 (e.g., control=list(maxit=2000))
#' 
#' @examples 
#' data("example_TPC_data")
#' sp1 <- example_TPC_data %>% filter(isolate.id=='CH30_4_RI_03' & dilution==1)
#' nbcurve.traits<-get.nbcurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',
#' plotQ=TRUE,conf.bandQ = TRUE,fpath=NA)
#' nbcurve.traits
#' 
#' @export
#' @importFrom bbmle mle2 vcov predict.mle2
#' @import emdbook
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#' @import mleTools
#' @importFrom stats na.omit uniroot AIC
get.nbcurve.tpc<-function(temperature,mu,method='grid.mle2',plotQ=FALSE,conf.bandQ=TRUE,
                          fpath=NA,id=NA,suppress.grid.mle2.warnings=TRUE,...){
  tpc.tmp<-na.omit(data.frame(temperature,mu))
  ntemps<-length(unique(tpc.tmp$temperature))
  
  if(ntemps<=4){
    print("Caution in get.nbcurve.tpc - focal data set has <=4 unique temperatures, risk of overfitting is high!")
  }
  id<-id[1]
  
  if(method=='grid.mle2'){
    
    # set up search of a grid of parameter guesses
    grids<-list(topt=seq(15,35,5),w=seq(10,40,5),a=seq(-0.5,-3,-0.5),b=c(-0.05,0,0.05))
    start<-list(topt=NA,w=NA,a=NA,b=NA,s=log(2))
    
    if(suppress.grid.mle2.warnings){
      fit0<-suppressWarnings(grid.mle2(minuslogl=mu~dnorm(mean=nbcurve(temperature,topt,w,a,b),
                                                          sd=exp(s)),
                                       grids=grids,start=start,data=tpc.tmp,...))
    }else{
      fit0<-grid.mle2(minuslogl=mu~dnorm(mean=nbcurve(temperature,topt,w,a,b),sd=exp(s)),
                                       grids=grids,start=start,data=tpc.tmp,...)
    }
    cfg<-coef(fit0$res.best) # this seemed to be throwing problems b/c of an issue with accessing mle2...?

    # polish best fit model, using formula interface:
    guesses<-as.list(cfg)
    fit<-bbmle::mle2(mu~dnorm(mean=nbcurve(temperature,topt,w,a,b),sd=exp(s)),
              start=guesses,data=tpc.tmp)
  }
  
  if(method=='mle2'){
    print("Caution: this option not fully beta tested... 7/30/18")
    topt.guess <- tpc.tmp$temperature[tpc.tmp$mu==max(tpc.tmp$mu)]
    w.guess <- diff(range(tpc.tmp$temperature))
    a.guess <- -1.11
    b.guess <- 0.05
    fit<-bbmle::mle2(mu~dnorm(mean=nbcurve(temperature,topt,w,a,b),sd=exp(s)),
              start=list(topt=topt.guess,w=w.guess,a=a.guess,b=b.guess,s=log(2)),data=tpc.tmp)
  }
  
  # pull out parameters
  cf<-as.list(coef(fit))
  vcov.mat<-bbmle::vcov(fit)
  
  # additional responses/traits:
  tmin<-get.tlim(cf$topt,cf$w,cf$b,type='tmin')
  tmax<-get.tlim(cf$topt,cf$w,cf$b,type='tmax')
  
  # calculate R2
  rsqr<-get.R2(predict(fit),tpc.tmp$mu)

  # Calculate umax and confidence interval using the delta method (see Bolker book, pg 255)
  pd.umax<-predict(fit,newdata=data.frame(temperature=cf$topt))
  st.umax<-paste("nbcurve(c(",paste(cf$topt,collapse=','),"),topt,w,a,b)",sep='')
  dvs0.umax<-suppressWarnings(deltavar2(fun=parse(text=st.umax),meanval=cf,Sigma=vcov.mat))
  
  # simple Fisher confidence intervals:
  ciF<-ci.FI(fit)
  
  # save output:
  vec<-as.list(c(cf,rsqr=rsqr,tmin=tmin,tmax=tmax))
  vec$umax<-pd.umax
  vec$umax.ci<-1.96*sqrt(dvs0.umax)
  vec$ciF<-ciF
  vec$cf<-cf
  vec$vcov<-vcov.mat
  vec$nobs<-nrow(tpc.tmp)
  vec$ntemps<-ntemps
  vec$logLik<-logLik(fit)
  vec$aic<-AIC(fit)
  vec$data<-tpc.tmp
  
  # Plot results:
  if(plotQ){
    
    # use function to generate a single curve plot:
    p1<-plot.nbcurve(vec,xlim = c(min(tpc.tmp$temperature),
                                  max(tpc.tmp$temperature)),
                         plot.ci = conf.bandQ,plot.obs = TRUE)
    p1<-p1+scale_x_continuous('Temperature (C)')+
      scale_y_continuous('Growth rate (1/day)')+
      theme(panel.grid = element_blank())
    
    if(!is.na(id)){
      p1<-p1+ggtitle(id)
    }

    if(!is.na(fpath)){
      if(is.na(id)){
        time<-Sys.time()
        time<-gsub(x = time,pattern = ":",replacement = "_")
        full.path<-paste(fpath,"TPC_fit_",time,".pdf",sep='')        
      }else{
        full.path<-paste(fpath,"TPC_fit_",id,".pdf",sep='') 
      }
      ggsave(device = pdf(),filename = full.path,p1,width = 5,height = 4)
      dev.off()
    }else{
      print(p1)
    }
  }

  # Finished, return relevant stats.  
  return(vec)
}




#' Fit Double Exponential curve to growth rate vs. temperature data
#' 
#' @param temperature Temperature
#' @param mu Exponential growth rate
#' @param method Specify which fitting algorithm to use, 'mle2' or 'grid.mle2'
#' @param start.method Specify method for generating starting grid for 'grid.mle2' option
#' @param plotQ Should regression be visualized?
#' @param conf.bandQ Should we calculate a confidence band around the regression? logical.
#' @param fpath If visual requested, and valid file path provided here, plot will be saved as a .pdf file. Default is NA.
#' @param id Character string providing any information ID'ing the specifc curve being fit; used to label plots, if any are requested. Default is NA.
#' @param suppress.grid.mle2.warnings logical; should warnings arising from grid.mle2 invocation be suppressed (TRUE), or displayed (FALSE)? Default is TRUE.
#' @param ... Additional arguments passed to grid.mle2 (e.g., control=list(maxit=2000))
#' 
#' @export
#' @importFrom bbmle mle2 vcov predict.mle2
#' @import dplyr
#' @import emdbook
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#' @import mleTools
#' @import mgcv
#' @importFrom stats na.omit uniroot AIC
get.decurve.tpc<-function(temperature,mu,method='grid.mle2',start.method='general.grid',
                          plotQ=FALSE,conf.bandQ=TRUE,fpath=NA,id=NA,suppress.grid.mle2.warnings=TRUE,...){
  tpc.tmp<-na.omit(data.frame(temperature,mu))
  ntemps<-length(unique(tpc.tmp$temperature))
  
  # is this still necessary?
  id<-id[1]
  
  if(ntemps<=5){
    print("Caution in get.decurve.tpc - focal data set has <=5 unique temperatures, risk of overfitting is high!")
  }
  
  if(method=='grid.mle2'){
    
    if(start.method=='general.grid'){
      tpc.tmp.tb<-tpc.tmp %>% group_by(temperature) %>% summarise(mu=mean(mu))
      topt.guess<-mean(tpc.tmp.tb$temperature[tpc.tmp.tb$mu==max(tpc.tmp.tb$mu)])
      
      # set up search of a grid of parameter guesses
      grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
                  d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
      start<-list(topt=topt.guess,b1=NA,b2=NA,d0=NA,d2=NA,s=log(2))
      
      if(suppress.grid.mle2.warnings){
        fit0<-suppressWarnings(grid.mle2(minuslogl=mu~dnorm(mean=decurve(temperature,topt,exp(b1),exp(b2),exp(d0),exp(d2)),sd=exp(s)),grids=grids,start=start,data=tpc.tmp,...))
      }else{
        fit0<-grid.mle2(minuslogl=mu~dnorm(mean=decurve(temperature,topt,exp(b1),exp(b2),exp(d0),exp(d2)),sd=exp(s)),grids=grids,start=start,data=tpc.tmp,...)
      }
      cfg<-as.list(coef(fit0$res.best))
      
      # extract parameters for polished fit
      guesses<-as.list(cfg)
      guesses<-list(topt=cfg$topt,b1=exp(cfg$b1),b2=exp(cfg$b2),d0=exp(cfg$d0),d2=exp(cfg$d2),s=exp(cfg$s))
    }
    
    if(start.method=='smart.grid'){
      
      # guess topt
      tpc.tmp.tb<-tpc.tmp %>% group_by(temperature) %>% summarise(mu=mean(mu))
      topt.guess<-mean(tpc.tmp.tb$temperature[tpc.tmp.tb$mu==max(tpc.tmp.tb$mu)])
      
      # guess phi
      gam.tmp<-gam(mu~s(temperature,k=4),data=tpc.tmp)
      h<-0.1
      pd<-predict(gam.tmp,newdata=data.frame(temperature=c(topt.guess-h,topt.guess,topt.guess+h)))
      phi.guess<-fd2.central(pd,h)
      
      # could introduce a d0.guess, if negative growth rates are present at 
      # multiple low temperatures
      
      # set up search of a grid of parameter guesses
      #   - could consider the merits of moving expand.grid outside of grid.mle2, to allow thinning based on know relationships (ie, d2 > b2)?
      #   - or, add new capacity 'rules' that are passed to grid.mle2, and applied internally to the expand.grid list... hmm. Pass these as 'formula' object
      grids<-list(phi=c(0.8*phi.guess,phi.guess,1.2*phi.guess),b2=log(seq(0.1,0.5,0.2)),
                  d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
      start<-list(topt=topt.guess,phi=NA,b2=NA,d0=NA,d2=NA,s=log(2))
      
      # execute fit
      if(suppress.grid.mle2.warnings){
        fit0<-suppressWarnings(grid.mle2(minuslogl=mu~dnorm(mean=decurve2(temperature,topt,phi,exp(b2),exp(d0),exp(d2)),sd=exp(s)),grids=grids,start=start,data=tpc.tmp))
      }else{
        fit0<-grid.mle2(minuslogl=mu~dnorm(mean=decurve2(temperature,topt,phi,exp(b2),exp(d0),exp(d2)),sd=exp(s)),grids=grids,start=start,data=tpc.tmp)
      }
      cfg<-as.list(coef(fit0$res.best))

      # extract parameters for polish fit, converting back to decurve parameters:
      b1.g<-cfg$phi/(exp(cfg$b2)*(exp(cfg$b2)-exp(cfg$d2))*exp(exp(cfg$b2)*cfg$topt))
      guesses<-list(topt=cfg$topt,b1=b1.g,b2=exp(cfg$b2),d0=exp(cfg$d0),d2=exp(cfg$d2),s=exp(cfg$s))
    }
    
    # polish best fit model using formula interface
    fit<-bbmle::mle2(mu~dnorm(mean=decurve(temperature,topt,b1,b2,d0,d2),sd=s),
              start=guesses,data=tpc.tmp,control=list(maxit=0))
  }

  if(method=='mle2'){
    print("Caution: this option not fully beta tested... 7/30/18")
    topt.guess <- tpc.tmp$temperature[tpc.tmp$mu==max(tpc.tmp$mu)]
    b1.guess <- 0.09
    b2.guess <- 0.15
    d0.guess <- 0.03
    d2.guess <- 0.18
    fit<-bbmle::mle2(mu~dnorm(mean=decurve(temperature,topt,b1,b2,d0,d2),sd=exp(s)),
              start=list(topt=topt.guess,b1=b1.guess,b2=b2.guess,d0=d0.guess,d2=d2.guess,s=log(2)),data=tpc.tmp)
  }
  
  # pull out parameters
  cf<-as.list(coef(fit))
  
  # save vcov matrix
  vcov.mat<-bbmle::vcov(fit)
  
  # additional responses/traits:
  objective<-function(x){
    decurve(x,cf$topt,cf$b1,cf$b2,cf$d0,cf$d2)
  }
  #tmax<-uniroot(f = objective,interval = c(cf$topt,1.5*(cf$topt-log(cf$d2/cf$b2)/(cf$b2-cf$d2))))$root
  #tmin<-uniroot(f = objective,interval = c(1.5*(log(cf$d0/cf$b1)/cf$b2),cf$topt))$root  
  rt1<-try(uniroot(f = objective,interval = c(cf$topt,100))$root)
  if(inherits(rt1,'try-error')){
    tmax<-NA
  }else{
    tmax<-rt1
  }
  
  rt2<-try(uniroot(f = objective,interval = c(-400,cf$topt))$root)
  if(inherits(rt2,'try-error')){
    tmin<-NA
  }else{
    tmin<-rt2
  }

  # calculate R2
  rsqr<-get.R2(predict(fit),tpc.tmp$mu)
  
  # Calculate umax and confidence interval using the delta method (see Bolker book, pg 255)
  pd.umax<-predict(fit,newdata=data.frame(temperature=cf$topt))
  st.umax<-paste("decurve(c(",paste(cf$topt,collapse=','),"),topt,b1,b2,d0,d2)",sep='')
  dvs0.umax<-suppressWarnings(deltavar2(fun=parse(text=st.umax),meanval=cf,Sigma=vcov.mat))
  
  # simple Fisher confidence intervals:
  ciF<-ci.FI(fit)
  
  # save output:
  vec<-as.list(c(cf,rsqr=rsqr,tmin=tmin,tmax=tmax))
  vec$umax<-pd.umax
  vec$umax.ci<-1.96*sqrt(dvs0.umax)
  vec$ciF<-ciF
  vec$cf<-cf
  vec$vcov<-vcov.mat
  vec$n<-nrow(tpc.tmp)
  vec$ntemps<-ntemps
  vec$logLik<-logLik(fit)
  vec$aic<-AIC(fit)
  vec$data<-tpc.tmp
  
  # Plot results:
  if(plotQ){
    
    # use function to generate a single curve plot:
    p1<-plot.decurve(vec,xlim = c(min(tpc.tmp$temperature),
                                  max(tpc.tmp$temperature)),
                     plot.ci = conf.bandQ,plot.obs = TRUE)
    p1<-p1+scale_x_continuous('Temperature (C)')+
      scale_y_continuous('Growth rate (1/day)')+
      theme(panel.grid = element_blank())
    
    if(!is.na(id)){
      p1<-p1+ggtitle(id)
    }
    
    if(!is.na(fpath)){
      if(is.na(id)){
        time<-Sys.time()
        time<-gsub(x = time,pattern = ":",replacement = "_")
        full.path<-paste(fpath,"TPC_fit_",time,".pdf",sep='')        
      }else{
        full.path<-paste(fpath,"TPC_fit_",id,".pdf",sep='') 
      }
      ggsave(device = pdf(),filename = full.path,p1,width = 5,height = 4)
      dev.off()
    }else{
      print(p1)
    }
  }
  
  # Finished, return relevant stats.  
  return(vec)
}





#' Calculate R2
#' 
#' 1-sum((obs-pds)^2)/sum((obs-mean(obs))^2)
#' 
#' @param pds Predicted y values (could arise from predict(model) for an mle2 model)
#' @param obs Observed y values
#' 
#' @return R2 value
#' 
#' @export
get.R2<-function(pds,obs){
  #pds<-predict(model)
  rsqr<-1-sum((obs-pds)^2)/sum((obs-mean(obs))^2)
  rsqr
}

#' Delta method function (modified from emdbook)
#' 
#' See `?deltavar()` in package emdbook. Only change was to replace `expr <- as.expression(substitute(fun))` with `expr<-fun`. This was necessary b/c of problems that arose when
#' invoking `deltavar()` inside of another function. Specifically, I wanted to call 
#' `deltavar()` each time an outer function was executed, but apply deltavar to a unique
#' set of x-values upon each execution of the outer function. For reasons I don't 
#' entirely understand, this range of x-values defined within the outer function's 
#' environment was not being adequately passed to the internal environment of deltavar.
#' There's probably a classier/more resilient way of fixing this, but in the short 
#' term, I patched it over by creating `deltavar2()`, and using parse/a custom string in
#' `get.nbcurve.tpc()` when `deltavar2()` is called.
#' 
#' @param fun Function to calculate the variance of, given parameter estimates in meanvals
#' @param vars 'list of variable names: needed if params does not have names, or if some of the values specified in params should be treated as constant'
#' @param meanval 'possibly named vector of mean values of parameters'
#' @param Sigma 'numeric vector of variances or variance-covariance matrix'
#' @param verbose 'print details?'
#' 
#' @export
#' @importFrom stats D
deltavar2<-function (fun, meanval = NULL, vars, Sigma, verbose = FALSE) 
{
  #expr <- as.expression(substitute(fun))    
  expr<-fun  # Changed by CTK
  nvals <- length(eval(expr, envir = as.list(meanval)))
  vecexp <- nvals > 1
  if (missing(vars)) {
    if (missing(meanval) || is.null(names(meanval))) 
      stop("must specify either variable names or named values for means")
    vars <- names(meanval)
  }
  derivs <- try(lapply(vars, D, expr = expr), silent = TRUE)
  symbderivs <- TRUE
  if (inherits(derivs, "try-error")) {
      symbderivs <- FALSE
      warning("some symbols not in derivative table, using numeric derivatives")
      nderivs <- with(as.list(meanval), numericDeriv(expr[[1]], 
                                                     theta = vars))
      nderivs <- attr(nderivs, "gradient")
  }
  else {
    nderivs <- sapply(derivs, eval, envir = as.list(meanval))
  }
  if (verbose) {
    if (symbderivs) {
      cat("symbolic derivs:\n")
      print(derivs)
    }
    cat("value of derivs:\n")
    print(nderivs)
  }
  if (!is.matrix(Sigma) && length(Sigma) > 1) 
    Sigma <- diag(Sigma)
  if (vecexp && is.list(nderivs)) 
    nderivs <- do.call("cbind", nderivs)
  if (is.matrix(nderivs)) {
    r <- apply(nderivs, 1, function(z) c(z %*% Sigma %*% 
                                           matrix(z)))
  }
  else r <- c(nderivs %*% Sigma %*% matrix(nderivs))
  r
}



#' Approximate 2nd derivative (finite difference)
#' 
#' @param fx A vector of 3 function values, f(x-h), f(x), f(x+h)
#' @param h The step size used in the finite difference calculation
#' 
#' @export
fd2.central<-function(fx,h){
  (fx[3]-2*fx[2]+fx[1])/(h^2)
}


#' Predict values from nbcurve fit
#' 
#' @param fit.info The result of a single TPC curve fit from get.nbcurve.tpc
#' @param newdata A new data frame, containing a sequence of `temperature` values at which model predictions should be made; defaults to (-2,40)
#' @param se.fit logical; should standard error values be returned?
#' 
#' @export
predict.nbcurve<-function(fit.info,newdata=data.frame(temperature=seq(-2,40,0.1)),se.fit=FALSE){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(fit.info)==1){
    fit.info<-fit.info[[1]]
  }
  
  # generate predictions across a range of temperatures
  mu<-nbcurve(newdata$temperature, topt = fit.info$topt, w = fit.info$w, a = fit.info$a, b = fit.info$b)
  newdata$mu<-mu
  
  if(se.fit){
    st<-paste("nbcurve(c(",paste(newdata$temperature,collapse=','),"),topt,w,a,b)",sep='')
    dvs0<-suppressWarnings(deltavar2(fun=parse(text=st),meanval=fit.info$cf,Sigma=fit.info$vcov))
    newdata$se.fit<-sqrt(dvs0)
    
    # Better approach? Pass correct local environment to deltavar... BUSTED
    #dvs0<-deltavar3(fun=nbcurve(xs,topt,w,a,b),meanval=cf,Sigma = vcov(fit),cenv=current.env)
    #dvs0<-suppressWarnings(deltavar(fun=nbcurve(xs,topt,w,a,b),meanval=cf,Sigma=vcov(fit)))
  }
  
  return(newdata)
}



#' Plotting function for single Norberg curve
#' 
#' @param fit.info The result of a single TPC curve fit from get.nbcurve.tpc
#' @param plot.ci logical, should the resulting plot include 95\% confidence bands
#' @param plot.obs logical, should resulting plot include raw data
#' @param xlim x-axis range (temperature)
#' @param ylim y-axis range (adjusts internally to -0.2 to slightly above umax+CI)
#' 
#' @export
plot.nbcurve<-function(fit.info,plot.ci=TRUE,plot.obs=TRUE,xlim=c(-2,40),ylim=c(-0.2,5)){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(fit.info)==1){
    fit.info<-fit.info[[1]]
  }
  
  # adjust plotting window to specific curve?
  ylim<-c(ylim[1],1.1*(fit.info$umax+fit.info$umax.ci))
  
  # generate predictions along a range of temperatures
  if(plot.ci){
    preds<-predict.nbcurve(fit.info,se.fit = TRUE)    
  }else{
    preds<-predict.nbcurve(fit.info)    
  }
  
  # generate basic plot
  cplot<-ggplot(preds,aes(x=.data$temperature,y=.data$mu))+
    geom_hline(yintercept = 0)+
    geom_line()+
    coord_cartesian(xlim=xlim,ylim=ylim)+
    theme_bw()
  
  # add confidence bands?
  if(plot.ci){
    cplot<-cplot+geom_ribbon(aes(ymin=.data$mu-1.96*.data$se.fit,ymax=.data$mu+1.96*.data$se.fit),alpha=0.2)
  }
  
  # add observations?
  if(plot.obs){
    cplot<-cplot+geom_point(data=fit.info$data)
  }
  
  return(cplot) 
}


#' Predict values from decurve fit
#' 
#' @param fit.info The result of a single TPC curve fit from get.decurve.tpc
#' @param newdata A new data frame, containing a sequence of `temperature` values at which model predictions should be made; defaults to (-2,40)
#' @param se.fit logical; should standard error values be returned?
#' 
#' @export
predict.decurve<-function(fit.info,newdata=data.frame(temperature=seq(-2,40,0.1)),se.fit=FALSE){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(fit.info)==1){
    fit.info<-fit.info[[1]]
  }
  
  # generate predictions across a range of temperatures
  mu<-decurve(newdata$temperature, topt = fit.info$topt, b1 = fit.info$b1, b2 = fit.info$b2, d0 = fit.info$d0, d2 = fit.info$d2)
  newdata$mu<-mu
  
  if(se.fit){
    st<-paste("decurve(c(",paste(newdata$temperature,collapse=','),"),topt,b1,b2,d0,d2)",sep='')
    dvs0<-suppressWarnings(deltavar2(fun=parse(text=st),meanval=fit.info$cf,Sigma=fit.info$vcov))
    newdata$se.fit<-sqrt(dvs0)
  }
  
  return(newdata)
}


#' Plotting function for single Double Exponential curve
#' 
#' @param fit.info The result of a single TPC curve fit from get.decurve.tpc
#' @param plot.ci logical, should the resulting plot include 95\% confidence bands
#' @param plot.obs logical, should resulting plot include raw data
#' @param xlim x-axis range (temperature)
#' @param ylim y-axis range (adjusts internally to -0.2 to slightly above umax+CI)
#' 
#' @export
#' @import ggplot2
plot.decurve<-function(fit.info,plot.ci=TRUE,plot.obs=TRUE,xlim=c(-2,40),ylim=c(-0.2,5)){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(fit.info)==1){
    fit.info<-fit.info[[1]]
  }
  
  # adjust plotting window to specific curve?
  ylim<-c(ylim[1],1.1*(fit.info$umax+fit.info$umax.ci))
  
  # generate predictions along a range of temperatures
  if(plot.ci){
    preds<-predict.decurve(fit.info,se.fit = TRUE)    
  }else{
    preds<-predict.decurve(fit.info)    
  }
  
  # generate basic plot
  cplot<-ggplot(preds,aes(x=.data$temperature,y=.data$mu))+
    geom_hline(yintercept = 0)+
    geom_line()+
    coord_cartesian(xlim=xlim,ylim=ylim)+
    theme_bw()
  
  # add confidence bands?
  if(plot.ci){
    cplot<-cplot+geom_ribbon(aes(ymin=.data$mu-1.96*.data$se.fit,ymax=.data$mu+1.96*.data$se.fit),alpha=0.2)
  }
  
  # add observations?
  if(plot.obs){
    cplot<-cplot+geom_point(data=fit.info$data)
  }
  
  return(cplot) 
}
