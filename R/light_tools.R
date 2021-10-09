

#' Eiler-Peeters curve equation
#' 
#' See Eilers & Peeters 1988, Ecological Modeling. Function allows for linearly increasing growth rates at low light levels, 
#' which saturate at maximum growth rate of umax at the optimum light level Lopt, and then decline as light levels increase
#' due to photoinhibition. 
#' 
#' @param x light level, typically umol photons m^-2 s^-1
#' @param a slope of initial increase in growth rate at low light levels, aka affinity
#' @param umax maximum growth rate
#' @param Lopt optimal light level where maximum growth occurs
#' @param log.a logical, is a specified on natural log scale
#' @param log.umax logical, is umax specified on natural log scale
#' @param log.Lopt logical, is Lopt specified on natural log scale
#' 
#' @return Predicted exponential growth rate at light level x
#' 
#' @export
ep_curve<-function(x,a,umax,Lopt,log.a=FALSE,log.umax=FALSE,log.Lopt=FALSE){
  if(log.a){
    a<-exp(a)
  }
  if(log.umax){
    umax<-exp(umax)
  }
  if(log.Lopt){
    Lopt<-exp(Lopt)
  }
  
  res <- umax*(x/((umax*x^2)/(a*Lopt^2) + (1 - 2*(umaxL/(a*Lopt)))*x + umax/a))
  
  res
}


#' Light performance curve (LPC) class creation
#' 
#' @export
new_lpc<-function(){
  lpcObj<-list()
  class(lpcObj)<-"lpc"
  
  # npc type
  lpcObj$type<-character()
  
  # estimated parameters from fit
  lpcObj$model_cf<-list()
  lpcObj$display_cf<-list()
  lpcObj$display_ci<-matrix()
  lpcObj$vcov<-matrix()
  
  # fit diagnostics
  lpcObj$rsqr<-double()
  lpcObj$nlight<-integer()
  lpcObj$nobs<-integer()
  lpcObj$logLik<-double()
  lpcObj$aic<-double()
  
  # data
  lpcObj$data<-data.frame()
  
  lpcObj
}


#' Print method for LPC objects
#' 
#' @param x Object of class lpc
#' @param \dots Additional arguments (not used)
#' 
#' @export
print.lpc<-function(x,...){
  cat("Light performance Curve fit\n")
  cat("Model type: ",x$type,"\n")
  cat("\nEstimated parameters:\n")
  tmp<-cbind(x$display_cf,x$display_ci)
  colnames(tmp)<-c('parameter',colnames(x$display_ci))
  print(tmp)
  # figure out how to make this look tidier? Maybe model after:
  # Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - 
  #digmin)), digits = digits)
  
  cat("\nFit diagnostics:\n")
  cat("R2 = ",round(x$rsqr,4))
  cat(", logLik = ",x$logLik,"(df = ",attr(x$logLik,'df'),"), nobs = ",x$nobs)
  invisible(x)
}


#' Predict values from an LPC fit
#' 
#' @param object The result of a single LPC curve fit, an object of class lpc
#' @param newdata A new data frame, containing a sequence of `light` values at which model predictions should be made; defaults to (0,400)
#' @param se.fit logical; should standard error values be returned?
#' @param \dots Additional arguments (not used)
#' 
#' @export
#' @method predict lpc
predict.lpc<-function(object,newdata,se.fit=FALSE,...){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(object)==1){
    object<-object[[1]]
  }
  
  # instead of providing a default newdata, instead return values from observed 
  # temperatures as default if newdata is missing
  # should also check format of newdata is correct...
  if (missing(newdata) || is.null(newdata)) {
    newdata<-object$data
  }
  
  # generate predictions across a range of temperatures
  switch(object$type,
         ep={mu<-ep_curve(newdata$light, a=object$model_cf$a, umax=object$model_cf$umax, Lopt=object$model_cf$Lopt,
                                log.a=TRUE, log.umax=TRUE, log.Lopt = TRUE)},
         stop(print("unrecognized lpc model type in predict.lpc!")))
  newdata$mu<-mu
  
  if(se.fit){
    insert<-paste(newdata$light,collapse=',')
    
    switch(object$type,
           ep={st<-paste("ep_curve(c(",insert,"),a,umax,Lopt,log.a=TRUE,log.umax=TRUE,log.Lopt=TRUE)",sep='')},
           stop(print("unrecognized lpc model type in predict.lpc se calcs!")))
    dvs0<-suppressWarnings(deltavar2(fun=parse(text=st),meanval=object$model_cf,Sigma=object$vcov))
    newdata$se.fit<-sqrt(dvs0)
  }
  
  return(newdata)
}


#' Plotting function for single light performance curve
#' 
#' @param x The result of a single LPC curve fit from get.ep or similar
#' @param plot_ci logical, should the resulting plot include 95\% confidence bands
#' @param plot_obs logical, should resulting plot include raw data
#' @param xlim x-axis range (nutrients)
#' @param ylim y-axis range (adjusts internally from -0.2 to slightly above umax+CI)
#' @param main Character string providing plot title (usually id info for the plotted curve)
#' @param fpath If visual requested, and valid file path provided, save plot as a .pdf file
#' @param \dots Additional arguments (not used)
#' 
#' @export
#' @import ggplot2
plot.lpc<-function(x,plot_ci=TRUE,plot_obs=TRUE,xlim=NULL,ylim=NULL,main=NA,fpath=NA,...){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(x)==1){
    x<-x[[1]]
  }
  
  # adjust plotting window to specific curve?
  if(is.null(ylim)) ylim <- c(1.3*min(c(min(x$data$mu,na.rm=T)[1],0)),1.3*max(x$data$mu,na.rm=T)[1])
  if(is.null(xlim)) xlim <- c(0,1.2*max(x$data$light,na.rm=T)[1])
  
  # generate predictions along a range of nutrients
  preds<-predict(x,newdata=data.frame(light=seq(xlim[1],xlim[2],length.out = 200)),se.fit = plot_ci)    
  
  # generate basic plot
  cplot<-ggplot(preds,aes(x=.data$light,y=.data$mu))+
    geom_hline(yintercept = 0)+
    geom_line()+
    scale_x_continuous('Light')+
    scale_y_continuous('Growth rate')+
    coord_cartesian(xlim=xlim,ylim=ylim)+
    theme_bw()+
    ggtitle(main)
  
  # add confidence bands?
  if(plot_ci){
    cplot<-cplot+geom_ribbon(aes(ymin=.data$mu-1.96*.data$se.fit,ymax=.data$mu+1.96*.data$se.fit),alpha=0.2)
  }
  
  # add observations?
  if(plot_obs){
    cplot<-cplot+geom_point(data=x$data)
  }
  
  # if function received a file path
  if(!is.na(fpath)){
    if(is.na(main)){ # and a plot ID
      time<-Sys.time()
      time<-gsub(x = time,pattern = ":",replacement = "_")
      full.path<-paste(fpath,"LPC_fit_",time,".pdf",sep='')        
    }else{ # or not
      full.path<-paste(fpath,"LPC_fit_",main,".pdf",sep='') 
    }
    
    # save plot at full.path with given file name.
    ggsave(device = grDevices::pdf(),filename = full.path,cplot,width = 5,height = 4)
    grDevices::dev.off()
    invisible(cplot)
  }else{
    return(cplot) 
  }
}


#' Fit Eiler Peeters curve to growth rate vs. light data
#' 
#' Note: this function currently does not use grid.mle2, as the regressions are usually pretty stable; add this 
#' feature later.
#' 
#' @param light Light level
#' @param mu Exponential growth rate
#' @param method Specify which fitting algorithm to use, 'mle2' or 'grid.mle2'
#' @param ... Additional arguments passed to grid.mle2 (e.g., control=list(maxit=2000))
#' 
#' @examples 
#' 
#' @export
#' @import emdbook
#' @import ggplot2
get.ep<-function(light,mu,method='mle2',...){
  ep.tmp<-stats::na.omit(data.frame(light,mu))
  nlight<-length(unique(ep.tmp$light))
  
  if(nlight<=3){
    print("Caution in get.ep - focal data set has <=3 unique light levels, risk of overfitting is high!")
  }
  
  if(method=='grid.mle2'){
    print("Error: this option not yet implemented in get.ep ... 10/7/21")
  }
  
  if(method=='mle2'){
    #a.guess <- max(monod.tmp$nutrients)[1]/3
    umax.guess <- max(c(max(ep.tmp$mu)[1],0.01))
    Lopt.guess <- max(c(ep.tmp$light[ep.tmp$mu==max(ep.tmp$mu)],10))
    fit<-bbmle::mle2(mu~dnorm(mean=ep_curve(light,a,umax,Lopt,log.a=TRUE,log.umax=TRUE,log.Lopt = TRUE),sd=exp(s)),
                     start=list(a=log(a.guess),umax=log(umax.guess),Lopt=log(Lopt.guess),s=log(2)),data=ep.tmp)
    
    # reframe result on non-logged scale
    model_cf<-as.list(coef(fit))
    tmp.cfs<-model_cf
    tmp.cfs$a<-exp(tmp.cfs$a)
    tmp.cfs$umax<-exp(tmp.cfs$umax)
    tmp.cfs$Lopt<-exp(tmp.cfs$Lopt)
    fit0<-bbmle::mle2(mu~dnorm(mean=ep_curve(light,a,umax,Lopt,log.k=FALSE,log.umax=FALSE),sd=exp(s)),
                      start=tmp.cfs,control=list(maxit=0),
                      data=ep.tmp)
  }
  
  # pull out parameters on regular scale
  display_cf<-as.list(coef(fit0))
  
  # simple Fisher confidence intervals on regular scale:
  ci<-mleTools::ci.FI(fit0)
  
  # vcov on log scale
  vcov_mat<-vcov(fit)
  
  # calculate R2
  rsqr<-get.R2(predict(fit),ep.tmp$mu)
  
  # create empty object of class 'npc'
  vec<-new_lpc()
  
  # populate npc object
  vec$type<-'ep'
  vec$model_cf<-model_cf
  vec$display_cf<-display_cf
  vec$display_ci<-ci
  vec$vcov<-vcov_mat
  vec$nobs<-nrow(ep.tmp)
  vec$nlight<-nlight
  vec$rsqr<-rsqr
  vec$logLik<-logLik(fit)
  vec$aic<-stats::AIC(fit)
  vec$data<-ep.tmp
  
  # Finished, return relevant stats.  
  return(vec)
}
