

#' Monod curve equation
#' 
#' Details here.
#' 
#' @param x nutrients
#' @param umax maximum growth rate
#' @param k half saturation constant
#' @param log.k logical, is k specified on natural log scale
#' @param log.umax logical, is umax specified on natural log scale
#' 
#' @return Predicted exponential growth rate at nutrient concentration x
#' 
#' @export
monod_curve<-function(x,umax,k,log.k=FALSE,log.umax=FALSE){
  if(log.k){
    res<-x/(x+exp(k))
  }else{
    res<-x/(x+k)
  }
  if(log.umax){
    res<-exp(umax)*res
  }else{
    res<-umax*res
  }
  res
}

#' Monod curve with loss term equation
#' 
#' Details here.
#' 
#' @param x nutrients
#' @param umax maximum growth rate
#' @param k half saturation constant
#' @param z ln(intercept)
#' @param log.pars logical, are parameters specified on natural log scale
#' 
#' @return Predicted exponential growth rate at nutrient concentration x
#' 
#' @export
monod_curve_type2<-function(x,umax,k,z,log.pars=FALSE){
  if(log.pars){
    res<-exp(umax)*x/(x+exp(k))-exp(z)
  }else{
    res<-umax*x/(x+k)-z
  }
  res
}


#' Nutrient performance curve (NPC) class creation
#' 
#' @export
new_npc<-function(){
  npcObj<-list()
  class(npcObj)<-"npc"
  
  # npc type
  npcObj$type<-character()
  
  # estimated parameters from fit
  npcObj$model_cf<-list()
  npcObj$display_cf<-list()
  npcObj$display_ci<-matrix()
  npcObj$vcov<-matrix()
  
  # fit diagnostics
  npcObj$rsqr<-double()
  npcObj$nnutr<-integer()
  npcObj$nobs<-integer()
  npcObj$logLik<-double()
  npcObj$aic<-double()
  
  # data
  npcObj$data<-data.frame()
  
  npcObj
}


#' Print method for NPC objects
#' 
#' @param x Object of class npc
#' @param \dots Additional arguments (not used)
#' 
#' @export
print.npc<-function(x,...){
  cat("Nutrient performance Curve fit\n")
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


#' Predict values from an NPC fit
#' 
#' @param object The result of a single NPC curve fit from get.monod, an object of class npc
#' @param newdata A new data frame, containing a sequence of `nutrient` values at which model predictions should be made; defaults to (0,100)
#' @param se.fit logical; should standard error values be returned?
#' @param \dots Additional arguments (not used)
#' 
#' @export
#' @method predict npc
predict.npc<-function(object,newdata,se.fit=FALSE,...){
  
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
         monod={mu<-monod_curve(newdata$nutrients,umax=object$model_cf$umax,k=object$model_cf$k,
                                log.k=TRUE,log.umax=TRUE)},
         monod2={mu<-monod_curve_type2(newdata$nutrients,umax=object$model_cf$umax,k=object$model_cf$k,
                                       z=object$model_cf$z,log.pars=TRUE)},
         stop(print("unrecognized npc model type in predict.npc!")))
  newdata$mu<-mu
  
  if(se.fit){
    insert<-paste(newdata$nutrients,collapse=',')
    
    switch(object$type,
           monod={st<-paste("monod_curve(c(",insert,"),umax,k,log.k=TRUE,log.umax=TRUE)",sep='')},
           monod2={st<-paste("monod_curve_type2(c(",insert,"),umax,k,z,log.pars=TRUE)",sep='')},
           stop(print("unrecognized npc model type in predict.npc se calcs!")))
    dvs0<-suppressWarnings(deltavar2(fun=parse(text=st),meanval=object$model_cf,Sigma=object$vcov))
    newdata$se.fit<-sqrt(dvs0)
  }
  
  return(newdata)
}


#' Plotting function for single nutrient performance curve
#' 
#' @param x The result of a single NPC curve fit from get.monod or similar
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
#' 
plot.npc<-function(x,plot_ci=TRUE,plot_obs=TRUE,xlim=NULL,ylim=NULL,main=NA,fpath=NA,...){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(x)==1){
    x<-x[[1]]
  }
  
  # adjust plotting window to specific curve?
  if(is.null(ylim)) ylim <- c(1.3*min(c(min(x$data$mu,na.rm=T)[1],0)),1.3*max(x$data$mu,na.rm=T)[1])
  if(is.null(xlim)) xlim <- c(0,1.2*max(x$data$nutrients,na.rm=T)[1])
  
  # generate predictions along a range of nutrients
  preds<-predict(x,newdata=data.frame(nutrients=seq(xlim[1],xlim[2],length.out = 200)),se.fit = plot_ci)    
  
  # generate basic plot
  cplot<-ggplot(preds,aes(x=.data$nutrients,y=.data$mu))+
    geom_hline(yintercept = 0)+
    geom_line()+
    scale_x_continuous('Nutrients')+
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
      full.path<-paste(fpath,"NPC_fit_",time,".pdf",sep='')        
    }else{ # or not
      full.path<-paste(fpath,"NPC_fit_",main,".pdf",sep='') 
    }
    
    # save plot at full.path with given file name.
    ggsave(device = grDevices::pdf(),filename = full.path,cplot,width = 5,height = 4)
    grDevices::dev.off()
    invisible(cplot)
  }else{
    return(cplot) 
  }
}

#' logLik extraction method for NPC objects
#' 
#' @param x Object of class npc
#' @param \dots Additional arguments (not used)
#' 
#' @export
logLik.npc<-function(x,...){
  x$logLik
}

#' Fit Monod curve to growth rate vs. nutrient concentration data
#' 
#' Note: this function currently does not use grid.mle2, as the regressions are usually pretty stable; add this 
#' feature later.
#' 
#' @param nutrients Nutrient concentration
#' @param mu Exponential growth rate
#' @param method Specify which fitting algorithm to use, 'mle2' or 'grid.mle2'
#' @param ... Additional arguments passed to grid.mle2 (e.g., control=list(maxit=2000))
#' 
#' @examples 
#' data("example_monod_data")
#' fit<-get.monod(example_monod_data$phosphate,example_monod_data$growth.rate)
#' fit
#' 
#' # not run:
#' #plot(fit)
#' 
#' @export
#' @import emdbook
#' @import ggplot2
get.monod<-function(nutrients,mu,method='mle2',...){
  monod.tmp<-stats::na.omit(data.frame(nutrients,mu))
  nnutr<-length(unique(monod.tmp$nutrients))
  
  if(nnutr<=3){
    print("Caution in get.monod - focal data set has <=3 unique nutrients, risk of overfitting is high!")
  }

  if(method=='grid.mle2'){
    print("Error: this option not yet implemented... 2/13/20")
  }
  
  if(method=='mle2'){
    k.guess <- max(monod.tmp$nutrients)[1]/3
    umax.guess <- max(c(max(monod.tmp$mu)[1],0.01))
    fit<-bbmle::mle2(mu~dnorm(mean=monod_curve(nutrients,umax,k,log.k=TRUE,log.umax=TRUE),sd=exp(s)),
                       start=list(umax=log(umax.guess),k=log(k.guess),s=log(2)),data=monod.tmp)

    # reframe result on non-logged scale
    model_cf<-as.list(coef(fit))
    tmp.cfs<-model_cf
    tmp.cfs$umax<-exp(tmp.cfs$umax)
    tmp.cfs$k<-exp(tmp.cfs$k)
    fit0<-bbmle::mle2(mu~dnorm(mean=monod_curve(nutrients,umax,k,log.k=FALSE,log.umax=FALSE),sd=exp(s)),
                     start=tmp.cfs,control=list(maxit=0),
                     data=monod.tmp)
  }
  
  # pull out parameters on regular scale
  display_cf<-as.list(coef(fit0))
  
  # simple Fisher confidence intervals on regular scale:
  ci<-mleTools::ci.FI(fit0)
  
  # vcov on log scale
  vcov_mat<-vcov(fit)
  
  # calculate R2
  rsqr<-get.R2(predict(fit),monod.tmp$mu)
  
  # create empty object of class 'npc'
  vec<-new_npc()
  
  # populate npc object
  vec$type<-'monod'
  vec$model_cf<-model_cf
  vec$display_cf<-display_cf
  vec$display_ci<-ci
  vec$vcov<-vcov_mat
  vec$nobs<-nrow(monod.tmp)
  vec$nnutr<-nnutr
  vec$rsqr<-rsqr
  vec$logLik<-logLik(fit)
  vec$aic<-stats::AIC(fit)
  vec$data<-monod.tmp
  
  # Finished, return relevant stats.  
  return(vec)
}


#' Fit Monod curve with loss term to growth rate vs. nutrient concentration data
#' 
#' Note: this function currently does not use grid.mle2, as the regressions are usually pretty stable; add this 
#' feature later.
#' 
#' @param nutrients Nutrient concentration
#' @param mu Exponential growth rate
#' @param method Specify which fitting algorithm to use, 'mle2' or 'grid.mle2'
#' @param ... Additional arguments passed to grid.mle2 (e.g., control=list(maxit=2000))
#' 
#' @examples 
#' 
#' @export
#' @import emdbook
#' @import ggplot2
get.monod2<-function(nutrients,mu,method='mle2',...){
  monod.tmp<-stats::na.omit(data.frame(nutrients,mu))
  nnutr<-length(unique(monod.tmp$nutrients))
  
  if(nnutr<=3){
    print("Caution in get.monod - focal data set has <=3 unique nutrients, risk of overfitting is high!")
  }
  
  if(method=='grid.mle2'){
    print("Error: this option not yet implemented... 2/13/20")
  }
  
  if(method=='mle2'){
    
    k.guess <- max(monod.tmp$nutrients)[1]/3
    z.guess <- max(c(0.000001,abs(min(monod.tmp$mu)[1])))
    umax.guess <- max(c(max(monod.tmp$mu)[1],0.01))+z.guess
    fit<-bbmle::mle2(mu~dnorm(mean=monod_curve_type2(nutrients,umax,k,z,log.pars=TRUE),sd=exp(s)),
                       start=list(umax=log(umax.guess),k=log(k.guess),z=log(z.guess),s=log(2)),
                       data=monod.tmp)

    # reframe result on non-logged scale
    model_cf<-as.list(coef(fit))
    tmp.cfs<-model_cf
    tmp.cfs$umax<-exp(tmp.cfs$umax)
    tmp.cfs$k<-exp(tmp.cfs$k)
    tmp.cfs$z<-exp(tmp.cfs$z)
    fit0<-bbmle::mle2(mu~dnorm(mean=monod_curve_type2(nutrients,umax,k,z,log.pars=FALSE),sd=exp(s)),
                      start=tmp.cfs,control=list(maxit=0),
                      data=monod.tmp)
  }
  
  # pull out parameters on regular scale
  display_cf<-as.list(coef(fit0))
  
  # simple Fisher confidence intervals on regular scale:
  ci<-mleTools::ci.FI(fit0)
  
  # vcov on log scale
  vcov_mat<-vcov(fit)
  
  # calculate R2
  rsqr<-get.R2(predict(fit),monod.tmp$mu)
  
  # create empty object of class 'npc'
  vec<-new_npc()
  
  # populate npc object
  vec$type<-'monod2'
  vec$model_cf<-model_cf
  vec$display_cf<-display_cf
  vec$display_ci<-ci
  vec$vcov<-vcov_mat
  vec$nobs<-nrow(monod.tmp)
  vec$nnutr<-nnutr
  vec$rsqr<-rsqr
  vec$logLik<-logLik(fit)
  vec$aic<-stats::AIC(fit)
  vec$data<-monod.tmp
  
  # Finished, return relevant stats.  
  return(vec)
}
