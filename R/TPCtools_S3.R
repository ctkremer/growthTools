


#' Thermal performance curve class creation
#' 
#' @export
new_tpc<-function(){
  tpcObj<-list()
  class(tpcObj)<-"tpc"
  
  # tpc type
  tpcObj$type<-character()
  
  # estimated parameters from fit
  tpcObj$cf<-list()
  tpcObj$cf_ciFI<-matrix()
  tpcObj$vcov<-matrix()
  
  # cardinal traits (could format as table instead?)
  tpcObj$umax<-double()
  tpcObj$umax_ci<-NA
  tpcObj$topt<-double()
  tpcObj$topt_ci<-NA
  tpcObj$tmin<-double()
  tpcObj$tmin_ci<-NA
  tpcObj$tmax<-double()
  tpcObj$tmax_ci<-NA
  
  # fit diagnostics
  tpcObj$rsqr<-double()
  tpcObj$ntemps<-integer()
  tpcObj$nobs<-integer()
  tpcObj$logLik<-double()
  tpcObj$aic<-double()
  
  # data
  tpcObj$data<-data.frame()
  
  tpcObj
}


#' Print method for TPC objects
#' 
#' @param x Object of class tpc
#' @param \dots Additional arguments (not used)
#' 
#' @export
print.tpc<-function(x,...){
  cat("Thermal Performance Curve fit\n")
  cat("Model type: ",x$type,"\n")
  cat("\nEstimated parameters:\n")
  tmp<-cbind(x$cf,x$cf_ciFI)
  colnames(tmp)<-c('parameter',colnames(x$cf_ciFI))
  print(tmp)
  # figure out how to make this look tidier? Maybe model after:
  # Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - 
  #digmin)), digits = digits)
  
  cat("\nFit diagnostics:\n")
  cat("R2 = ",round(x$rsqr,4))
  cat(", logLik = ",x$logLik,"(df = ",attr(x$logLik,'df'),"), nobs = ",x$nobs)
  invisible(x)
}

#' Predict values from a TPC fit
#' 
#' @param object The result of a single TPC curve fit from get.nbcurve.tpc, an object of class tpc
#' @param newdata A new data frame, containing a sequence of `temperature` values at which model predictions should be made; defaults to (-2,40)
#' @param se.fit logical; should standard error values be returned?
#' @param \dots Additional arguments (not used)
#' 
#' @export
#' @method predict tpc
predict.tpc<-function(object,newdata,se.fit=FALSE,...){
  
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
         nbcurve={mu<-nbcurve(newdata$temperature, topt = object$cf$topt, w = object$cf$w, 
                              a = object$cf$a, b = object$cf$b)},
         decurve={mu<-decurve(newdata$temperature, topt = object$cf$topt, b1 = object$cf$b1, 
                              b2 = object$cf$b2, d0 = object$cf$d0, d2 = object$cf$d2)},
         droopcurve={mu<-droopcurve(newdata$temperature, b1 = object$cf$b1, 
                              b2 = object$cf$b2, d0 = object$cf$d0, d1 = object$cf$d1,
                              d2 = object$cf$d2, X = object$cf$X)},
         stop(print("unrecognized tpc model type in predict.tpc!")))
  newdata$mu<-mu
  
  if(se.fit){
    insert<-paste(newdata$temperature,collapse=',')
    
    switch(object$type,
           nbcurve={st<-paste("nbcurve(c(",insert,"),topt,w,a,b)",sep='')},
           decurve={st<-paste("decurve(c(",insert,"),topt,b1,b2,d0,d2)",sep='')},
           droopcurve={st<-paste("droopcurve(c(",insert,"),b1,b2,d0,d1,d2,X)",sep='')},
           stop(print("unrecognized tpc model type in predict.tpc!")))
    dvs0<-suppressWarnings(deltavar2(fun=parse(text=st),meanval=object$cf,Sigma=object$vcov))
    newdata$se.fit<-sqrt(dvs0)
    
    # Better approach? Pass correct local environment to deltavar... BUSTED
    #dvs0<-deltavar3(fun=nbcurve(xs,topt,w,a,b),meanval=cf,Sigma = vcov(fit),cenv=current.env)
    #dvs0<-suppressWarnings(deltavar(fun=nbcurve(xs,topt,w,a,b),meanval=cf,Sigma=vcov(fit)))
  }
  
  return(newdata)
}



#' Plotting function for single thermal performance curve
#' 
#' @param x The result of a single TPC curve fit from get.nbcurve.tpc or similar
#' @param plot_ci logical, should the resulting plot include 95\% confidence bands
#' @param plot_obs logical, should resulting plot include raw data
#' @param xlim x-axis range (temperature)
#' @param ylim y-axis range (adjusts internally from -0.2 to slightly above umax+CI)
#' @param main Character string providing plot title (usually id info for the plotted curve)
#' @param fpath If visual requested, and valid file path provided, save plot as a .pdf file
#' @param \dots Additional arguments (not used)
#' 
#' @export
#' @import ggplot2
plot.tpc<-function(x,plot_ci=TRUE,plot_obs=TRUE,xlim=c(-2,40),ylim=c(-0.2,5),main=NA,fpath=NA,...){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(x)==1){
    x<-x[[1]]
  }
  
  # adjust plotting window to specific curve?
  #ylim<-c(ylim[1],1.1*(x$umax+x$umax_ci[2]))
  
  # generate predictions along a range of temperatures
  preds<-predict(x,newdata=data.frame(temperature=seq(xlim[1],xlim[2],length.out = 100)),se.fit = plot_ci)    
  
  # generate basic plot
  cplot<-ggplot(preds,aes(x=.data$temperature,y=.data$mu))+
    geom_hline(yintercept = 0)+
    geom_line()+
    scale_x_continuous('Temperature')+
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
      full.path<-paste(fpath,"TPC_fit_",time,".pdf",sep='')        
    }else{ # or not
      full.path<-paste(fpath,"TPC_fit_",main,".pdf",sep='') 
    }
    
    # save plot at full.path with given file name.
    ggsave(device = grDevices::pdf(),filename = full.path,cplot,width = 5,height = 4)
    grDevices::dev.off()
    invisible(cplot)
  }else{
    return(cplot) 
  }
}

