


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
#' @export
print.tpc<-function(object){
  cat("Thermal Performance Curve fit\n")
  cat("Model type: ",object$type,"\n")
  cat("\nEstimated parameters:\n")
  tmp<-cbind(object$cf,object$cf_ciFI)
  colnames(tmp)<-c('parameter',colnames(object$cf_ciFI))
  print(tmp)
  # figure out how to make this look tidier? Maybe model after:
  # Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - 
  #digmin)), digits = digits)
  
  cat("\nFit diagnostics:\n")
  cat("R2 = ",round(object$rsqr,4))
  cat(", logLik = ",object$logLik,"(df = ",attr(object$logLik,'df'),"), nobs = ",object$nobs)
  invisible(object)
}

#' Predict values from a TPC fit
#' 
#' @param fit The result of a single TPC curve fit from get.nbcurve.tpc
#' @param newdata A new data frame, containing a sequence of `temperature` values at which model predictions should be made; defaults to (-2,40)
#' @param se.fit logical; should standard error values be returned?
#' 
#' @export
#' @method predict tpc
predict.tpc<-function(fit,newdata,se.fit=FALSE){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(fit)==1){
    fit<-fit[[1]]
  }
  
  # instead of providing a default newdata, instead return values from observed 
  # temperatures as default if newdata is missing
  # should also check format of newdata is correct...
  if (missing(newdata) || is.null(newdata)) {
    newdata<-fit$data
  }
  
  # generate predictions across a range of temperatures
  switch(fit$type,
         nbcurve={mu<-nbcurve(newdata$temperature, topt = fit$cf$topt, w = fit$cf$w, 
                              a = fit$cf$a, b = fit$cf$b)},
         decurve={mu<-decurve(newdata$temperature, topt = fit$cf$topt, b1 = fit$cf$b1, 
                              b2 = fit$cf$b2, d0 = fit$cf$d0, d2 = fit$cf$d2)},
         stop(print("unrecognized tpc model type in predict.tpc!")))
  newdata$mu<-mu
  
  if(se.fit){
    insert<-paste(newdata$temperature,collapse=',')
    
    switch(fit$type,
           nbcurve={st<-paste("nbcurve(c(",insert,"),topt,w,a,b)",sep='')},
           decurve={st<-paste("decurve(c(",insert,"),topt,b1,b2,d0,d2)",sep='')},
           stop(print("unrecognized tpc model type in predict.tpc!")))
    dvs0<-suppressWarnings(deltavar2(fun=parse(text=st),meanval=fit$cf,Sigma=fit$vcov))
    newdata$se.fit<-sqrt(dvs0)
    
    # Better approach? Pass correct local environment to deltavar... BUSTED
    #dvs0<-deltavar3(fun=nbcurve(xs,topt,w,a,b),meanval=cf,Sigma = vcov(fit),cenv=current.env)
    #dvs0<-suppressWarnings(deltavar(fun=nbcurve(xs,topt,w,a,b),meanval=cf,Sigma=vcov(fit)))
  }
  
  return(newdata)
}



#' Plotting function for single thermal performance curve
#' 
#' @param fit The result of a single TPC curve fit from get.nbcurve.tpc or similar
#' @param plot.ci logical, should the resulting plot include 95\% confidence bands
#' @param plot.obs logical, should resulting plot include raw data
#' @param xlim x-axis range (temperature)
#' @param ylim y-axis range (adjusts internally from -0.2 to slightly above umax+CI)
#' @param main Character string providing plot title (usually id info for the plotted curve)
#' @param fpath If visual requested, and valid file path provided here, plot will be saved as a .pdf file. Default is NA.
#' 
#' @export
#' @import ggplot2
plot.tpc<-function(fit,plot_ci=TRUE,plot_obs=TRUE,xlim=c(-2,40),ylim=c(-0.2,5),main=NA,fpath=NA){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(fit)==1){
    fit<-fit[[1]]
  }
  
  # adjust plotting window to specific curve?
  #ylim<-c(ylim[1],1.1*(fit$umax+fit$umax_ci[2]))
  
  # generate predictions along a range of temperatures
  preds<-predict(fit,newdata=data.frame(temperature=seq(xlim[1],xlim[2],length.out = 100)),se.fit = plot_ci)    
  
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
    cplot<-cplot+geom_point(data=fit$data)
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
  }
  
  return(cplot) 
}

