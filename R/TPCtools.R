
#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Conduct fits of thermal reaction norms to various parametric equations

# - currently focus is on Norberg fits
# - expand to include double exponential
# - expand using dplyr approaches to make automation over multiple species easier
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
#' growth rate. See Thomas et al. 2012 for more on this, and \code{nbcurve2()} for an
#' alternative formulation.
#' 
#' @param x Temperature
#' @param opt (Competitive) Optimum temperature
#' @param w Thermal niche width
#' @param a ln(Exponential intercept)
#' @param b Exponential scaling
#' 
#' @return Predicted exponential growth rate at temperature x
#' 
#' @export
nbcurve<-function(x,opt,w,a,b){
  res<-exp(a)*exp(b*x)*(1-((x-opt)/(w/2))^2)
  res
}

#' Re-parameterized Norberg functions
#' 
#' With this formulation, optimum temperature is explicitly the temperature at which
#' a species achieves its highest growth rate.
#' 
#' @param x Temperature
#' @param opt Optimum temperature
#' @param w Thermal niche width
#' @param a ln(Exponential intercept)
#' @param b Exponential scaling
#' 
#' @return Predicted exponential growth rate at temperature x
#' 
#' @export
nbcurve2<-function(x,opt,w,a,b){
  num <- -2-2*b*opt+2*b*x+sqrt(4+(b*w)^2)
  res<-exp(a)*exp(b*x)*(1-(num/(b*w))^2)
  res
}

#' Get true optimum temperature from nbcurve()
#' 
#' Given the competitive optimum temperature and other Norberg parameters, calculate the
#' true optimum temperature.
#' 
#' @param topt Competitive optimum temperature
#' @param w Thermal niche width
#' @param b Exponential scaling
#' 
#' @return Optimum temperature
#' 
#' @export
get.topt<-function(topt,w,b){
  (1/(2*b))*(-2+2*b*topt+sqrt(4+b^2*w^2))
}

#' Calculate thermal niche limits
#' 
#' Note: check whether this operates on nbcurve() or nbcurve2() parameters!!!
#' 
#' @param topt Competitive optimum temperature
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


#' Fit Norberg curve to growth rate vs. temperature data
#' 
#' @param temp Temperature
#' @param mu Exponential growth rate
#' @param method Specify which fitting algorithm to use, 'mle2' or 'grid.mle2'
#' @param plotQ Should regression be visualized?
#' @param fpath If visual requested, and valid file path provided here, plot will be saved as a .pdf file. Default is NA.
#' 
#' @export
#' @import bbmle
#' @import mleTools
#' @import emdbook
get.nbcurve.tpc<-function(temp,mu,method='grid.mle2',plotQ=F,fpath=NA){
  tmp<-data.frame(mu,temp)
  
  if(method=='grid.mle2'){
    
    # set up search of a grid of parameter guesses
    grids<-list(o=seq(15,25,5),w=seq(15,40,5))
    start<-list(o=NA,w=NA,a=-1.11,b=0.05,s=10)
    
    fit0<-grid.mle2(minuslogl=mu~dnorm(mean=nbcurve2(temp,o,w,a,b),sd=s),
                    grids=grids,start=start,data=tmp)
    cfg<-coef(fit0$res.best) # this seemed to be throwing problems b/c of an issue with accessing mle2...?

    # polish best fit model, using formual interface:
    guesses<-as.list(cfg)
    fit<-mle2(mu~dnorm(mean=nbcurve2(temp,o,w,a,b),sd=s),
              start=guesses,data=tmp)
  }
  
  if(method=='mle2'){
    o.guess <- tmp$temp[tmp$mu==max(tmp$mu)]
    w.guess <- diff(range(tmp$temp))
    a.guess <- -1.11
    b.guess <- 0.05
    fit<-mle2(mu~dnorm(mean=nbcurve2(temp,o,w,a,b),sd=s),
              start=list(o=o.guess,w=w.guess,a=a.guess,b=b.guess,s=2),data=tmp)
  }
  
  # pull out parameters
  cf<-as.list(coef(fit))
#  cf$a<-exp(cf$a)
  
  # additional responses/traits:
  tmin<-get.tlim(cf$o,cf$w,cf$b,type='tmin')
  tmax<-get.tlim(cf$o,cf$w,cf$b,type='tmax')
  
  # calculate R2
  rsqr<-get.R2(fit,tmp$mu)

  # Calculate confidence intervals using the delta method (see Bolker book, pg 255)
  #   ... why is this useful? really, we'd like CI's on the x-axis position of these traits
  #focal.dvs<-suppressWarnings(deltavar(fun=nbcurve2(c(cf$o,tmin,tmax),o,w,a,b),meanval=cf,Sigma=vcov(fit)))
  #focal.ci<-1.96*sqrt(focal.dvs)
  
  # save output:
  #vec<-c(cf,rsqr,tmin,tmax,focal.ci)
  vec<-as.list(c(cf,rsqr=rsqr,tmin=tmin,tmax=tmax))
  
  # Plot results:
  if(plotQ){
    xs<-seq(min(tmp$temp),max(tmp$temp),0.1)
    new.data<-data.frame(temp=xs)
    new.data$mu<-predict(fit,newdata=new.data)
    
    # confidence bands via the delta method (see Bolker book, pg. 255)
    #dvs<-suppressWarnings(deltavar(fun=nbcurve2(xs,o,w,a,b),meanval=cf,Sigma=vcov(fit)))
    #new.data$ml.ci<-1.96*sqrt(dvs)
    
    # plot it out
    p1<-ggplot(tmp,aes(x=temp,y=mu))+
        geom_point()+
        #geom_ribbon(data=new.data,aes(ymin=mu-ml.ci,ymax=mu+ml.ci),alpha=0.2)+
        geom_line(data=new.data)+
        geom_hline(yintercept = 0)+
        theme_bw()
    if(!is.na(fpath)){
      ggsave(paste(fpath,"TPC_fit.pdf",sep=""),p1)
    }else{
      print(p1)
    }
  }

  # Finished, return relevant stats.  
  return(vec)
}

#' Calculate R2 for an mle2 model
#' 
#' Assuming a normal error distribution and an mle2 fit that uses the 'formula' syntax
#' (this function will not work if mle2 is handed a black-boxed NLL function)
#' 
#' @param model An mle2 model object
#' @param obs Observed y values
#' 
#' @return R2 value
#' 
#' @export
get.R2<-function(model,obs){
  expc<-predict(model)
  rsqr<-1-sum((obs-expc)^2)/sum((obs-mean(obs))^2)
  rsqr
}



