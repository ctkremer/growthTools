
# Need to design some better graphing functions, which can take key outputs from 
# get.nbcurve.tpc() and generate nice-looking, easy to combine or customize
# plots of tpc functions...

nbcurve.traits<-get.nbcurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',
                                plotQ=F,conf.bandQ = F,fpath=NA)
#get.nbcurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',
#                plotQ=T,conf.bandQ = T,fpath=NA)

# take the tpcs entry from:
nb.res <- sp1b %>% group_by(isolate.id,dilution) %>% 
  do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',
                          plotQ=F,conf.bandQ=F,fpath=NA,id=.$dilution))


### Alternative? develop a predict.nbcurve function

# Predict function:
predict.nbcurve<-function(fit.info,newdata=data.frame(temperature=seq(-2,40,0.1)),se.fit=FALSE){
  
  # Check level of nesting for fit.info, and reduce if necessary
  if(length(fit.info)==1){
    fit.info<-fit.info[[1]]
  }
  
  # generate predictions across a range of temperatures
  mu<-nbcurve(newdata$temperature, topt = fit.info$topt, w = fit.info$w, a = fit.info$a, b = fit.info$b)
  newdata$mu<-mu
  
  if(se.fit){
    st<-paste("nbcurve(c(",paste(ts,collapse=','),"),topt,w,a,b)",sep='')
    dvs0<-suppressWarnings(deltavar2(fun=parse(text=st),meanval=fit.info$cf,Sigma=fit.info$vcov))
    newdata$se.fit<-sqrt(dvs0)  
  }
  
  return(newdata)
}


# generate predictions for a single curve:
fit.info<-nb.res[1,]$tpcs
tmp<-predict.nbcurve(fit.info,se.fit=T)
head(tmp)

# generate predictions for multiple curves:
tmp2 <- nb.res %>% group_by(isolate.id,dilution) %>% do(predict.nbcurve(.$tpcs,se.fit=T))
head(tmp2)

# plot multiple curves simultaneously
ggplot(tmp2,aes(x=temperature,y=mu))+
  geom_hline(yintercept = 0)+
  geom_line(aes(colour=factor(dilution)))+
  coord_cartesian(ylim=c(-0.2,1),xlim=c(-2,35))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)')+
  theme_bw()


#fit.infoX<-fit.info
#fit.info<-fit.infoX

# Norberg curve plotting function

plot.nbcurve.tpc<-function(fit.info,plot.ci=TRUE,plot.obs=TRUE,xlim=c(-2,40),ylim=c(-0.2,5)){
  
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
  cplot<-ggplot(preds,aes(x=temperature,y=mu))+
    geom_hline(yintercept = 0)+
    geom_line()+
    coord_cartesian(xlim=xlim,ylim=ylim)+
    theme_bw()

  # add confidence bands?
  if(plot.ci){
    cplot<-cplot+geom_ribbon(aes(ymin=mu-1.96*se.fit,ymax=mu+1.96*se.fit),alpha=0.2)
  }
  
  # add observations?
  if(plot.obs){
    cplot<-cplot+geom_point(data=fit.info$data)
  }
  
  return(cplot) 
}


# use function to generate a single curve plot:
c1<-plot.nbcurve.tpc(fit.info)
c1

# demonstrate changes in layout
c1+scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)')+
  theme(panel.grid = element_blank())




### Plot multiple curves with observations?


# generate predictions for multiple curves:
tmp2 <- nb.res %>% group_by(isolate.id,dilution) %>% do(predict.nbcurve(.$tpcs,se.fit=T))
head(tmp2)

# recapitulate corresponding data frame of observations:
tmp2.obs <- nb.res %>% group_by(isolate.id,dilution) %>% do(.$tpcs[[1]]$data)
head(tmp2.obs)

# plot multiple curves simultaneously
ggplot(tmp2,aes(x=temperature,y=mu))+
  geom_hline(yintercept = 0)+
  geom_line(aes(colour=factor(dilution)))+
  geom_point(data=tmp2.obs,aes(colour=factor(dilution)),alpha=0.2)+
  coord_cartesian(ylim=c(-0.2,1),xlim=c(-2,35))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)')+
  theme_bw()
