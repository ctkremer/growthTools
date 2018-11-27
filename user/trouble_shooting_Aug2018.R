
### Trouble-shooting for Joshua

library(mleTools)
library(growthTools)
library(lubridate)
library(dplyr)

tmp<-read.csv('/Users/colin/Downloads/Dimensions_July_thermal_assay_all.csv')

# format
tmp$Time_Day <- mdy_hm(tmp$Time_Day)
tmp$Dilution <- as.factor(tmp$Dilution)
tmp$Temp <- as.factor(tmp$Temp) 
tmp$Notes <- as.character(tmp$Notes)

# create a difftime column
tmp2<-tmp %>% group_by(Taxonomic_ID,Sample_Code,Replicate,Dilution,Temp) %>% summarise(minD=min(Time_Day),minF=min(Fluor),maxF=max(Fluor))
tmp2<-merge(tmp,tmp2)
tmp2$dtime<-as.numeric(difftime(tmp2$Time_Day,tmp2$minD,units = 'days'))

# remove all rows with notes, these indicate contamination or observations with a problem
tmp3 <- tmp2[tmp2$Notes == '',] # dropping 416 rows
dat <- tmp3

dat$ln.Fluor<-log(dat$Fluor)

# try to run in bulk:
#devtools::load_all()

dat$id<-paste(dat$Sample_Code,dat$Temp,dat$Replicate,dat$Dilution)
str(dat)

dat2<-dat[!is.na(dat$ln.Fluor),]

gdat <- dat %>% group_by(Sample_Code,Temp,Replicate,Dilution) %>% do(grs=get.growth.rate(x=.$dtime,y=.$ln.Fluor,id=.$id,fpath=NA,methods=c('linear'),verbose=F))
gdat


################

# Take 2:

# Calculating growthrates/tools using Colins' package
# updated 8-29-18

library(mleTools) # installed 8-25-18
library(growthTools) # installed 8-25-18
library(dplyr)
library(ggplot2)
library(lubridate)

packageVersion("mleTools") #0.1.0
packageVersion("growthTools") # 0.0.0.9000
packageVersion("dplyr") # 0.7.6
packageVersion("ggplot2") # 3.0.0

# R version
R.Version() # 3.4.4 Someone to lean on

# To view a document outlining how this package works and what it contains, try:
vignette("growthTools_vignette",package="growthTools")

tmp <- read.csv('/Users/colin/Dropbox/DJ_Field_Campaign/DJ_Data/Pico_Syn_Assays/Dimensions_July_thermal_assay_all.csv')
#tmp <- read.csv('/Users/joshuakling/Dropbox/Grad_School/Research/DJ_Field_Campaign/DJ_Data/Pico_Syn_Assays/Dimensions_July_thermal_assay_all.csv')

# format
tmp$Time_Day <- mdy_hm(tmp$Time_Day)
tmp$Dilution <- as.factor(tmp$Dilution)
tmp$Temp <- as.factor(tmp$Temp) 
tmp$Notes <- as.character(tmp$Notes)
tmp$ln.Fluor <- log(tmp$Fluor)

# create a difftime column
tmp2<-tmp %>% group_by(Taxonomic_ID,Sample_Code,Replicate,Dilution,Temp) %>% summarise(minD=min(Time_Day),minF=min(ln.Fluor),maxF=max(ln.Fluor))
tmp2<-merge(tmp,tmp2)
tmp2$dtime<-as.numeric(difftime(tmp2$Time_Day,tmp2$minD,units = 'days'))

# remove all rows with notes, these indicate contamination or observations with a problem
tmp3 <- tmp2[tmp2$Notes == '',] # dropping 416 rows

dat <- tmp3

# test run with a single isolate
t.dat <- dat[dat$Sample_Code == 'LA3' & dat$Replicate == 'A' & dat$Temp == 29,]
head(t.dat)

t.res <- get.growth.rate(t.dat$dtime,t.dat$ln.Fluor,plot.best.Q = T,id = 'LA3 29',methods = 'linear')
t.res

# run on all files
gdat <- dat %>% group_by(Sample_Code,Temp,Replicate,Dilution) %>% 
  do(grs=get.growth.rate(x=.$dtime,y=.$ln.Fluor,fpath=NA,methods= c('linear')))
#%>% gdat

# get mean growth rates for each strain at each temperature
tmp <- gdat %>% group_by(Sample_Code,Temp,Replicate,Dilution) %>% 
  summarise(mu=as.numeric(unlist(grs)[1]),best.model=as.character(unlist(gdat$grs)[3]))
head(tmp)

# averaging across dilution periods
rdy_for_fits <- tmp %>% group_by(Sample_Code,Temp,Replicate) %>% 
  summarise(means=mean(mu,na.rm = T)) %>% ungroup()
head(rdy_for_fits)

# Set directory for plots
#path <- '/Users/colin/temp_plots/'
path <- '/Users/joshuakling/Dropbox/Grad_School/Research/DJ_Field_Campaign/DJ_Data/Pico_Syn_Assays/Plots_from_Colins_nbc_package/'
path

rdy_for_fits$Temp2 <- as.numeric(as.character(rdy_for_fits$Temp))

length(unique(rdy_for_fits$Sample_Code))
ggplot(rdy_for_fits,aes(x=Temp2,y=means))+
  geom_point()+
  geom_hline(yintercept = 0)+
  facet_wrap(~Sample_Code)+
  theme_bw()

###############################################################

### both of the below produce errors ###

###############################################################

data.frame(rdy_for_fits) %>% filter(Sample_Code=='LA103') %>% do(tpcs=get.nbcurve.tpc(.$Temp2,.$means,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=path,id=.$Sample_Code))

tdat<-data.frame(rdy_for_fits) %>% filter(Sample_Code=='LA103')
temp<-tdat$Temp2
mu<-tdat$means

write.csv(x = tdat,'/Users/colin/Desktop/tmp_tdat.csv',row.names=F)


plot(mu~temp)
curve(nbcurve2(x,opt = 23,w = 12,a = log(0.12),b=0.05),0,35,col='red',add=T)


nbc_dat <- rdy_for_fits %>% group_by(Sample_Code) %>% 
  do(tpcs=get.nbcurve.tpc(.$Temp2,.$means,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=path,id=.$Sample_Code))

res2 <- nbc_dat %>% summarise(Sample_Code,topt=tpcs$o,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr,a=exp(tpcs$a),b=tpcs$b,w=tpcs$w,tmin.lw=)
data.frame(res2)

nbc_dat[nbc_dat$Sample_Code=='LA117',2][[1]]

save(file = "/Users/colin/Desktop/josh_fits.RData",nbc_dat)



