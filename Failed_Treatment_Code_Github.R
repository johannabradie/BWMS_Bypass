rm(list=ls())
set.seed(2935)

## This code loads required data and generates necessary distributions.  It is done within its own code to keep the code cleaner. 
source("../Bradie_etal2021_Curves.R")

## Identify column numbers for arrival and source environemetnal information 
rtrips_Arr_cols=35:37
rtrips_Source_cols=40:42

library('plyr')

########## FUNCTIONS #######
# Load custom R functions ### 

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

meanfun <- function (data, i){
  d<-data [i,]
  return (mean (d))   }

euc_d<-function(x,y) {
  dist=sum((x-y)^2)^0.5
  return(dist)
}

estab<-function(alpha,N,c=1) {
  est=1-exp(1)^-(alpha*N^c)
  return(est)
}

establish<-function(est) {
  random.probs=runif(length(est),0,1)
  establishments=length(which(est>random.probs))
  return(establishments)}


########  DATA READ-IN AND SET-UP  #####################
sal_sens=1  ## Sets baseline for salinity sensitivity analysis 
date="2022_03-17"


## This iterates through different sensitivity analyses.  s=1 is the main analysis, with the other s describing different sensitivity analyses that were run.
## Normally just run s=1 and then run the s=2-->9 once everything was finished with the main analyses. 
for (s in 1:1) {
  if (s==1) {
    current_date=paste(date,'_s',s,sep="")
    percent.reduction=0.95
    alpha_par1=0.00005
    alpha_par2=5
    fail_threshold=10
    
  }
  
  if (s==2) {
    current_date=paste(date,'_s',s,sep="")
    percent.reduction=0.95
    alpha_par1=0.0001
    alpha_par2=5
    fail_threshold=10
  } 
  
  if (s==3) {
    current_date=paste(date,'_s',s,sep="")
    percent.reduction=0.95
    alpha_par1=0.0005
    alpha_par2=5
    fail_threshold=10
  }

if (s==4) {
  current_date=paste(date,'_s',s,sep="")
  percent.reduction=0.95
  alpha_par1=0.00005
  alpha_par2=5
  fail_threshold=25
}

if (s==5) {
  current_date=paste(date,'_s',s,sep="")
  percent.reduction=0.995
  alpha_par1=0.00005
  alpha_par2=5
  fail_threshold=10
}
if (s==6) {
  current_date=paste(date,'_s',s,sep="")
  percent.reduction=0.9
  alpha_par1=0.00005
  alpha_par2=5
  fail_threshold=10
}
  if (s==7) {
    current_date=paste(date,'_s',s,sep="")
    percent.reduction=0.95
    alpha_par1=0.00005
    alpha_par2=5
    fail_threshold=10
    sal_sens=0.75
  }
  if (s==8) {
    current_date=paste(date,'_s',s,sep="")
    percent.reduction=0.95
    alpha_par1=0.00005
    alpha_par2=5
    fail_threshold=10
    sal_sens=2
  }
  if (s==9) {
    current_date=paste(date,'_s',s,sep="")
    percent.reduction=0.95
    alpha_par1=0.00005
    alpha_par2=5
    fail_threshold=10
    sal_sens=5
  }

ntrips=10000 # NUmber of trips in each simulation 

dens_change=0 #This remains from a different project where we sensitivity tested changes in tank concentratinons. 
# alpha_sensitivity=FALSE
c=1
# current_date=Sys.Date()
# current_date="2021-12-15alpha5e-5.5"
# percent.reduction=0.995

exc_Atlantic=exchange_data[which(exchange_data$Region=="Atlantic"),]  # generate exchange locations available for ships transiting each region 
exc_Pacific=exchange_data[which(exchange_data$Region=="Pacific"),]
exc_GL=exchange_data[which(exchange_data$Region=="GLSLR"),]


for (p in 1:3) {  # This loop covers each of the different pathways we considered
  if (p==1) {pathway="Coast.Int"}
  if (p==2) {pathway="GL.Int"}
  if (p==3) {pathway="GL.Dom"}

for (t in 1:100) {  # We ran the simulation 100 times with 10,000 trips in each because running 1,000,000 at once caused memory issues.  All results were saved and then merged in the analysis script. 
  
if (pathway=="GL.Int") { # Each of these loads appropriate data per pathway
  trips=subset(rtrips[rtrips$ArrivalRegionCode==3 & rtrips$SourceCode==7,])  # Only consider ships transiting the right pathway
  exc_data=which(exchange_data$Region=="GLSLR") # Use exchange data for ships going to that region 

  for (i in 1:length(trips[,1])) { # Calculate envt distance between source and arrival regions 
      trips$ENV_DIST[i]=euc_d(x=trips[i,rtrips_Arr_cols], y=rtrips[i,rtrips_Source_cols])}
    
    rows=sample(1:length(trips[,1]),size=ntrips,replace=TRUE)  # sample with replacement from historical data to generate trips for this simulation 
    trips=trips[rows,]  #trips table is simulated trip data 
    GL_exchange=sample(1:length(exc_GL[,1]),size=ntrips,replace=T) # simulate exchange locations from historical exchange data 
    
    # CALCULATE environmental distance values between exchange location and destination ports using simulated exchange location data
    trips$ENV_DIST_EXC=NA
    exc_data=trips[,rtrips_Arr_cols]
    exc_data$Exc.Scaled.Min.Temp=NA
    exc_data$Exc.Scaled.Max.Temp=NA
    exc_data$Exc.Scaled.Annual.Temp=NA
    exc_data[,4:6]=exc_GL[GL_exchange,14:16]
    
    for (i in 1:length(exc_data[,1])) {
      trips$ENV_DIST_EXC[i]=euc_d(x=exc_data[i,1:3],y=exc_data[i,4:6])
    }
    
    trips=trips[,c(6,22,23,44:46,34)]  # drop columns we don't need for this analysis 
    
    trips$Exc_Type=sample(c("ER","FT"),size=ntrips,replace=T,prob=c((926/(926+545)),(545/(926+545))))  #determine exchange type based on TC database prevalence for this pathway 
    
    # select concentration of ZP in tank sample based on theoretical curve 
    trips$sampledens.ZP<-rnbinom(length(trips$ArrivalRegion=="GLSLR"), size = glisz.ZP, mu = (glimu.ZP*(1+dens_change))) 
    trips$sampledens.ZP[which(trips$sampledens.ZP>GLInt.maxZP)]=GLInt.maxZP
    #next, using the same jurisdictions, determine the proportion of sample densities that are NIS, and sample for each trip#
    trips$propNIS.ZP<-rbeta(length(trips$ArrivalRegion=="GLSLR"), shape1 = glishape1.ZP, shape2 = (glishape2.ZP))
    
    ## Same, but for phytopankton 
    trips$sampledens.phyto<-rnbinom(length(trips$ArrivalRegion=="GLSLR"), size = glisz.phyto, mu = (glimu.phyto*(1+dens_change)))
    trips$sampledens.phyto[which(trips$sampledens.phyto>GLInt.maxphyto)]=GLInt.maxphyto
    #next, using the same jurisdictions, determine the proportion of sample densities that are NIS, and sample for each trip#
    trips$propNIS.phyto<-rbeta(length(trips$ArrivalRegion=="GLSLR"), shape1 = glishape1.phyto, shape2 = (glishape2.phyto))
    
    # Determine sample concentration after exchange by sampling a new concentratino value from theoretical curve 
    trips$excdens.ZP<-rnbinom(length(trips$ArrivalRegion=="GLSLR"), size = glisz.ZP, mu = (glimu.ZP*(1+dens_change)))
    trips$excdens.ZP[which(trips$excdens.ZP>GLInt.maxZP)]=GLInt.maxZP
    # Same, for phyto
    trips$excdens.phyto<-rnbinom(length(trips$ArrivalRegion=="GLSLR"), size = glisz.phyto, mu = (glimu.phyto*(1+dens_change)))
    trips$excdens.phyto[which(trips$excdens.phyto>GLInt.maxphyto)]=GLInt.maxphyto
    
    ## Determine colonization pressure for ZP and phyto from historical sample data
    GLI_CP.ZP=nisdens.ZP$CPNIS[which(nisdens.ZP$ArrivalRegionCode==3 & nisdens.ZP$SourceRegionCode==7)]
    trips$NIS_CP.ZP=rep(NA,length=ntrips)
    trips$NIS_CP.ZP=sample(GLI_CP.ZP,size=length(trips[,1]),replace = TRUE)
    GLI_CP.phyto=nisdens.phyto$CPNIS[which(nisdens.phyto$ArrivalRegionCode==3 & nisdens.phyto$SourceRegionCode==7)]
    trips$NIS_CP.phyto=sample(GLI_CP.phyto,size=ntrips,replace = TRUE)
    
  }

if (pathway=="GL.Dom") {  # next two repeat as above but for different pathways. 
  trips=subset(rtrips[rtrips$ArrivalRegionCode==3 & rtrips$SourceCode!=7,])
  exc_data=which(exchange_data$Region=="GLSLR")
  
    for (i in 1:length(trips[,1])) {
      trips$ENV_DIST[i]=euc_d(x=trips[i,rtrips_Arr_cols], y=rtrips[i,rtrips_Source_cols])}
    
    rows=sample(1:length(trips[,1]),size=ntrips,replace=TRUE)
    trips=trips[rows,]
    
    GL_exchange=sample(1:length(exc_GL[,1]),size=ntrips,replace=T)
    
    # CALCULATE environmental distance values between exchange location and destination ports using simulated exchange location data
    trips$ENV_DIST_EXC=NA
    exc_data=trips[,rtrips_Arr_cols]
    exc_data$Exc.Scaled.Min.Temp=NA
    exc_data$Exc.Scaled.Max.Temp=NA
    exc_data$Exc.Scaled.Annual.Temp=NA
    exc_data[,4:6]=exc_GL[GL_exchange,14:16]
    
    
    for (i in 1:length(exc_data[,1])) {
      trips$ENV_DIST_EXC[i]=euc_d(x=exc_data[i,1:3],y=exc_data[i,4:6])
    }
    
    trips=trips[,c(6,22,23,44:46,34)]
    
    trips$Exc_Type=sample(c("ER","FT"),size=ntrips,replace=T,prob=c((926/(926+545)),(545/(926+545))))
    
    trips$sampledens.ZP<-rnbinom(length(trips$ArrivalRegion=="GLSLR"), size = glglsz.ZP, mu = (glglmu.ZP*(1+dens_change)))
    trips$sampledens.ZP[which(trips$sampledens.ZP>GLDom.maxZP)]=GLDom.maxZP
    
    #next, using the same jurisdictions, determine the proportion of sample densities that are NIS, and sample for each trip#
    trips$propNIS.ZP<-rbeta(length(trips$ArrivalRegion=="GLSLR"), shape1 = gldomshape1.ZP, shape2 = (gldomshape2.ZP))
    
    
    trips$sampledens.phyto<-rnbinom(length(trips$ArrivalRegion=="GLSLR"), size = glisz.phyto, mu = (glimu.phyto*(1+dens_change)))
    trips$sampledens.phyto[which(trips$sampledens.phyto>GLInt.maxphyto)]=GLInt.maxphyto
    
    #next, using the same jurisdictions, determine the proportion of sample densities that are NIS, and sample for each trip#
    trips$propNIS.phyto<-rbeta(length(trips$ArrivalRegion=="GLSLR"), shape1 = glishape1.phyto, shape2 = (glishape2.phyto))
    
    ## Generate exchange densities for trips that are considered high propagule load ##
    trips$excdens.ZP<-rnbinom(length(trips$ArrivalRegion=="GLSLR"), size = glglsz.ZP, mu = (glglmu.ZP*(1+dens_change)))
    trips$excdens.ZP[which(trips$excdens.ZP>GLDom.maxZP)]=GLDom.maxZP
    
    trips$excdens.phyto<-rnbinom(length(trips$ArrivalRegion=="GLSLR"), size = glisz.phyto, mu = (glimu.phyto*(1+dens_change)))
    trips$excdens.phyto[which(trips$excdens.phyto>GLInt.maxphyto)]=GLInt.maxphyto
    
    GLGL_CP.ZP=nisdens.ZP$CPNIS[which(nisdens.ZP$ArrivalRegionCode==3 & nisdens.ZP$SourceRegionCode!=7)]
    trips$NIS_CP.ZP=rep(NA,length=ntrips)
    trips$NIS_CP.ZP=sample(GLGL_CP.ZP,size=length(trips[,1]),replace = TRUE)
    
    GLGL_CP.phyto=nisdens.phyto$CPNIS[which(nisdens.phyto$ArrivalRegionCode==3 & nisdens.phyto$SourceRegionCode==7)]
    trips$NIS_CP.phyto=sample(GLGL_CP.phyto,size=ntrips,replace = TRUE)
    
  }
  

if (pathway=="Coast.Int") {
  rtrips_1=subset(rtrips[rtrips$ArrivalRegionCode==2 & rtrips$SourceCode==7,])
  rtrips_2=subset(rtrips[rtrips$ArrivalRegionCode==4 & rtrips$SourceCode==7,])
  trips=rbind(rtrips_1,rtrips_2)
  rm(rtrips_1,rtrips_2)
  
  for (i in 1:length(trips[,1])) {
    trips$ENV_DIST[i]=euc_d(x=trips[i,rtrips_Arr_cols], y=rtrips[i,rtrips_Source_cols])}
  
  rows=sample(1:length(trips[,1]),size=ntrips,replace=TRUE)
  
  trips=trips[rows,]
  
  Atl_rows=which(trips$ArrivalRegion=="Atlantic")
  Pac_rows=which(trips$ArrivalRegion=="Pacific")
  
  Atl_exchange=sample(1:length(exc_Atlantic[,1]),size=length(Atl_rows),replace=T)
  Pac_exchange=sample(1:length(exc_Pacific[,1]),size=length(Pac_rows),replace=T)
  
  # CALCULATE environmental distance values between exchange location and destination ports using simulated exchange location data
  trips$ENV_DIST_EXC=NA
  exc_data=trips[,rtrips_Arr_cols]
  exc_data$Exc.Scaled.Min.Temp=NA
  exc_data$Exc.Scaled.Max.Temp=NA
  exc_data$Exc.Scaled.Annual.Temp=NA
  exc_data[Atl_rows,4:6]=exc_Atlantic[Atl_exchange,14:16]
  exc_data[Pac_rows,4:6]=exc_Pacific[Pac_exchange,14:16]
  
  for (i in 1:length(exc_data[,1])) {
    trips$ENV_DIST_EXC[i]=euc_d(x=exc_data[i,1:3],y=exc_data[i,4:6])
  }
  
  trips=trips[,c(6,22,23,44:46,34)]
  
  trips$Exc_Type=NA
  trips$Exc_Type[Pac_rows]=sample(c("ER","FT"),size=length(Pac_rows),replace=T,prob=c((5706/(5706+6405)),(6405/(5706+6405))))
  trips$Exc_Type[Atl_rows]=sample(c("ER","FT"),size=length(Atl_rows),replace=T,prob=c((2874/(2874+2100)),(2100/(2874+2100))))
  
  for (i in 1:length(trips[,1])) {
    trips$sampledens.ZP[i]<-ifelse(trips$ArrivalRegion[i]=="Pacific", rnbinom(1, size = pacsz.ZP, mu = (pacmu.ZP*(1+dens_change))), 
                              rnbinom(1, size = atlsz.ZP, mu = (atlmu.ZP*(1+dens_change))))}
    trips$sampledens.ZP[which(trips$sampledens.ZP>PacInt.maxZP & trips$ArrivalRegion=="Pacific")]=PacInt.maxZP
    trips$sampledens.ZP[which(trips$sampledens.ZP>AtlInt.maxZP & trips$ArrivalRegion=="Atlantic")]=AtlInt.maxZP
    
  #next, using the same jurisdictions, determine the proportion of sample densities that are NIS, and sample for each trip#
  trips$propNIS.ZP<-ifelse(trips$ArrivalRegion=="Pacific", rbeta(length(trips$ArrivalRegion=="Pacific"), shape1 = pacshape1.ZP, shape2 = (pacshape2.ZP)), 
                      rbeta(length(trips$ArrivalRegion=="Atlantic"), shape1 = atlshape1.ZP, shape2 = (atlshape2.ZP)))


  for (i in 1:length(trips[,1])) {
      trips$sampledens.phyto[i]<-ifelse(trips$ArrivalRegion[i]=="Pacific", rnbinom(1, size = pacsz.phyto, mu = (pacmu.phyto*(1+dens_change))), 
                                     rnbinom(1, size = atlsz.phyto, mu = (atlmu.phyto*(1+dens_change))))}
      trips$sampledens.phyto[which(trips$sampledens.phyto>PacInt.maxphyto & trips$ArrivalRegion=="Pacific")]=PacInt.maxphyto
      trips$sampledens.phyto[which(trips$sampledens.phyto>AtlInt.maxphyto & trips$ArrivalRegion=="Atlantic")]=AtlInt.maxphyto
      
  ## Generate exchange densities for trips that are considered high propagule load ##
  
  trips$excdens.ZP<-ifelse(trips$ArrivalRegion=="Pacific", rnbinom(length(trips$ArrivalRegion=="Pacific"), size = pacsz.ZP, mu = (pacmu.ZP*(1+dens_change))), 
                              rnbinom(length(trips$ArrivalRegion=="Atlantic"), size = atlsz.ZP, mu = (atlmu.ZP*(1+dens_change))))
  trips$excdens.ZP[which(trips$excdens.ZP>PacInt.maxZP & trips$ArrivalRegion=="Pacific")]=PacInt.maxZP
  trips$excdens.ZP[which(trips$excdens.ZP>AtlInt.maxZP & trips$ArrivalRegion=="Atlantic")]=AtlInt.maxZP
  
  trips$excdens.phyto<-ifelse(trips$ArrivalRegion=="Pacific", rnbinom(length(trips$ArrivalRegion=="Pacific"), size = pacsz.phyto, mu = (pacmu.phyto*(1+dens_change))), 
                                 rnbinom(length(trips$ArrivalRegion=="Atlantic"), size = atlsz.phyto, mu = (atlmu.phyto*(1+dens_change))))
  trips$excdens.phyto[which(trips$excdens.phyto>PacInt.maxphyto & trips$ArrivalRegion=="Pacific")]=PacInt.maxphyto
  trips$excdens.phyto[which(trips$excdens.phyto>AtlInt.maxphyto & trips$ArrivalRegion=="Atlantic")]=AtlInt.maxphyto
  
  
  #next, using the same jurisdictions, determine the proportion of sample densities that are NIS, and sample for each trip#
  trips$propNIS.phyto<-ifelse(trips$ArrivalRegion=="Pacific", rbeta(length(trips$ArrivalRegion=="Pacific"), shape1 = pacshape1.phyto, shape2 = (pacshape2.phyto)), 
                           rbeta(length(trips$ArrivalRegion=="Atlantic"), shape1 = atlshape1.phyto, shape2 = (atlshape2.phyto)))
  
  Atl_CP.ZP=nisdens.ZP$CPNIS[which(nisdens.ZP$ArrivalRegionCode==2 & nisdens.ZP$SourceRegionCode==7)]
  Pac_CP.ZP=nisdens.ZP$CPNIS[which(nisdens.ZP$ArrivalRegionCode==4 & nisdens.ZP$SourceRegionCode==7)]
  
  trips$NIS_CP.ZP=rep(NA,length=ntrips)
  trips$NIS_CP.ZP[Atl_rows]=sample(Atl_CP.ZP,size=length(Atl_rows),replace = TRUE)
  trips$NIS_CP.ZP[Pac_rows]=sample(Pac_CP.ZP,size=length(Pac_rows),replace = TRUE)
  
  Atl_CP.phyto=nisdens.phyto$CPNIS[which(nisdens.phyto$ArrivalRegionCode==2 & nisdens.phyto$SourceRegionCode==7)]
  Pac_CP.phyto=nisdens.phyto$CPNIS[which(nisdens.phyto$ArrivalRegionCode==4 & nisdens.phyto$SourceRegionCode==7)]
  
  trips$NIS_CP.phyto=rep(NA,length=ntrips)
  trips$NIS_CP.phyto[Atl_rows]=sample(Atl_CP.phyto,size=length(Atl_rows),replace = TRUE)
  trips$NIS_CP.phyto[Pac_rows]=sample(Pac_CP.phyto,size=length(Pac_rows),replace = TRUE)
}

  ## Determine popdens and exc dens by considering that tank population concentrations may differ from sample concentrations 
  trips$popdens.ZP<-ifelse(trips$sampledens.ZP==0,abs(rnorm(length(trips$sampledens.ZP==0),0,1)),
                           ifelse(trips$sampledens.ZP>0,abs(rnorm(length(trips$sampledens.ZP>0),trips$sampledens.ZP,sqrt(trips$sampledens.ZP))),0))
  
  trips$popdens.phyto<-ifelse(trips$sampledens.phyto==0,abs(rnorm(length(trips$sampledens.phyto==0),0,1)),
                              ifelse(trips$sampledens.phyto>0,abs(rnorm(length(trips$sampledens.phyto>0),trips$sampledens.phyto,sqrt(trips$sampledens.phyto))),0))

  trips$excdens.ZP<-ifelse(trips$excdens.ZP==0,abs(rnorm(length(trips$excdens.ZP==0),0,1)),
                           ifelse(trips$excdens.ZP>0,abs(rnorm(length(trips$excdens.ZP>0),trips$excdens.ZP,sqrt(trips$excdens.ZP))),0))
  
  trips$excdens.phyto<-ifelse(trips$excdens.phyto==0,abs(rnorm(length(trips$excdens.phyto==0),0,1)),
                              ifelse(trips$excdens.phyto>0,abs(rnorm(length(trips$excdens.phyto>0),trips$excdens.phyto,sqrt(trips$excdens.phyto))),0))
  
  ## Assign concentratinos after successful IMO treatment based on sampling from theoretical distribution from empirical data
  trips$imodens.ZP=rpois(length(trips[,1]),fit_pass.ZP$estimate[1])
  trips$imodens.ZP[which(trips$imodens.ZP>trips$popdens.ZP)]=trips$popdens.ZP[which(trips$imodens.ZP>trips$popdens.ZP)]
  trips$imodens.ZP[which(trips$imodens.ZP>=10)]=9.9
  
  trips$imodens.phyto=rpois(length(trips[,1]),fit_pass.phyto$estimate[1])
  trips$imodens.phyto[which(trips$imodens.phyto>trips$popdens.phyto)]=trips$popdens.phyto[which(trips$imodens.phyto>trips$popdens.phyto)]
  trips$imodens.phyto[which(trips$imodens.phyto>=10)]=9.9
  
  ## Determine fail concentration by sampling from real fail data for ZP
  trips$faildens.ZP=rlnorm(length(trips[,1]),fit_fail$estimate[1],fit_fail$estimate[2])
  trips$faildens.phyto=rlnorm(length(trips[,1]),fit_fail$estimate[1],fit_fail$estimate[2])
  ## If generated fail concentratino is greater than pop concuentration, replace with population concentration 
  trips$faildens.ZP[which(trips$faildens.ZP>trips$popdens.ZP)]=trips$popdens.ZP[which(trips$faildens.ZP>trips$popdens.ZP)]
  trips$faildens.phyto[which(trips$faildens.phyto>trips$popdens.phyto)]=trips$popdens.phyto[which(trips$faildens.phyto>trips$popdens.phyto)]
  ## fail threshold was something that was evaluated as a side sensitivity but is no longer relevant
  trips$faildens.ZP[which(trips$faildens.ZP<fail_threshold)]=fail_threshold
  trips$faildens.phyto[which(trips$faildens.phyto<fail_threshold)]=fail_threshold
 
  #Calculate survival probability for species based on environmental distance between source and recipient port (i.e. no exchange)
  logOddsS = model$coefficients[1] +model$coefficients[2]*trips$ENV_DIST
  trips$pEst_reg = exp(logOddsS)/(1+ exp(logOddsS))
  
  #Calculate survival probability for species based on environmental distance between exchange location and recipient port (i.e. mid-ocean exchange)
  logOddsE = model$coefficients[1] +model$coefficients[2]*trips$ENV_DIST_EXC
  trips$pEst_exc = exp(logOddsE)/(1+ exp(logOddsE))
  
  ## Select a row from empirical data to CP in data matches CP in our simulation to use relative abundance information from 
  trips$NIS_ZP_row=NA
  trips$NIS_Phyto_row=NA
  for (i in 1:length(trips[,1])) {
    
    if(length(which(nisdens.ZP$CPNIS==trips$NIS_CP.ZP[i]))>1) {
      trips$NIS_ZP_row[i]=sample(which(nisdens.ZP$CPNIS==trips$NIS_CP.ZP[i]),1)}
    
    if(length(which(nisdens.ZP$CPNIS==trips$NIS_CP.ZP[i]))==1) {
      trips$NIS_ZP_row[i]=which(nisdens.ZP$CPNIS==trips$NIS_CP.ZP[i])}
    
    if(length(which(nisdens.phyto$CPNIS==trips$NIS_CP.phyto[i]))>1) {
      trips$NIS_Phyto_row[i]=sample(which(nisdens.phyto$CPNIS==trips$NIS_CP.phyto[i]),1)}
    
    if(length(which(nisdens.phyto$CPNIS==trips$NIS_CP.phyto[i]))==1) {
      trips$NIS_Phyto_row[i]=which(nisdens.phyto$CPNIS==trips$NIS_CP.phyto[i])}
  }
  trips$NIS_ZP_row=unlist(trips$NIS_ZP_row)
  trips$NIS_Phyto_row=unlist(trips$NIS_Phyto_row)
  
  #Create array with assemblage data for ZP and phyto 
  assemblage.ZP=array(rep(NA, ntrips*110*10), dim=c(ntrips, 110, 11))
  assemblage.phyto=array(rep(NA, ntrips*179*10), dim=c(ntrips, 51, 11))
  dimnames(assemblage.ZP) <- list(1:ntrips, 1:110,c("rel.abund","alpha.draw","T0","Exc","alpha.un","alpha.ex","P.Env.No.Exc","P.Env.Exc","Surv.No.Exc","Surv.Exc","red.abund"))
  dimnames(assemblage.phyto) <- list(1:ntrips, 1:51,c("rel.abund","alpha.draw","T0","Exc","alpha.un","alpha.ex","P.Env.No.Exc","P.Env.Exc","Surv.No.Exc","Surv.Exc","red.abund"))
  
  # Set reduced abundance to 0 as a baseline to be replaced in next part of script 
  assemblage.ZP[,,11]=0
  assemblage.phyto[,,11]=0
  
  for (i in 1:ntrips) {
    ## Reduce trip percentrage reduction if it will put the concentration below the failure threshold 
    if(trips$popdens.ZP[i]*(1-percent.reduction)<fail_threshold) {trip.per.red=1-(fail_threshold/trips$popdens.ZP[i])
    #Otherwise, use the modelled trip reduction 
    } else {trip.per.red=percent.reduction}
    
    ## Put relative abundance data from randomly selected row into appropriate trip row
    assemblage.ZP[i,,1]=unlist(nisdens.ZP[trips$NIS_ZP_row[i],6:115])
    assemblage.phyto[i,,1]=unlist(nisdens.phyto[trips$NIS_Phyto_row[i],4:54])
    
    ## Put reduced abundance into red.abund column 
    a=as.vector(which(assemblage.ZP[i,,1]!=0))
    if(length(a)==1) {assemblage.ZP[i,a,11]=(1-trip.per.red)}
    if(length(a)>1) {assemblage.ZP[i,a,11]=assemblage.ZP[i,a,1]*(1-trip.per.red)}

    a=as.vector(which(assemblage.phyto[i,,1]!=0))
    if(length(a)==1) assemblage.phyto[i,a,11]=(1-trip.per.red)
    if(length(a)>1) {assemblage.phyto[i,a,11]=assemblage.phyto[i,a,1]*(1-trip.per.red)}
    
    ## Put actual concentration in T0 column 
    assemblage.ZP[i,,3]=assemblage.ZP[i,,1]*trips$popdens.ZP[i]
    assemblage.phyto[i,,3]=assemblage.phyto[i,,1]*trips$popdens.phyto[i]

  }
  # Determine how much of assemblage will change based on exchange type 
  mean_change=NA
  N_changed.ZP=NA
  N_changed.phyto=NA
  for (i in 1:ntrips) {
    if (trips$Exc_Type[i]=="ER") {
      mean_change[i]=rnorm(1,mean=97.9,sd=2/3)}
    if (trips$Exc_Type[i]=="FT") {
      mean_change[i]=rnorm(1,mean=70.1,sd=29.9/3)}
  }
    mean_change[mean_change>100]=100
    mean_change[mean_change<0]=0
    assemblage.ZP[,,4]=0
    assemblage.phyto[,,4]=0
    
    # Round to nearest whole species and determine which species will be changed 
    for (i in 1:ntrips) {
          N_changed.ZP[i]=round(mean_change[i]*0.01*trips$NIS_CP.ZP[i],digits=0)
          # print(i)
          # print(trips$NIS_CP.ZP[i])
          # print(N_changed.ZP[i])
          N_changed.phyto[i]=round(mean_change[i]*0.01*trips$NIS_CP.phyto[i],digits=0)
          # print(length(which(assemblage.ZP[i,,1]>0)))
          # print(N_changed.ZP[i])
          # print(length(which(assemblage.phyto[i,,1]>0)))
          # print(N_changed.phyto[i])
          
          if(length(which(assemblage.ZP[i,,1]>0))>1) {
              assemblage.ZP[i,sample(which(assemblage.ZP[i,,1]>0),N_changed.ZP[i],replace=FALSE),4]=1}
          if(length(which(assemblage.phyto[i,,1]>0))>1) {
              assemblage.phyto[i,sample(which(assemblage.phyto[i,,1]>0),N_changed.phyto[i],replace=FALSE),4]=1}
          if(length(which(assemblage.ZP[i,,1]>0))==1) {
            assemblage.ZP[i,which(assemblage.ZP[i,,1]>0),4]=1}
          if(length(which(assemblage.phyto[i,,1]>0))==1) {
            assemblage.phyto[i,which(assemblage.phyto[i,,1]>0),4]=1}
          
          
      }
    
    #### Assign environmental survival probabilities with and without exchange 
    for (i in 1:ntrips) {
    assemblage.ZP[i,,7]=trips$pEst_reg[i]
    assemblage.phyto[i,,7]=trips$pEst_reg[i]
    assemblage.ZP[i,,8]=trips$pEst_reg[i]
    assemblage.phyto[i,,8]=trips$pEst_reg[i]
    assemblage.ZP[i,which(assemblage.ZP[i,,4]==1),8]=trips$pEst_exc[i]
    assemblage.phyto[i,which(assemblage.phyto[i,,4]==1),8]=trips$pEst_exc[i]
    
    }
  
    ## Determine alpha values for species in tank 
    assemblage.ZP[,,2]<-rbeta(length(assemblage.ZP[,,2]), shape1 = alpha_par1, shape2 = alpha_par2)
    assemblage.phyto[,,2]<-rbeta(length(assemblage.phyto[,,2]), shape1 = alpha_par1, shape2 = alpha_par2)
    
    ## Alter alpha values based on salinity match 
    for (i in 1:ntrips) {
      #DHigher alpha values
      if (trips$sal_dif_unexchanged[i]==1) {
        assemblage.ZP[i,,5]=assemblage.ZP[i,,2]/(2*sal_sens)
        assemblage.phyto[i,,5]=assemblage.phyto[i,,2]/(2*sal_sens)
        }
      if (trips$sal_dif_unexchanged[i]>1) {
        assemblage.ZP[i,,5]=assemblage.ZP[i,,2]/(10*sal_sens)
        assemblage.phyto[i,,5]=assemblage.phyto[i,,2]/(10*sal_sens)
      }
      if (trips$sal_dif_unexchanged[i]==0) {
        assemblage.ZP[i,,5]=assemblage.ZP[i,,2]
        assemblage.phyto[i,,5]=assemblage.phyto[i,,2]
      }
      if (trips$sal_dif_exchanged[i]==1) {
        assemblage.ZP[i,,6]=assemblage.ZP[i,,2]/(2*sal_sens)
        assemblage.phyto[i,,6]=assemblage.phyto[i,,2]/(2*sal_sens)
      }
      if (trips$sal_dif_exchanged[i]>1) {
        assemblage.ZP[i,,6]=assemblage.ZP[i,,2]/(10*sal_sens)
        assemblage.phyto[i,,6]=assemblage.phyto[i,,2]/(10*sal_sens)
      }
      if (trips$sal_dif_exchanged[i]==0) {
        assemblage.ZP[i,,6]=assemblage.ZP[i,,2]
        assemblage.phyto[i,,6]=assemblage.phyto[i,,2]
      }
    }
    
    ## Calculated densities for discharge at arrival
    trips$disc.dens.norm.ZP=trips$popdens.ZP*trips$propNIS.ZP
    
    trips$disc.dens.norm.phyto=trips$popdens.phyto*trips$propNIS.phyto
    
    trips$disc.dens.exc.ZP=trips$excdens.ZP*trips$propNIS.ZP
    
    trips$disc.dens.exc.phyto=trips$excdens.phyto*trips$propNIS.phyto
    
    trips$disc.dens.imo.ZP=trips$imodens.ZP*trips$propNIS.ZP
    
    trips$disc.dens.imo.phyto=trips$imodens.phyto*trips$propNIS.phyto
    
    trips$disc.dens.faildist.ZP=trips$faildens.ZP*trips$propNIS.ZP
    
    trips$disc.dens.faildist.phyto=trips$faildens.phyto*trips$propNIS.phyto
    
    ## Calculate surviving NIS propagules with environment 
    
    assemblage.ZP[,,9]=0
    assemblage.ZP[,,10]=0
    assemblage.phyto[,,9]=0
    assemblage.phyto[,,10]=0
   
    for (i in 1:ntrips){
        random.probs=runif(length(assemblage.ZP[i,,1]),0,1)
        assemblage.ZP[i,which(assemblage.ZP[i,,7]>random.probs),9]=1
        random.probs=runif(length(assemblage.ZP[i,,1]),0,1)
        assemblage.ZP[i,which(assemblage.ZP[i,,8]>random.probs),10]=1
        
        random.probs=runif(length(assemblage.phyto[i,,1]),0,1)
        assemblage.phyto[i,which(assemblage.phyto[i,,7]>random.probs),9]=1
        random.probs=runif(length(assemblage.phyto[i,,1]),0,1)
        assemblage.phyto[i,which(assemblage.phyto[i,,8]>random.probs),10]=1
        
        assemblage.ZP[i,which(assemblage.ZP[i,,4]==0),10]=assemblage.ZP[i,which(assemblage.ZP[i,,4]==0),9]
        assemblage.phyto[i,which(assemblage.phyto[i,,4]==0),10]=assemblage.phyto[i,which(assemblage.phyto[i,,4]==0),9]
  
        trips$disc.dens.partial.ZP.noE[i]=trips$popdens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,11])
        trips$disc.dens.partial.ZP.exc[i]=trips$excdens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,11])
        trips$disc.dens.partial.phyto.noE[i]=trips$popdens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,11])
        trips$disc.dens.partial.phyto.exc[i]=trips$excdens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,11])
        
        trips$surv.dens.norm.ZP.noE[i]=trips$popdens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,1]*assemblage.ZP[i,,9])
        trips$surv.dens.norm.ZP.exc[i]=trips$excdens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,1]*assemblage.ZP[i,,10])
        trips$surv.dens.imo.ZP.noE[i]=trips$imodens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,1]*assemblage.ZP[i,,9])
        trips$surv.dens.imo.ZP.exc[i]=trips$imodens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,1]*assemblage.ZP[i,,10])
        trips$surv.dens.partial.ZP.noE[i]=trips$popdens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,11]*assemblage.ZP[i,,9])
        trips$surv.dens.partial.ZP.exc[i]=trips$excdens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,11]*assemblage.ZP[i,,10])
        trips$surv.dens.faildist.ZP.noE[i]=trips$faildens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,1]*assemblage.ZP[i,,9])
        trips$surv.dens.faildist.ZP.exc[i]=trips$faildens.ZP[i]*trips$propNIS.ZP[i]*sum(assemblage.ZP[i,,1]*assemblage.ZP[i,,10])
        
        trips$surv.dens.norm.phyto.noE[i]=trips$popdens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,1]*assemblage.phyto[i,,9])
        trips$surv.dens.norm.phyto.exc[i]=trips$excdens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,1]*assemblage.phyto[i,,10])
        trips$surv.dens.imo.phyto.noE[i]=trips$imodens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,1]*assemblage.phyto[i,,9])
        trips$surv.dens.imo.phyto.exc[i]=trips$imodens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,1]*assemblage.phyto[i,,10])
        trips$surv.dens.partial.phyto.noE[i]=trips$popdens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,11]*assemblage.phyto[i,,9])
        trips$surv.dens.partial.phyto.exc[i]=trips$excdens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,11]*assemblage.phyto[i,,10])
        trips$surv.dens.faildist.phyto.noE[i]=trips$faildens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,1]*assemblage.phyto[i,,9])
        trips$surv.dens.faildist.phyto.exc[i]=trips$faildens.phyto[i]*trips$propNIS.phyto[i]*sum(assemblage.phyto[i,,1]*assemblage.phyto[i,,10])

        
        trips$est.pop.noEx.wV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,5],
                                            N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,9]*trips$popdens.ZP[i]*trips$DisVol[i])))
        
        trips$est.pop.noEx.noV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,5],
                                               N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,9]*trips$popdens.ZP[i])))
        
        trips$est.pop.Ex.wV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,6],
                                               N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,10]*trips$excdens.ZP[i]*trips$DisVol[i])))
        
        trips$est.pop.Ex.noV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,6],
                                               N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,10]*trips$excdens.ZP[i])))
        
        trips$est.imo.noEx.wV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,5],
                                                    N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,9]*trips$imodens.ZP[i]*trips$DisVol[i])))
        
        trips$est.imo.noEx.noV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,5],
                                                     N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,9]*trips$imodens.ZP[i])))
        
        trips$est.imo.Ex.wV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,6],
                                                  N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,10]*trips$imodens.ZP[i]*trips$DisVol[i])))
        
        trips$est.imo.Ex.noV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,6],
                                                   N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,10]*trips$imodens.ZP[i])))
        
        trips$est.partial.noEx.wV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,5],
                                                    N=(assemblage.ZP[i,,11]*assemblage.ZP[i,,9]*trips$popdens.ZP[i]*trips$DisVol[i])))
        
        trips$est.partial.noEx.noV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,5],
                                                     N=(assemblage.ZP[i,,11]*assemblage.ZP[i,,9]*trips$popdens.ZP[i])))
        
        trips$est.partial.Ex.wV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,6],
                                                  N=(assemblage.ZP[i,,11]*assemblage.ZP[i,,10]*trips$excdens.ZP[i]*trips$DisVol[i])))
        
        trips$est.partial.Ex.noV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,6],
                                                   N=(assemblage.ZP[i,,11]*assemblage.ZP[i,,10]*trips$excdens.ZP[i])))
        
        trips$est.faildist.noEx.wV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,5],
                                                        N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,9]*trips$faildens.ZP[i]*trips$DisVol[i])))
        
        trips$est.faildist.noEx.noV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,5],
                                                         N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,9]*trips$faildens.ZP[i])))
        
        trips$est.faildist.Ex.wV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,6],
                                                      N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,10]*trips$faildens.ZP[i]*trips$DisVol[i])))
        
        trips$est.faildist.Ex.noV.ZP[i]=establish(estab(alpha=assemblage.ZP[i,,6],
                                                       N=(assemblage.ZP[i,,1]*assemblage.ZP[i,,10]*trips$faildens.ZP[i])))
        
        trips$est.pop.noEx.wV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,5],
                                                    N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,9]*trips$popdens.phyto[i]*trips$DisVol[i]*1000000)))
        
        trips$est.pop.noEx.noV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,5],
                                                     N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,9]*trips$popdens.phyto[i])))
        
        trips$est.pop.Ex.wV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,6],
                                                  N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,10]*trips$excdens.phyto[i]*trips$DisVol[i]*1000000)))
        
        trips$est.pop.Ex.noV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,6],
                                                   N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,10]*trips$excdens.phyto[i])))
        
        trips$est.imo.noEx.wV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,5],
                                                    N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,9]*trips$imodens.phyto[i]*trips$DisVol[i]*1000000)))
        
        trips$est.imo.noEx.noV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,5],
                                                     N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,9]*trips$imodens.phyto[i])))
        
        trips$est.imo.Ex.wV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,6],
                                                  N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,10]*trips$imodens.phyto[i]*trips$DisVol[i]*1000000)))
        
        trips$est.imo.Ex.noV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,6],
                                                   N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,10]*trips$imodens.phyto[i])))
        
        trips$est.partial.noEx.wV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,5],
                                                       N=(assemblage.phyto[i,,11]*assemblage.phyto[i,,9]*trips$popdens.phyto[i]*trips$DisVol[i]*1000000)))
        
        trips$est.partial.noEx.noV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,5],
                                                        N=(assemblage.phyto[i,,11]*assemblage.phyto[i,,9]*trips$popdens.phyto[i])))
        
        trips$est.partial.Ex.wV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,6],
                                                     N=(assemblage.phyto[i,,11]*assemblage.phyto[i,,10]*trips$excdens.phyto[i]*trips$DisVol[i]*1000000)))
        
        trips$est.partial.Ex.noV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,6],
                                                      N=(assemblage.phyto[i,,11]*assemblage.phyto[i,,10]*trips$excdens.phyto[i])))
        
        trips$est.faildist.noEx.wV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,5],
                                                           N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,9]*trips$faildens.phyto[i]*trips$DisVol[i]*1000000)))
        
        trips$est.faildist.noEx.noV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,5],
                                                            N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,9]*trips$faildens.phyto[i])))
        
        trips$est.faildist.Ex.wV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,6],
                                                         N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,10]*trips$faildens.phyto[i]*trips$DisVol[i]*1000000)))
        
        trips$est.faildist.Ex.noV.phyto[i]=establish(estab(alpha=assemblage.phyto[i,,6],
                                                          N=(assemblage.phyto[i,,1]*assemblage.phyto[i,,10]*trips$faildens.phyto[i])))
    }
  
    write.csv(trips,file=paste("trips_p",p,"_t",t,"_",current_date,".csv",sep=""))
    print(p)
    print(t)
  } #t
    
  } #p

} #s

