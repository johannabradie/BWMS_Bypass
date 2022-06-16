########## FUNCTIONS #######
library(fitdistrplus)
    
### This scripts fits the theoretical curves based on empirical sample data used to generate total tank concentrations, NIS concentrations
  
    ##Read rtrips data which contains one year of ship transits to Canada
    rtrips<-read.csv("./rtrips_May2020.csv", header=T, quote="\"")
    #remove the two records with ID 3729 and 2550 - these are atl to pac and gl to pac which we are not considering#
    rtrips<-rtrips[!(rtrips$ArrivalRegionCode==4 & rtrips$SourceCode==2),]
    rtrips<-rtrips[!(rtrips$ArrivalRegionCode==3 & rtrips$SourceCode==4),]
    rtrips$DisVol=as.numeric(rtrips$DisVol)
    
    ## Reads in port environmental data
    portdata<-read.csv("portdataSep2021.csv")
    portdata_cols=c((which(colnames(portdata)=="Min.Temp")),which(colnames(portdata)=="Max.Temp"),which(colnames(portdata)=="Annual.Temp"),which(colnames(portdata)=="Salinity"))
    #scale port data
    portdata=cbind(portdata,scale(portdata[,portdata_cols]))
    colnames(portdata)[15:18]=c("Scaled.Min.Temp","Scaled.Max.Temp","Scaled.Annual.Temp","Scaled.Salinity")
    rm(portdata_cols)
    
    ## Merge tables with port data and environmental data 
    rtrips<-merge(x=rtrips,y=portdata[,1:2],by.x="ArrivalPortOnly",by.y="PortName",all.x=T)
    rtrips<-merge(x=rtrips,y=portdata[,1:2],by.x="SourcePortOnly",by.y="PortName",all.x=T)
    colnames(rtrips)[32:33]=c("Arr_PortID","Source_PortID")
    
    # GL_Ports_Class=read.csv("GLports_US_CAN_Classification.csv")
    # rtrips=merge(rtrips,GL_Ports_Class,by.x="ArrivalPort",by.y="Port",all.x=T)
    
    #Read in data with exchange locations and environment in those locations 
    exchange_data <-read.csv("exchange_data_Sep11_2019.csv")
    ## Drop any exchange locations where environmental data was not available nearby
    exchange_data<-subset(exchange_data,exchange_data$Distance<50)
    Salinity=rep(NA,length(exchange_data[,1]))
    Scaled.Min.Temp=rep(NA,length(exchange_data[,1]))
    Scaled.Max.Temp=rep(NA,length(exchange_data[,1]))
    Scaled.Annual.Temp=rep(NA,length(exchange_data[,1]))
    Scaled.Salinity=rep(NA,length(exchange_data[,1]))
    exchange_data<-cbind(exchange_data,Salinity,Scaled.Min.Temp,Scaled.Max.Temp,Scaled.Annual.Temp,Scaled.Salinity)
    rm(Salinity,Scaled.Min.Temp,Scaled.Max.Temp,Scaled.Annual.Temp,Scaled.Salinity)
    # ### Adds scaled temp data and salinity data to data with exchange locations
    
    
    for (i in 1:length(exchange_data[,1])) {
      replacement=intersect(which(portdata$Latitude==exchange_data$closest_lat[i]),which(portdata$Longitude==exchange_data$closest_long[i]))
      if(length(replacement)==0) {print(paste(i,", Error port not found"))}
      exchange_data[i,c(13:17)]=portdata[replacement[1],c(11,15:18)]
    }
    # 
    # ### Adds scaled temperature data to rtrips 
    rtrips<-merge(x=rtrips,y=portdata[,c(1,11,15:18)],by.x="Arr_PortID",by.y="Port.ID",all.x=T)
    names(rtrips)[34:38]=c("Arr.Salinity","Arr.Scaled.Min.Temp.","Arr.Scaled.Max.Temp","Arr.Scaled.Annual.Temp","Arr.Scaled.Salinity")
    rtrips<-merge(x=rtrips,y=portdata[,c(1,11,15:18)],by.x="Source_PortID",by.y="Port.ID",all.x=T)
    names(rtrips)[39:43]=c("Source.Salinity","Source.Scaled.Min.Temp.","Source.Scaled.Max.Temp","Source.Scaled.Annual.Temp","Source.Scaled.Salinity")
    
    ##If you are getting errors for ports not found, check the text fields for the port data and rtrips file.  They must contain '.' instead of a space.
    ## Also make sure port data file has actual port data (by port name) and lat/long point data from WOA.

########################################################################################
## Calculate salinity difference between source/arrival and exchange/arrival
    ## All exchange locations are marine, so exact identity of exchange location is not needed to compute this 

    
    tempA=rtrips$Arr.Salinity
    tempA[tempA<=5]=1  # Freshwater is 1
    tempA[tempA>18]=3  # Marine is 3
    tempA[tempA!=1&tempA!=3]=2  # Brackish is 2
    tempS=rtrips$Source.Salinity
    tempS[tempS<=5]=1
    tempS[tempS>18]=3
    tempS[tempS!=1&tempS!=3]=2
    tempE=rep(3,length=length(tempS))
    
    rtrips$sal_dif_unexchanged=abs(tempS-tempA)
    rtrips$sal_dif_exchanged=abs(tempE-tempA)
    
    rm(i,replacement,tempA,tempS,tempE)



## Distributions for ZP ####

    alldens<-read.table("./ALLDENS.txt", header=T, quote="\"")
    nisdens<-read.table("./ASSEMB.txt", header = T, quote="\"")
    nisdens.ZP<-read.table("./Assemb_ZP_TC_Failed.txt", header = T, quote="\"")
    GLdens<-read.csv("GL_DOM_Zoo_Data.csv")
    colnames(GLdens)[1]="ArrivalRegion"
    
    #Fit tank density distributions for ZP
    x<-round(alldens[alldens$ArrivalRegion==4 & alldens$SourceRegion==7,4]) # Pacific
    PacInt.maxZP=max(x)*1.1
    fit=fitdistr(x, "negative binomial")
    pacsz.ZP<-as.numeric(fit$estimate[1])
    pacmu.ZP<-as.numeric(fit$estimate[2])
    
    x<-round(alldens[alldens$ArrivalRegion==2 & alldens$SourceRegion==7,4]) # Atlantic
    AtlInt.maxZP=max(x)*1.1
    fit=fitdistr(x, "negative binomial")
    atlsz.ZP<-as.numeric(fit$estimate[1])
    atlmu.ZP<-as.numeric(fit$estimate[2])
    
    x<-round(alldens[alldens$ArrivalRegion==3 & alldens$SourceRegion==7,4]) # GL
    GLInt.maxZP=max(x)*1.1
    fit=fitdistr(x, "negative binomial")
    glisz.ZP<-as.numeric(fit$estimate[1])
    glimu.ZP<-as.numeric(fit$estimate[2])
    
    # x<-round(GLdens[GLdens$ArrivalRegion==3 & GLdens$SourceRegion!=7,4]) # GL domestic in all Canada
   
    # # x=x[-which(GLdens$Ship.ID=="L022-01")]
    # fit=fitdistr(x, "lognormal")
    # gldommean.ZP<-as.numeric(fit$estimate[1])
    # gldomsdlog.ZP<-as.numeric(fit$estimate[2])
    # # ad.test(rlnorm(n=100,meanlog=fit$estimate[1],sdlog=as.numeric(fit$estimate[2])),x)
    
    x<-round(alldens[alldens$ArrivalRegion==3& alldens$SourceRegion==3,4]) #GL Domestic with GL data.
    GLDom.maxZP=max(x)*1.1
    # fit=fitdistr(x, "negative binomial")
    glglsz.ZP<-0.4034  #WAS fit by hand by Andrew
    glglmu.ZP<-123550.7
    
  # ##Fit proportion NIS curves
    #to handle beta fitting, remove 0 and 1 for propNIS, replace with 0.00001, 0.99999#
    alldens$PropNIS[alldens$PropNIS==0]<-0.00001
    alldens$PropNIS[alldens$PropNIS==1]<-0.99999
    
    x<-alldens[alldens$ArrivalRegion==4 & alldens$SourceRegion==7,5]
    fit=fitdistr(x, "beta", start = list(shape1 = 1, shape2 = 5))
    pacshape1.ZP<-as.numeric(fit$estimate[1])
    pacshape2.ZP<-as.numeric(fit$estimate[2])
    
    x<-alldens[alldens$ArrivalRegion==2 & alldens$SourceRegion==7,5]
    fit=fitdistr(x, "beta", start = list(shape1 = 1, shape2 = 5))
    atlshape1.ZP<-as.numeric(fit$estimate[1])
    atlshape2.ZP<-as.numeric(fit$estimate[2])
    
    x<-alldens[alldens$ArrivalRegion==3 & alldens$SourceRegion==7,5]
    fit=fitdistr(x, "beta", start = list(shape1 = 1, shape2 = 2))
    glishape1.ZP<-as.numeric(fit$estimate[1])
    glishape2.ZP<-as.numeric(fit$estimate[2])
    
    x<-alldens[alldens$ArrivalRegion==3 & alldens$SourceRegion!=7,5]
    fit=fitdistr(x, "beta", start = list(shape1 = 1, shape2 = 2))
    # ad.test(rbeta(n=25,shape1=fit$estimate[1],shape2=as.numeric(fit$estimate[2])),x)
    gldomshape1.ZP<-as.numeric(fit$estimate[1])
    gldomshape2.ZP<-as.numeric(fit$estimate[2])
  
    treatment_compliance=read.csv('Post_Treatment_Densities_ZP1.csv')
    fit_pass.ZP=fitdistr(treatment_compliance$Density[which(treatment_compliance$Outcome=="Pass")], "poisson")
    fit_fail=fitdistr(treatment_compliance$Density[which(treatment_compliance$Outcome=="Fail")], "lognormal")

### Phyto distributions ### 

    alldens<-read.csv("./ALLDENS_allphyto1_OCM.csv")
    nisdens<-read.csv("./ASSEMB_allphyto_OCM.csv") 
    nisdens.phyto<-read.csv("./ASSEMB_phyto_TCFailed.csv") 
  
  # Fit tank density distributions for phytoplankton
    x<-round(alldens[alldens$ArrivalRegion==4 & alldens$SourceRegion==7,4])
    PacInt.maxphyto=max(x)*1.1
    fit=fitdistr(x, "negative binomial")
    pacsz.phyto<-as.numeric(fit$estimate[1])
    pacmu.phyto<-as.numeric(fit$estimate[2])
    
    x<-round(alldens[alldens$ArrivalRegion==2 & alldens$SourceRegion==7,4])
    AtlInt.maxphyto=max(x)*1.1
    fit=fitdistr(x, "negative binomial")
    atlsz.phyto<-as.numeric(fit$estimate[1])
    atlmu.phyto<-as.numeric(fit$estimate[2])
    
    x<-round(alldens[alldens$ArrivalRegion==3 & alldens$SourceRegion==7,4])
    GLInt.maxphyto=max(x)*1.1
    fit=fitdistr(x, "negative binomial")
    glisz.phyto<-as.numeric(fit$estimate[1])
    glimu.phyto<-as.numeric(fit$estimate[2])
  
    treatment_compliance=read.csv('Post_Treatment_Densities_ZP1.csv')
    fit_fail.phyto=fitdistr(treatment_compliance$Density[which(treatment_compliance$Outcome=="Fail")], "lognormal")
    
    treatment_compliance=read.csv('Post_Treatment_Densities_Phyto1.csv')
    fit_pass.phyto=fitdistr(round(treatment_compliance$Density[which(treatment_compliance$Outcome=="Pass")]), "poisson")
  
  #Fit NIS density distributions for all phytoplankton
    #to handle beta fitting, remove 0 and 1 for propNIS, replace with 0.00001, 0.99999#
    alldens$PropNIS[alldens$PropNIS==0]<-0.00001
    alldens$PropNIS[alldens$PropNIS==1]<-0.99999
    
    x<-alldens[alldens$ArrivalRegion==4 & alldens$SourceRegion==7,5]
    fit=fitdistr(x, "beta", start = list(shape1 = 0.9, shape2 = 4))
    pacshape1.phyto<-as.numeric(fit$estimate[1])
    pacshape2.phyto<-as.numeric(fit$estimate[2])
    
    x<-alldens[alldens$ArrivalRegion==2 & alldens$SourceRegion==7,5]
    fit=fitdistr(x, "beta", start = list(shape1 = 0.9, shape2 = 4))
    atlshape1.phyto<-as.numeric(fit$estimate[1])
    atlshape2.phyto<-as.numeric(fit$estimate[2])
    
    x<-alldens[alldens$ArrivalRegion==3 & alldens$SourceRegion==7,5]
    fit=fitdistr(x, "beta", start = list(shape1 = .9, shape2 = 4))
    glishape1.phyto<-as.numeric(fit$estimate[1])
    glishape2.phyto<-as.numeric(fit$estimate[2])


### Environmental distance model ###

    ### Read in a list of presence distances and absences distances 
    full_P_dist<-read.csv("P_dist.csv")  #Table with presence distances for NIS
    full_A_dist<-read.csv("A_dist.csv")  # Table with background (absence) distances for NIS, computed in different project. See bradie et al 2021
    full_P_dist=full_P_dist[,-1]
    full_A_dist=full_A_dist[,-1]
    
    P_pts=cbind(full_P_dist,rep(1,length(full_P_dist)))
    b=sample(full_A_dist,length(full_P_dist),replace=FALSE)
    A_pts=cbind(b,rep(0,length(full_P_dist)))
    all_data=rbind(P_pts,A_pts)
    colnames(all_data)=list("dist","status")
    all_data=as.data.frame(all_data)
    sample.dat=sample(1:length(all_data[,1]),size=(0.8*length(all_data[,1])),replace=F)
    fitting=as.data.frame(all_data[sample.dat,])
    testing=as.data.frame(all_data[-sample.dat,])
    fitting=as.data.frame(fitting)
    model <- glm(status ~dist,family=binomial(link='logit'),data=fitting)

rm(A_pts,all_data,alldens,fit,fitting,GLdens,nisdens,P_pts,testing,treatment_compliance,b,full_A_dist,full_P_dist,sample.dat,x)
