#RMJOC II Statistics Metrics and Comparison
#Volume Frequency
#March 2016
# require(e1071)
# require(ggplot2)
# require(hydroTSM)
# require(lfstat)
# require(lubridate)
# require(plyr)
# require(reshape)
# require(scales)
# require(stats)
# library(gtable)
# library(grid)
#########################################################################################################################################################
#Log-Pearson Type III Analysis
#k=1 (base 10) and k=0.4343 (base e)
#Parameters:
# Location(y0)
# Scale(alpha)
# Shape(beta)
#alpha>0 f(x) positively skewed, has lower bound and varies in the range exp(y0)<=x<inf
#alpha<0 f(x) positively or negatively skewed depending on alpha and beta, has upper bound and varies in the range 0<x<=exp(y0)
#rth moment always exists for alpha<0; exists only if r<k/alpha
#If X is log-Pearson type III then Y=log(X) is gamma-3 distributed with parameters alpha, beta, & y0
#NOTE: Data should only be the flow values of interest
###############################################################################################################################
# #STAND ALONE SETUP
# locationName="ALBO" #this will need to be inset in a loop for all the locations desired
# dataPath=paste0(getwd(),"/Data/")
# modelName=data.frame(read.csv(paste0(dataPath,"modelNames.csv"),header=TRUE))
# rcp=list("45","85")
##############################################################################################################################
#READ frequency factor table
#source is http://streamflow.engr.oregonstate.edu/analysis/floodfreq/skew.htm
Kt=read.table(paste0(mainPath,"Rscripts/LPIII_freqFact_main.csv"),header=FALSE,sep=",",skip=7)
#rename columns
colnames(Kt)=c("Cs","yr1.0101","yr2","yr5","yr10","yr20","yr50","yr100","yr200","yr500")
#head(Kt)
########################################################################################################################################################
#Baseline data
#Get the baseline data file but this is only needed if the calling this script as stand-alone and not through the 
#the main script named 'statsAnalysis_main.R'
# baseData=data.frame(read.table(paste0(dataPath,"baseline.csv"),header=TRUE,sep=","))
# #Change dates and add informaiton useful for sorting
# baseData$Date=as.POSIXct(x=baseData$Date,tz="GMT",format="%m/%d/%Y")
# baseData$WY=water_year(baseData$Date,origin="usgs",as.POSIX=TRUE)
# baseData$Day=day(baseData$Date)
# baseData$Mon=month(baseData$Date)
# subIndex=data.frame(read.table(paste0(dataPath,"julianDayLookup.csv"),header=TRUE,sep=","))
# baseData=join(newData,subIndex,by=c("Mon","Day"))
baseData$SeasonalFlow=with(baseData,ifelse(baseData$Mon>=4 & baseData$Mon<=7,baseData$Flow,0))
#Create Daily flow data frame for baseline
baseVolAnn=data.frame(baseData$WY,baseData$Flow)
baseVolAnn=ddply(baseData,.(WY),function(z)sum(z$Flow)*(86400/43560)/1000)
colnames(baseVolAnn)=c("WY","AnnualVol_kaf")
baseVolAJ=data.frame(baseData$WY,baseData$Flow)
baseVolAJ=ddply(baseData,.(WY),function(z)sum(z$SeasonalFlow)*(86400/43560)/1000)
colnames(baseVolAJ)=c("WY","AJ_Vol_kaf")
#-----------------------------------------------------------------------------------------------------------------------
#Determine mean, SD and skewness coefficient for Baseline LPIII analysis
#Compute the log10 values of the seaonal volume data
baseVolAJ$y=log10(baseVolAJ$AJ_Vol_kaf)
#Estimation of sample stats
mu=mean(baseVolAJ$y)
sumSq=ddply(baseVolAJ,.(WY),function(z)(z$y-mu)^2)
SD=sqrt(1/(nrow(baseVolAJ)-1)*sum(sumSq[,2]))
sumCube=ddply(baseVolAJ,.(WY),function(z)(z$y-mu)^3)
gammaEst=nrow(baseVolAJ)*sum(sumCube[,2])/((nrow(baseVolAJ)-1)*(nrow(baseVolAJ)-2)*SD^3)
#Estimation of weighted 
gammaG=0.0 #this shoudld be adjusted for regional skew but put here for testing purposes
MSE_gammaG=0.302
#Estimate of skew
if (gammaEst<=0.9){
  A=-0.33+0.08*abs(gammaEst)
}else{
  A=-0.52+0.30*abs(gammaEst)
}
if (gammaEst<=1.5){
  B=0.94-0.26*abs(gammaEst)
}else{
  B=0.55
}
MSE_gamma=10^(A-B*log10(nrow(baseVolAJ)/10))
#Weighted skew for parameter estimation  
gammaW=(gammaEst*MSE_gammaG+gammaG*MSE_gamma)/(MSE_gammaG+MSE_gamma)
#-----------------------------------------------------------------------------------------------------------------------
#LPIII Quantile estimates for Baseline
#Find frequency factors for weighted skew and various return peiods
#Loop through all the return periods 1-yr through 500-yr
eventFlowBase=list()
for (event in 2:ncol(Kt)){
  x=Kt[,1]
  y=Kt[,event]
  Cs_interp=approx(x,y,gammaW,ties=mean)
  eventFlowBase[event-1]=mu+SD*Cs_interp[[2]]
  #eventFlowFinal=rbind(append(eventFlowFinal,eventFlow[event-1]))
}
eventFlowBase=do.call(rbind.data.frame,eventFlowBase)
eventFlowBase[,2]=as.numeric(c("1.01","2","5","10","20","50","100","200","500"))
colnames(eventFlowBase)=c("Volume","ReturnPeriod")
########################################################################################################################################################
########################################################################################################################################################
#Climate Change datasets
# centerYr=2030 #can be input directly but will generally pull from main script
# UB=upperBound
# LB=lowerBound
# periodCenter=as.POSIXct(paste0(centerYr,"-10-01")) #The user would change the year only to center in another year
# #for (i in 1:3){
#   i=1 #testing only
#   fName=paste0(modelName[i,],"_RCP",rcp[2],"_BCSD_VIC_P1-ALBO-biascorrected_streamflow-provisional_0.1.csv")
#   tempData=data.frame(read.table(paste0(dataPath,fName),header=TRUE,sep=","))
#   colnames(tempData)=c("Date","Flow")
#   #Subset the data to remove null values
#   newData=subset(tempData,Flow!=-9999,select=c(Date,Flow))
#   #Change dates and add informaiton useful for sorting
#   newData$Date=as.POSIXct(newData$Date,tz="GMT")
#   newData$WY=water_year(newData$Date,origin="usgs",as.POSIX=TRUE)
#   newData$Day=day(newData$Date)
#   newData$Mon=month(newData$Date)
#   subIndex=data.frame(read.table(paste0(dataPath,"julianDayLookup.csv"),header=TRUE,sep=","))
#   newData=join(newData,subIndex,by=c("Mon","Day"))
  newData$SeasonalFlow=with(newData,ifelse(newData$Mon>=4 & newData$Mon<=7,newData$Flow,0))
  #Subset data again for period of interest
  #This may need to be changed depending how the visualization is desired
  nd_vf=subset(newData,Date>=as.POSIXct(paste0(centerYr-LB,"-09-30")) & Date<=as.POSIXct(paste0(centerYr+UB,"-09-30")),select=c(Date,Flow,SeasonalFlow,WY,Day,Mon,JulianWY))
  #Create Daily volume data frame for cc dataset
  #"newData" isn't used since this includes all of the years in the CC dataset
  ccVolAnn=data.frame(nd_vf$WY,nd_vf$Flow)
  ccVolAnn=ddply(nd_vf,.(WY),function(z)sum(z$Flow)*(86400/43560)/1000)
  colnames(ccVolAnn)=c("WY","AnnualVol_kaf")
  ccVolAJ=data.frame(nd_vf$WY,nd_vf$SeasonalFlow)
  ccVolAJ=ddply(nd_vf,.(WY),function(z)sum(z$SeasonalFlow)*(86400/43560)/1000)
  colnames(ccVolAJ)=c("WY","AJ_Vol_kaf")
  #------------------------------------------------------------------------------------------------------------------------------------------
  #Determine mean, SD and skewness coefficient for Baseline LPIII analysis
  #Compute the log10 values of the flow data
  ccVolAJ$y=log10(ccVolAJ$AJ_Vol_kaf)
  #Estimation of sample stats
  mucc=mean(ccVolAJ$y)
  sumSqcc=ddply(ccVolAJ,.(WY),function(z)(z$y-mucc)^2)
  SDcc=sqrt(1/(nrow(ccVolAJ)-1)*sum(sumSqcc[,2]))
  sumCubecc=ddply(ccVolAJ,.(WY),function(z)(z$y-mucc)^3)
  gammaEstcc=nrow(ccVolAJ)*sum(sumCubecc[,2])/((nrow(ccVolAJ)-1)*(nrow(ccVolAJ)-2)*SDcc^3)
  #gammaEst #testing only
  #Estimation of weighted 
  gammaGcc=0.0 #this shoudld be adjusted for regional skew but put here for testing purposes
  MSE_gammaGcc=0.302
  #Estimate of skew
  if (gammaEstcc<=0.9){
    A=-0.33+0.08*abs(gammaEstcc)
  }else{
    A=-0.52+0.30*abs(gammaEstcc)
  }
  if (gammaEstcc<=1.5){
    B=0.94-0.26*abs(gammaEstcc)
  }else{
    B=0.55
  }
  MSE_gammacc=10^(A-B*log10(nrow(ccVolAJ)/10))
  #Weighted skew for parameter estimation  
  gammaWcc=(gammaEstcc*MSE_gammaGcc+gammaGcc*MSE_gammacc)/(MSE_gammaGcc+MSE_gammacc)
  #-----------------------------------------------------------------------------------------------------------------------
  #LPIII Quantile estimates for Baseline
  #Find frequency factors for weighted skew and various return peiods
  #Loop through all the return periods 1-yr through 200-yr
  eventFlowCC=list()
  for (event in 2:ncol(Kt)){
    x=Kt[,1]
    y=Kt[,event]
    Cs_interpCC=approx(x,y,gammaWcc,ties=mean)
    eventFlowCC[event-1]=mu+SD*Cs_interpCC[[2]]
    #eventFlowFinal=rbind(append(eventFlowFinal,eventFlow[event-1]))
  }
  eventFlowCC=do.call(rbind.data.frame,eventFlowCC)
  eventFlowCC[,2]=as.numeric(c("1.01","2","5","10","25","50","100","200","500"))
  colnames(eventFlowCC)=c("Flow","ReturnPeriod")
  #Combine baseline and CC datasets and put back in to "real" space for flows
  finalEventFlowLog=merge(eventFlowBase,eventFlowCC,by="ReturnPeriod") #all values still in log space
  colnames(finalEventFlowLog)=c("ReturnPeriod","Baseline","CC")
  finalEventFlow=data.frame(finalEventFlowLog$ReturnPeriod)
  colnames(finalEventFlow)=c("ReturnPeriod")
  finalEventFlow$Baseline=10^finalEventFlowLog$Baseline
  finalEventFlow$CC=10^finalEventFlowLog$CC
  finalEventFlow$PercentFrequency=1/finalEventFlow$ReturnPeriod
  colnames(finalEventFlow)=c("ReturnPeriod",baselineName,paste0(modelName[name,]),"PercentFrequency")
  #####################################################################################################
  #PLOTTING
  ####################################################################################################
  #Functions
  reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
  }
  #----------------------------------------------------------------------------------------------------
  #manipulate dates for plotting
  df_vf=melt(data=finalEventFlow,id=c("ReturnPeriod","PercentFrequency"))
  colnames(df_vf)=c("Return_Period","Percent_Frequency","Dataset","Flow")
  df_vf$Location=locationName
  df_vf_stacked=rbind(df_vf_stacked,df_vf) #create data frame for stacked plots
  colnames(df_vf_stacked)=c("Return_Period","Percent_Frequency","Dataset","Flow","Location")
  if(indPlotsSwitch==1){
    p3=ggplot(df_vf,aes(x=Return_Period,y=Flow/1000),log10="y")+geom_line(aes(color=Dataset,linetype=Dataset),cex=1)+
      scale_x_log10(breaks=c(1,10,50,100,200,500),labels=c(1,10,50,100,200,500))+#scale_y_log10()+
      xlab("Return Period [yrs]")+ylab("Volume [maf]")+
      ggtitle(label=paste0("April-July Volume Frequency Analysis for ",locationName,"\nDaily Flows (WY",centerYr-LB+1,"-",centerYr+UB,") \n",
      modelName[name,],"_rcp",rcp[emissionScen],"_",bcProcess,"_",hydroModelName,"_",parameterSetName," (",nameTag,"_v",version,")"))+
      theme(plot.title=element_text(size=11, face="bold",hjust=0.5),legend.position="none")
    
    p4=ggplot(df_vf,aes(x=Percent_Frequency*100,y=Flow/1000,log10="y"))+geom_line(aes(color=Dataset,linetype=Dataset),cex=1)+
      scale_x_continuous(trans=reverselog_trans(10),breaks=c(100,10,2,1,.5,.2),labels=c(100,10,2,1,0.5,0.2))+
      xlab("Percent Chance Exceedance [%]")+ylab("Volume [maf]")
    
    #grid.arrange(p3,p4,newpage = TRUE)
    # png(paste0(dataPath,"Volume Frequency/",bcProcess,"/",nameTag,"_v",version,"/volFreq_",modelName[name,],"_rcp",rcp[emissionScen],"_",bcProcess,"_",hydroModelName,"_",parameterSetName,"_",
    #            locationName,"_",centerYr,"_",baselineName,".png"), width=6,height=5,units="in",res=200)
    
    g34=arrangeGrob(p3, p4+theme(legend.position = "bottom"),nrow=2) 
    ggsave(paste0("volFreq_",modelName[name,],"_rcp",rcp[emissionScen],"_",bcProcess,"_",hydroModelName,"_",parameterSetName,"_",locationName,"_",centerYr,"_",
                  baselineName,".png"),plot=g34,path=paste0(dataPath,"Volume Frequency/",bcProcess,"/",nameTag,"_v",version),
                  width=6,height=4,units="in")
    #dev.off()
  }
  ######################################################################
  #Save Data
  if (csvswitch==1){
    write.csv(df,file=paste0(dataPath,"/Volume Frequency/volFreq_",modelName[name,],"_rcp",rcp[emissionScen],"_",bcProcess,"_",hydroModelName,"_",parameterSetName,"_",locationName,"_",centerYr,".csv"))
  }
#}
