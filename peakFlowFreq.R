#RMJOC II Statistics Metrics and Comparison
#Peak Flow Frequency
#Feb 2016
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
###############################################################################################################################
#READ frequency factor table
#source is http://streamflow.engr.oregonstate.edu/analysis/floodfreq/skew.htm
Kt=read.table(paste0(mainPath,"Rscripts/LPIII_freqFact_main.csv"),header=FALSE,sep=",",skip=7)
#rename columns
colnames(Kt)=c("Cs","yr1.0101","yr2","yr5","yr10","yr20","yr50","yr100","yr200","yr500")
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

#Create Daily flow data frame for baseline
basePeak=data.frame(baseData$WY,baseData$Flow)
basePeak=ddply(baseData,.(WY),function(z)max(z$Flow))
colnames(basePeak)=c("WY","Flow")
#clean up data to ensure no negatives; this may need to be evaluated in more detail if the baseline has negative peaks
basePeak$Flow=with(basePeak,ifelse(basePeak$Flow<0,mean(basePeak$Flow),basePeak$Flow))
#-----------------------------------------------------------------------------------------------------------------------
#Determine mean, SD and skewness coefficient for Baseline LPIII analysis
#Compute the log10 values of the flow data
basePeak$y=log10(basePeak$Flow)
#Estimation of sample stats
mu=mean(basePeak$y)
sumSq=ddply(basePeak,.(WY),function(z)(z$y-mu)^2)
SD=sqrt(1/(nrow(basePeak)-1)*sum(sumSq[,2]))
sumCube=ddply(basePeak,.(WY),function(z)(z$y-mu)^3)
gammaEst=nrow(basePeak)*sum(sumCube[,2])/((nrow(basePeak)-1)*(nrow(basePeak)-2)*SD^3)
#gammaEst #testing only
#Estimation of weighted 
gammaG=0.0 #assume no regional skew adjustment #this shoudld be adjusted for regional skew but put here for testing purposes
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
MSE_gamma=10^(A-B*log10(nrow(basePeak)/10))
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
colnames(eventFlowBase)=c("Flow","ReturnPeriod")
########################################################################################################################################################
########################################################################################################################################################
#Climate Change datasets
# yr=centerYr #can be input directly but will generally pull from main script
# UB=upperBound
# LB=lowerBound
# periodCenter=as.POSIXct(paste0(yr,"-10-01")) #The user would change the year only to center in another year
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
  
  #Subset data again for period of interest
  #This may need to be changed depending how the visualization is desired
  nd_pf=subset(newData,Date>=as.POSIXct(paste0(centerYr-LB,"-09-30")) & Date<=as.POSIXct(paste0(centerYr+UB,"-09-30")),select=c(Date,Flow,WY,Day,Mon,JulianWY))
  #Create Daily flow data frame for CC dataset
  ccPeak=data.frame(nd_pf$WY,nd_pf$Flow)
  ccPeak=ddply(nd_pf,.(WY),function(z)max(z$Flow))
  colnames(ccPeak)=c("WY","Flow")
  #------------------------------------------------------------------------------------------------------------------------------------------
  #Determine mean, SD and skewness coefficient for Baseline LPIII analysis
  #Compute the log10 values of the flow data
  ccPeak$y=log10(ccPeak$Flow)
  #Estimation of sample stats
  mucc=mean(ccPeak$y)
  sumSqcc=ddply(ccPeak,.(WY),function(z)(z$y-mucc)^2)
  SDcc=sqrt(1/(nrow(ccPeak)-1)*sum(sumSqcc[,2]))
  sumCubecc=ddply(ccPeak,.(WY),function(z)(z$y-mucc)^3)
  gammaEstcc=nrow(ccPeak)*sum(sumCubecc[,2])/((nrow(ccPeak)-1)*(nrow(ccPeak)-2)*SDcc^3)
  #gammaEst #testing only
  #Estimation of weighted 
  gammaGcc=0.0 #assume not regional skew adjustment #this shoudld be adjusted for regional skew but put here for testing purposes
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
  MSE_gammacc=10^(A-B*log10(nrow(ccPeak)/10))
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
  }
  eventFlowCC_pf=do.call(rbind.data.frame,eventFlowCC)
  eventFlowCC_pf[,2]=as.numeric(c("1.01","2","5","10","25","50","100","200","500"))
  colnames(eventFlowCC_pf)=c("Flow","ReturnPeriod")
  #Combine baseline and CC datasets and put back in to "real" space for flows
  finalEventFlowLog_pf=merge(eventFlowBase,eventFlowCC_pf,by="ReturnPeriod") #all values still in log space
  colnames(finalEventFlowLog_pf)=c("ReturnPeriod","Baseline",paste0(modelName[name,]))
  finalEventFlow_pf=data.frame(finalEventFlowLog_pf$ReturnPeriod)
  colnames(finalEventFlow_pf)=c("ReturnPeriod")
  finalEventFlow_pf$Baseline=10^finalEventFlowLog_pf$Baseline
  finalEventFlow_pf$CC=10^finalEventFlowLog_pf[,3]
  finalEventFlow_pf$PercentFrequency=1/finalEventFlow_pf$ReturnPeriod
  colnames(finalEventFlow_pf)=c("ReturnPeriod",baselineName,paste0(modelName[name,]),"PercentFrequency")
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
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  #----------------------------------------------------------------------------------------------------
  #manipulate dates for plotting
  df_pf=melt(data=finalEventFlow_pf,id=c("ReturnPeriod","PercentFrequency"))
  colnames(df_pf)=c("Return_Period","Percent_Frequency","Dataset","Flow")
  df_pf$Location=locationName
  df_pf_stacked=rbind(df_pf_stacked,df_pf) #create data frame for stacked plots
  colnames(df_pf_stacked)=c("Return_Period","Percent_Frequency","Dataset","Flow","Location")
  if(indPlotsSwitch==1){
    p1=ggplot(df_pf,aes(x=Return_Period,y=Flow/1000),log10="y")+geom_line(aes(color=Dataset,linetype=Dataset),cex=1)+
       scale_x_log10(breaks=c(1,10,50,100,200,500),labels=c(1,10,50,100,200,500))+#scale_y_log10()+
       xlab("Return Period [yrs]")+ylab("Flow [kcfs]")+
       ggtitle(label=paste0("Peak Flow Frequency Analysis for ",locationName,"\nDaily Flows (WY",centerYr-LB+1,"-",centerYr+UB,") \n",
       paste0(modelName[name,],"_rcp",rcp[emissionScen],"_",bcProcess,"_",hydroModelName,"_",parameterSetName," (",nameTag,"_v",version,")")))+
       theme(plot.title=element_text(size=10, face="bold",hjust=0.5),legend.position = "none")
    
    p2=ggplot(df_pf,aes(x=Percent_Frequency*100,y=Flow/1000,log10="y"))+geom_line(aes(color=Dataset,linetype=Dataset),cex=1)+
      scale_x_continuous(trans=reverselog_trans(10),breaks=c(100,10,2,1,.5,.2),labels=c(100,10,2,1,0.5,0.2))+
      xlab("Percent Chance Exceedance [%]")+ylab("Flow [kcfs]")+theme(legend.position = "none")
    #legend=g_legend(p1) 
    #grid.arrange(p1,p2)
    #grid.draw(rbind(ggplotGrob(p1),ggplotGrob(p2),size="last")) #for visualizing in Rstudio
    # png(paste0(dataPath,"/Peak Frequency/",bcProcess,"/",nameTag,"_v",version,"/peakFreq_",modelName[name,],"_rcp",rcp[emissionScen],"_",bcProcess,"_",hydroModelName,"_",parameterSetName,"_",
    #            locationName,"_",centerYr,"_",baselineName,".png"), width=6,height=5,units="in",res=200)
    g12=arrangeGrob(p1, p2+theme(legend.position = "bottom"),nrow=2) 
    ggsave(paste0("peakFreq_",modelName[name,],"_rcp",rcp[emissionScen],"_",bcProcess,"_",hydroModelName,"_",parameterSetName,"_",locationName,"_",centerYr,"_",
           baselineName,".png"),plot=g12,path=paste0(dataPath,"Peak Frequency/",bcProcess,"/",nameTag,"_v",version),
           width=6,height=4,units="in")
    #dev.off()
  }
  ######################################################################
  #Save Data
  if (csvswitch==1){
    write.csv(df_pf,file=paste0(dataPath,"/Peak Frequency/peakfreq_",modelName[name,],"_rcp",rcp[emissionScen],"_",bcProcess,"_",hydroModelName,"_",parameterSetName,"_",locationName,"_",centerYr,".csv"))
  }
#}
