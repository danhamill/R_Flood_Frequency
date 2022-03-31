##import peak flow data
Qpeaknenmrbk <- read_delim("~/OHWM Research/Rstats/HydroData/Qpeaknenmrbk.txt", "\t", 
                           escape_double = FALSE, na = "NA", trim_ws = TRUE)

pkyear<-as.Date(Qpeaknenmrbk$peak_tm, "%Y-%m-%d")

new.dat<-data.frame(pkyear, Qpeaknenmrbk$peak_cd, Qpeaknenmrbk$gage_ht_cd)
names(new.dat)<-c("pkyear","peakQ", "peakstage")
View(new.dat)

#Using package lmomco to do flood frequency analysis
input_data<-new.dat$peakQ

#choosing a distribution
# choose a distribution 
# (see ?dist.list for available distributions)
# log Pearson 3 is not one of them but is a log-transformed equivalent of pe3
# note: this script recognizes 'lp3' to stand for log Pearson Type 3
dist<-"gev"
fa <- FrequencyAnalysis(series=input_data, distribution=dist) 

# estimate 90% confidence intervals
ci <- BootstrapCI(series=input_data,   # flow data
                  distribution=dist,   # distribution
                  n.resamples = 2.5E4, # number of re-samples to conduct
                  ci = 0.90)           # confidence interval level

# generate frequency plot
frequencyPlot(series=input_data, ci$ci)+labs(title="Baker River USGS Gage 01076000, Generalized Extreme Value Distribution")


##Redo analysis using log pearson type 3 (lp3)
dist<-"lp3"
fa <- FrequencyAnalysis(series=input_data, distribution=dist)

# estimate 90% confidence intervals
ci <- BootstrapCI(series=input_data,   # flow data
                  distribution=dist,   # distribution
                  n.resamples = 2.5E4, # number of re-samples to conduct
                  ci = 0.90)           # confidence interval level

# generate frequency plot
frequencyPlot(series=input_data, ci$ci)
title("Baker River USGS Gage 01076000")