# In this file we use the Kubios Data imported

rm(list = ls()) # Clear environment
cat("\014") # Clear console
dev.off() # Clear plot window

library(yarrr)
library(lme4)
library(emmeans)
library(pander)

library(reshape)
library(lme4)
library(lmerTest)
library(pander)

library(ggpubr)
library(car)

library(arrow)
library(tibble)
library(reshape2)

nAGQ = 0 # Set to 1 for eventual analysis
plotPrefix = "Plots/"

# Set and Get directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location

BASE_PATH = "Z:\\shares\\ghep_lab\\2020_VanhollebekeKappen_EEGStressMatrices\\Data\\Processed\\Physiological (ECG, EDA, RSP)\\Kubios\\"
# DATA_PATH = paste0(BASE_PATH,"aligned_data")
# FEATURES_PATH = paste0(DATA_PATH,"/features_gemaps/")

mappeddata = read.csv("Z:\\shares\\ghep_lab\\2020_VanhollebekeKappen_EEGStressMatrices\\Data\\Raw\\Physiological (ECG, EDA, RSP)\\mapped_data.csv")

files <- Sys.glob(paste0(BASE_PATH, "\\*hrv.csv")) #Get all CSV files in Kubios directory

HRVData = NULL
for (p in files) {
  print(p)
  elements <- read.csv(p) # Read csv
  index = which(elements == "  RMSSD (ms):                 ", arr.ind = FALSE) # Find index and +1 +3 +5
  
  # Get all relevant values with indexes from csv
  baselineRMSSD = elements[index +1, 1]
  controlRMSSD = elements[index +3, 1]
  stressRMSSD = elements[index +5, 1]
  
  index = which(elements == "  HF (n.u.):                  ", arr.ind = FALSE) # Find index and +2 +4 +6 (Becuase AR ipv FFT)
  baselineHF = elements[index +2, 1]
  controlHF = elements[index +4, 1]
  stressHF = elements[index +6, 1]
  
  # Get edf_Id from filname
  strindex = unlist(gregexpr(pattern ='EEG_',p))
  pptnum_old = substr(p,strindex + 4,strindex + 5)
  
  if (length(which(mappeddata$edf_marker_id == strtoi(pptnum_old))) > 0) { # Check if this is present in mapping and if so: get the new participant number
    pptnum_new = mappeddata$Participant[which(mappeddata$edf_marker_id == strtoi(pptnum_old))]
  } else {
    pptnum_new = NaN
  }
  
  
  HRVData = rbind(HRVData, data.frame(pptnum_old,pptnum_new,baselineRMSSD,controlRMSSD,stressRMSSD, baselineHF, controlHF, stressHF))
  if (nrow(HRVData) > 100) {
    break
  }
  # break
}


HRVData <- melt(HRVData, id.vars = c("pptnum_old", "pptnum_new")) # Get it to long format

write.csv(HRVData, file = paste0(BASE_PATH,'KubiosHRVData.csv'))

################### Get Jonas DATA ###################

HRVData <- read.csv(paste0(BASE_PATH,'KubiosHRVData.csv'))
jonasData <-
  as.data.frame(read_parquet("user_data_per_block.parquet"))

jonasData <- jonasData[ , -which(names(jonasData) %in% c("start_ts","stop_ts"))]

jonasData$block <- as.factor(jonasData$block)
jonasData$ptcpt_id <- as.factor(jonasData$ptcpt_id)

fullData = NULL
for (p in 1:nrow(jonasData)) {
  pptnum = jonasData$ptcpt_id[p]
  moment = jonasData$block[p]
  RMSSDJonas = jonasData$RMSSD[p]
  
  index = which(HRVData$pptnum_new == pptnum & HRVData$variable == paste0(moment,'RMSSD')) # Get right row in HRVData for correct pptnum and block
  RMSSDKubios = HRVData$value[index]
  
  fullData = rbind(fullData, data.frame(pptnum,moment,RMSSDJonas, RMSSDKubios))
    if (nrow(fullData) > 400) { # Edit this when Sofie is finished with anaysis
    break
  }
}


########### CORR PLOT ###########
library("ggpubr")
ggscatter(fullData, x = "RMSSDJonas", y = "RMSSDKubios", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Jonas", ylab = "Kubios")


########### VIZ ####################
# fullData$RMSSDKubios <- as.factor(fullData$RMSSDKubios)

formulas = c('RMSSDKubios ~ moment', 'RMSSDJonas ~ moment')

plotTitles = c('HRV RMSSD Kubios Partial', 'HRV RMSSD Jonas Partial')

data <- fullData
for(i in 2) {
  data <- data[complete.cases(get(names(data)[i + 2], data)), ] # get only complete cases for this specific variable
  # data <- data[is.na(get(names(data)[i + 2], data)), i + 2] <- 0 # First, Turn NA into Zero
  data <- data[get(names(data)[i + 2], data) != 0, ] # Now remove all rows with zeroes
  
  formula <- paste0(formulas[i], ' + (1|pptnum)')
  # Model
  d0.1 <- lmer(formula,data=data)
  d0.2 <- glmer(formula,data=data, family = gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.3 <- glmer(formula,data=data, family = gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
  d0.4 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.5 <- glmer(formula,data=data, family = Gamma(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.6 <- glmer(formula,data=data, family = Gamma(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
  d0.7 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.8 <- glmer(formula,data=data, family = inverse.gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.9 <- glmer(formula,data=data, family = inverse.gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
  modelNames = c(d0.1,d0.2,d0.3,d0.4,d0.5,d0.6,d0.7,d0.8,d0.9)
  
  # Model Selection
  tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3), AIC(d0.4), AIC(d0.5), AIC(d0.6), AIC(d0.7), AIC(d0.8), AIC(d0.9))
  chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC
  
  Anova(chosenModel[[1]]) # Run Anova, double square brackets because of list properties
  print("Stats control vs stress:")
  print(emmeans(chosenModel[[1]], pairwise ~ moment , adjust ="fdr", type="response"))
  
  # Plotting
  dpi=600    #pixels per square inch
  # jpeg(paste0(plotPrefix, "Figure", "_", plotTitles[i], ".jpeg"), width=8*dpi, height=4*dpi, res=dpi)
  par(mfcol = c(1, 1))
  plot <- pirateplot(
    formula = formulas[i],
    data = data,
    theme = 1,
    pal = "info",
    main = plotTitles[i],
    bean.f.o = .6, # Bean fill
    point.o = .3,  # Points
    inf.f.o = .7,  # Inference fill
    inf.b.o = .8,  # Inference border
    avg.line.o = 1,  # Average line
    # bar.f.o = .5, # Bar
    inf.f.col = "white",  # Inf fill col
    inf.b.col = "black",  # Inf border col
    avg.line.col = "black",  # avg line col
    bar.f.col = gray(.8),  # bar filling color
    point.pch = 21,
    point.bg = "white",
    point.col = "black",
    point.cex = .7,
    
    xlab = "",
    ylim = c(0, 160)
  )
}

########### VIZ 2 : Full DATA ####################

formulas = c('value ~ variable')

plotTitles = c('HRV RMSSD Kubios')

data <- HRVData

data = data[- grep("RMSSD", data$variable),] # Use this to only check for HF values. Be sure to run data <- HRVData before it
data = data[- grep("HF", data$variable),] # Use this to only check for HF values. Be sure to run data <- HRVData before it


for(i in 1) {
  data <- data[complete.cases(get(names(data)[i + 2], data)), ] # get only complete cases for this specific variable
  # data <- data[is.na(get(names(data)[i + 2], data)), i + 2] <- 0 # First, Turn NA into Zero
  data <- data[get(names(data)[i + 2], data) != 0, ] # Now remove all rows with zeroes
  
  formula <- paste0(formulas[i], ' + (1|pptnum_old)')
  # Model
  d0.1 <- lmer(formula,data=data)
  d0.2 <- glmer(formula,data=data, family = gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.3 <- glmer(formula,data=data, family = gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
  d0.4 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.5 <- glmer(formula,data=data, family = Gamma(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.6 <- glmer(formula,data=data, family = Gamma(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
  # d0.7 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.8 <- glmer(formula,data=data, family = inverse.gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.9 <- glmer(formula,data=data, family = inverse.gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
  modelNames = c(d0.1,d0.2,d0.3,d0.4,d0.5,d0.6,d0.8,d0.9)
  
  # Model Selection
  tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3), AIC(d0.4), AIC(d0.5), AIC(d0.6), AIC(d0.8), AIC(d0.9))
  chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC
  
  Anova(chosenModel[[1]]) # Run Anova, double square brackets because of list properties
  print("Stats control vs stress:")
  print(emmeans(chosenModel[[1]], pairwise ~ variable , adjust ="fdr", type="response"))
  
  # Plotting
  dpi=600    #pixels per square inch
  # jpeg(paste0(plotPrefix, "Figure", "_", plotTitles[i], ".jpeg"), width=8*dpi, height=4*dpi, res=dpi)
  par(mfcol = c(1, 1))
  plot <- pirateplot(
    formula = formulas[i],
    data = data,
    theme = 1,
    pal = "info",
    main = plotTitles[i],
    bean.f.o = .6, # Bean fill
    point.o = .3,  # Points
    inf.f.o = .7,  # Inference fill
    inf.b.o = .8,  # Inference border
    avg.line.o = 1,  # Average line
    # bar.f.o = .5, # Bar
    inf.f.col = "white",  # Inf fill col
    inf.b.col = "black",  # Inf border col
    avg.line.col = "black",  # avg line col
    bar.f.col = gray(.8),  # bar filling color
    point.pch = 21,
    point.bg = "white",
    point.col = "black",
    point.cex = .7,
    
    xlab = "",
    ylim = c(0, 160)
  )
}


########## SUBBLOCK ANALYSIS ##################

jonasData <-
  as.data.frame(read_parquet("user_data_per_subblock.parquet"))

jonasData$block <- as.factor(jonasData$block)
jonasData$pptnum <- as.factor(jonasData$ptcpt_id)

formulas = c('HRV ~ block', 'HR ~ block')

plotTitles = c('HRV', 'HR')

data <- jonasData
for(i in 2) {
  data <- data[complete.cases(get(names(data)[i + 2], data)), ] # get only complete cases for this specific variable
  # data <- data[is.na(get(names(data)[i + 2], data)), i + 2] <- 0 # First, Turn NA into Zero
  data <- data[get(names(data)[i + 2], data) != 0, ] # Now remove all rows with zeroes
  
  formula <- paste0(formulas[i], ' + (1|pptnum)')
  # Model
  d0.1 <- lmer(formula,data=data)
  d0.2 <- glmer(formula,data=data, family = gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.3 <- glmer(formula,data=data, family = gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
  # d0.4 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.5 <- glmer(formula,data=data, family = Gamma(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.6 <- glmer(formula,data=data, family = Gamma(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
  d0.7 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.8 <- glmer(formula,data=data, family = inverse.gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  d0.9 <- glmer(formula,data=data, family = inverse.gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
  modelNames = c(d0.1,d0.2,d0.3,d0.5,d0.6,d0.7,d0.8,d0.9)
  
  # Model Selection
  tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3), AIC(d0.5), AIC(d0.6), AIC(d0.7), AIC(d0.8), AIC(d0.9))
  chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC
  
  Anova(chosenModel[[1]]) # Run Anova, double square brackets because of list properties
  print("Stats control vs stress:")
  print(emmeans(chosenModel[[1]], pairwise ~ block , adjust ="fdr", type="response"))
  
  # Plotting
  dpi=600    #pixels per square inch
  # jpeg(paste0(plotPrefix, "Figure", "_", plotTitles[i], ".jpeg"), width=8*dpi, height=4*dpi, res=dpi)
  par(mfcol = c(1, 1))
  plot <- pirateplot(
    formula = formulas[i],
    data = data,
    theme = 1,
    pal = "info",
    main = plotTitles[i],
    bean.f.o = .6, # Bean fill
    point.o = .3,  # Points
    inf.f.o = .7,  # Inference fill
    inf.b.o = .8,  # Inference border
    avg.line.o = 1,  # Average line
    # bar.f.o = .5, # Bar
    inf.f.col = "white",  # Inf fill col
    inf.b.col = "black",  # Inf border col
    avg.line.col = "black",  # avg line col
    bar.f.col = gray(.8),  # bar filling color
    point.pch = 21,
    point.bg = "white",
    point.col = "black",
    point.cex = .7,
    
    xlab = "",
    # ylim = c(0, 160)
  )
}
