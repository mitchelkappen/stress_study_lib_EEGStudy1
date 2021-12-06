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

nAGQ = 0 # Set to 1 for eventual analysis

BASEPATH <- "D:/Data/EEG_Study_1/uz_study/features/"
plotPrefix <- paste0(BASEPATH, "figures/")
plotPrefix <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/figures/")

data <-
  as.data.frame(read.csv(paste0(BASEPATH, "audioandresponsedata.csv"))) # This contains everything except for the ones that were too bad -- see Excel drive


data$participantNum <- as.factor(data$participantNum)

# Delete rows with zeroes
# sum(data$F0semitoneFrom27.5Hz_sma3nz_amean == 0)

# Add columns for subblock 1-2-3
data$time[data$fileNum == 0] = 0
data$time[data$fileNum == 1|data$fileNum == 5] = 1
data$time[data$fileNum == 2|data$fileNum == 6] = 2
data$time[data$fileNum == 3|data$fileNum == 7] = 3
data <- add_column(data[,1:ncol(data)-1], time = data$time, .after = "fileNum")

# Add columns for block 0-1-2
data$block[data$fileNum == 0] = 'baseline'
data$block[data$fileNum == 1|data$fileNum == 2|data$fileNum == 3] = 'control'
data$block[data$fileNum == 5|data$fileNum == 6|data$fileNum == 7] = 'stress'
data <- add_column(data[,1:ncol(data)-1], block = data$block, .after = "time")

data$time <- as.factor(data$time)
data$block <- as.factor(data$block)

# Take out all NA's values (Resting blocks)
data = data[is.na(data$block) == FALSE, ] # Take out all NA's

# Add F0 baseline corrected
for(i in unique(data$participantNum)){
  if(sum(data$participantNum == i & data$time == 0) > 0){
    if(sum(data$participantNum == i & data$time == 1) > 0){
      data$F0BaselineCorrected[data$participantNum == i & data$time == 1] = data$F0semitoneFrom27.5Hz_sma3nz_amean[data$participantNum == i & data$time == 1] - data$F0semitoneFrom27.5Hz_sma3nz_amean[data$participantNum == i & data$time == 0]
    }
    if(sum(data$participantNum == i & data$time == 2) > 0){
      data$F0BaselineCorrected[data$participantNum == i & data$time == 2] = data$F0semitoneFrom27.5Hz_sma3nz_amean[data$participantNum == i & data$time == 2] - data$F0semitoneFrom27.5Hz_sma3nz_amean[data$participantNum == i & data$time == 0]
    }
    if(sum(data$participantNum == i & data$time == 3) > 0){
      data$F0BaselineCorrected[data$participantNum == i & data$time == 3] = data$F0semitoneFrom27.5Hz_sma3nz_amean[data$participantNum == i & data$time == 3] - data$F0semitoneFrom27.5Hz_sma3nz_amean[data$participantNum == i & data$time == 0]
      # break
    }
  }
}

# Get relevant data
dataBackup = data # Backup data so we can go back to this whenever
data = data[data$block ==  'control' | data$block == 'stress',] # Get only block 1 and 2 data

data = data[data$HNRdBACF_sma3nz_amean > 0, ] # Kick out all the lines with negative HNR, that means the signal is too noisy to trust
######## STATS AND VIZ LETS GO ########

# formulas = c('F0semitoneFrom27.5Hz_sma3nz_amean ~ time', 'jitterLocal_sma3nz_amean ~ time', 'shimmerLocaldB_sma3nz_amean ~ time', 'HNRdBACF_sma3nz_amean ~ time', 'VoicedSegmentsPerSec ~ time', 'MeanVoicedSegmentLengthSec ~ time')
# formulas = c('F0semitoneFrom27.5Hz_sma3nz_amean ~ block', 'jitterLocal_sma3nz_amean ~ block', 'shimmerLocaldB_sma3nz_amean ~ block', 'HNRdBACF_sma3nz_amean ~ block', 'VoicedSegmentsPerSec ~ block', 'MeanVoicedSegmentLengthSec ~ block', 'valence ~ block', 'arousal ~ block')
# formulas = c('valence ~ time * block', 'arousal ~ time * block')
formulas = c('valence ~ fileNum', 'arousal ~ fileNum')

# plotTitles = c('F0', 'Jitter', 'Shimmer', 'HNR', 'VoicedPerSec', 'MeanVoicedLength', 'Valence', 'Arousal')
plotTitles = c('Valence', 'Arousal')

# for(i in 1:length(formulas)) {
for(i in 1) {
  formula <- paste0(formulas[i], ' + (1|participantNum)')
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
  print("Stats baseline vs control vs stress:")
  # print(emmeans(chosenModel[[1]], pairwise ~ condition , adjust ="fdr", type="response")) # This is the right one for physiological
  print(emmeans(chosenModel[[1]], pairwise ~ block , adjust ="fdr", type="response")) # This is the right one for self-reports
  
  # Plotting
  # dpi=600    #pixels per square inch
  # jpeg(paste0(plotPrefix, "Figure", "_", plotTitles[i], ".jpeg"), width=8*dpi, height=8*dpi, res=dpi)
  par(mfcol = c(1, 1))
  plotAROUSAL <- pirateplot(
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
  )
  # abline(lm(formulas[i], data=data), lwd=4, lty=2, col = "red")
  
  # Secondary points and axis:
  # mtext("fileNum",1,line=1,at=0.2,col="red")
  
  # mtext("Experiment phase",1,line=3,at=0.2,col="blue")
  # axis(side=1, at=c(1:9), line = 3, labels=c('prerest','control1','control2','control3','midrest','stress1','stress2','stress3','postrest' ))
  # axis(side=1, at=c(1:6), line = 3, labels=c('control1','control2','control3','stress1','stress2','stress3' ))
  
  dev.off()
}


############################# CARDIO DATA ######################

IBIdata <-
  as.data.frame(read_parquet("df_tot_merged.parquet"))



data <-
  as.data.frame(read.csv(paste0(BASEPATH, "SAMscompiled.csv")))

# Add columns for subblock 1-2-3
data$time[data$fileNum == 0] = 0
data$time[data$fileNum == 1|data$fileNum == 5] = 1
data$time[data$fileNum == 2|data$fileNum == 6] = 2
data$time[data$fileNum == 3|data$fileNum == 7] = 3
data <- add_column(data[,1:ncol(data)-1], time = data$time, .after = "fileNum")

# Add columns for block 0-1-2
data$block[data$fileNum == 0] = 'baseline'
data$block[data$fileNum == 1|data$fileNum == 2|data$fileNum == 3] = 'control'
data$block[data$fileNum == 5|data$fileNum == 6|data$fileNum == 7] = 'stress'
data <- add_column(data[,1:ncol(data)-1], block = data$block, .after = "time")

data$fileNum <- as.factor(data$fileNum)
data$participantNum <- as.factor(data$participantNum)

# formulas = c('valence ~ time * block', 'arousal ~ time * block')
formulas = c('valence ~ fileNum', 'arousal ~ fileNum')

plotTitles = c('Valence', 'Arousal')

for(i in 2) {
  formula <- paste0(formulas[i], ' + (1|participantNum)')
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
  print("Stats baseline vs control vs stress:")
  # print(emmeans(chosenModel[[1]], pairwise ~ condition , adjust ="fdr", type="response")) # This is the right one for physiological
  print(emmeans(chosenModel[[1]], pairwise ~ fileNum , adjust ="fdr", type="response")) # This is the right one for self-reports
  
  # Plotting
  # dpi=600    #pixels per square inch
  # jpeg(paste0(plotPrefix, "Figure", "_", plotTitles[i], ".jpeg"), width=8*dpi, height=8*dpi, res=dpi)
  par(mfcol = c(1, 1))
  plotAROUSAL <- pirateplot(
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
  )
  # dev.off()
}

