##### Set environment #####
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
library(effects)

library(ggpubr)
library(car)
library(ggplot2)

library(arrow)
library(tibble)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location
source("functions.R") # This is a file in the same directory where you can stash your functions so you can save them there and have them together

#####  General settings ##### 
nAGQ = 1 # Set to 1 for eventual analysis
pvalues = c() # Create a variable to store all p-values to correct later

BASEPATH <- "D:/Data/EEG_Study_1/aligned_data/features/"
plotPrefix <- paste0(BASEPATH, "figures/")
# plotPrefix <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/figures/")

data <-
  as.data.frame(read.csv(paste0(BASEPATH, "dataComplete.csv"))) # This contains everything except for the ones that were too bad -- see Excel drive

##### Clean data up a bit #####
data <- data[1:(length(data)-8)] # Remove the physiology measures here, are treated in other script

data$participantNum <- as.factor(data$participantNum)

agesex <-
  as.data.frame(read.csv(paste0(BASEPATH, "SexAge.csv")))

data <- merge(data, agesex, by = c("participantNum"))

# Add columns for subblock 1-2-3
data$time[data$fileNum == 0] = 0
data$time[data$fileNum == 1|data$fileNum == 5] = 1
data$time[data$fileNum == 2|data$fileNum == 6] = 2
data$time[data$fileNum == 3|data$fileNum == 7] = 3
data <- add_column(data[,1:ncol(data)-1], time = data$time, .after = "fileNum")

# Add columns for block 0-1-2
data$block[data$fileNum == 0] = 'baseline'
data$block[data$fileNum == 1|data$fileNum == 2|data$fileNum == 3] = 'Neutral'
data$block[data$fileNum == 5|data$fileNum == 6|data$fileNum == 7] = 'Negative'
data <- add_column(data[,1:ncol(data)-1], block = data$block, .after = "time")

data$time <- as.factor(data$time)
data$block <- as.factor(data$block)
data$block <- relevel(data$block, "Neutral")
data$block <- relevel(data$block, "baseline")
data$Sex <- as.factor(data$Sex)

# Take out all NA's values (Resting blocks)
fullData <- data
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

# Sample descriptives
t.first <- data[match(unique(data$participantNum), data$participantNum),] # Create dataframe with one line per unique participant 
sprintf("Number of participants: %.f", nrow(t.first))
sprintf("Number of Men: %.f. Number of Women: %.f.", sum(t.first$Sex == 'M') , sum(t.first$Sex == 'F')) 
sprintf("Age, Mean: %.2f, SD: %.2f.", mean(t.first$Age) , sd(t.first$Age))

# Get relevant data
dataBackup = data # Backup data so we can go back to this whenever
data = data[data$block ==  'Neutral' | data$block == 'Negative',] # Get only block 1 and 2 data

data = data[data$HNRdBACF_sma3nz_amean > 0, ] # Kick out all the lines with negative HNR, that means the signal is too noisy to trust
######## Analysis ########
cbPalette <- c("#56B4E9", "#E69F00") # Colorblind plotting colors
####### Speech features #######
###### Speech features: F0 ######
formula <- 'F0semitoneFrom27.5Hz_sma3nz_amean ~ block + Sex + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ block, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

F0_plot <- audio_pretty_plot(emm0.1, "F0")
F0_plot <- F0_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='**', size=7)
ggsave(F0_plot, file=paste0(plotPrefix, "Figure_F0.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: Jitter ######
formula <- 'jitterLocal_sma3nz_amean ~ block + Sex + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1) # Other models did not converge
tabel <- cbind(AIC(d0.1))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ block, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

Jitter_plot <- audio_pretty_plot(emm0.1, "Jitter")
Jitter_plot <- Jitter_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='', size=7)
ggsave(Jitter_plot, file=paste0(plotPrefix, "Figure_Jitter.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: Shimmer ######
formula <- 'shimmerLocaldB_sma3nz_amean ~ block + Sex + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1) # Other models did not converge
tabel <- cbind(AIC(d0.1))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ block, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

Shimmer_plot <- audio_pretty_plot(emm0.1, "Shimmer")
Shimmer_plot <- Shimmer_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='**', size=7)
ggsave(Shimmer_plot, file=paste0(plotPrefix, "Figure_Shimmer.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: HNR ######
formula <- 'HNRdBACF_sma3nz_amean ~ block + Sex + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ block, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

HNR_plot <- audio_pretty_plot(emm0.1, "HNR")
HNR_plot <- HNR_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='**', size=7)
ggsave(HNR_plot, file=paste0(plotPrefix, "Figure_HNR.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: Voiced per sec ######
formula <- 'VoicedSegmentsPerSec ~ block + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ block, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

VoicedSeg_plot <- audio_pretty_plot(emm0.1, "Voiced segments per sec")
VoicedSeg_plot <- VoicedSeg_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='', size=7)
ggsave(VoicedSeg_plot, file=paste0(plotPrefix, "Figure_VoicedPerSec.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: Mean voiced segment length ######
formula <- 'MeanVoicedSegmentLengthSec ~ block + Sex + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1) # Other models did not converge
tabel <- cbind(AIC(d0.1))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ block, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

SegLength_plot <- audio_pretty_plot(emm0.1, "Mean voiced segment length")
SegLength_plot <- SegLength_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='', size=7)
ggsave(SegLength_plot, file=paste0(plotPrefix, "Figure_MeanVoicedLength.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

####### Self-reports #######
###### Self-reports: Valence ######
formula <- 'valence ~ block + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ block, adjust ="fdr", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts

valence_plot <- audio_pretty_plot(emm0.1, "Valence")
valence_plot <- valence_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 4, label='***', size=7)
ggsave(valence_plot, file=paste0(plotPrefix, "Figure_Valence.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Self-reports: Arousal ######
formula <- 'arousal ~ block + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ block, adjust ="fdr", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts

arousal_plot <- audio_pretty_plot(emm0.1, "Arousal")
arousal_plot <- arousal_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='*', size=7)
ggsave(arousal_plot, file=paste0(plotPrefix, "Figure_Arousal.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### p adjust #####
names = c('F0', 'Jitter', 'Shimmer', 'HNR', 'VoicedperSec', 'MeanVoicedSegLength')
ps = list()
ps[names] = p.adjust(pvalues, method = "fdr", length(pvalues)) # Create list containing fdr corrected pvalues
