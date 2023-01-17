######################################
#                                    #
# Analysis of Audio and Self-reports #
#       Social stressor data         #
#                                    #
######################################
# This code uses premade csv for speech variables and self-reports
# Here we perform data cleanup, analysis, and visualisation
# Author: Mitchel Kappen
# 12-4-2022
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
library(effectsize)

library(ggpubr)
library(car)
library(ggplot2)

library(arrow)
library(tibble)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location
source("functions.R") # This is a file in the same directory where you can stash your functions so you can save them there and have them together
options(contrasts = c("contr.sum","contr.poly")) #use this for the p value of the t test

#####  General settings ##### 
nAGQ = 1 # Glmer setting
pvalues = c() # Create a variable to store all p-values to correct later

plotPrefix <- "../figures/" # Define directory to store visualisations

data <- as.data.frame(read_parquet("../loc_data/df_gemaps_func.parquet")) # Read dataframe containing audio features and self-reports

##### Clean data up a bit #####
data$participantNum <- as.factor(data$participantNum)

agesex <- as.data.frame(read.csv("../loc_data/SexAge.csv")) # Read dataframe containing participants' Sex and Age
data <- merge(data, agesex, by = c("participantNum")) # Add demographics to main dataframe

# Add final condition names including 'baseline'
data$condition[data$fileNum == 0] = 'baseline'
data$condition[data$fileNum == 1|data$fileNum == 2|data$fileNum == 3] = 'Control'
data$condition[data$fileNum == 5|data$fileNum == 6|data$fileNum == 7] = 'Negative'
data$condition[data$fileNum == 4|data$fileNum == 8] = 'Rest'

# Get relevant data
dataBackup = data # Backup data so we can go back to this whenever
data = data[data$condition ==  'Control' | data$condition == 'Negative',] # Get only control and negative feedback data (kick out resting blocks)
data = data[data$HNRdBACF_sma3nz_amean > 0, ] # Kick out all the lines with negative HNR, that means the signal is too noisy to trust

# Factorize final relevant variables
data$condition <- as.factor(data$condition)
data$Sex <- as.factor(data$Sex)

# Audio Sample descriptives
t.first <- data[match(unique(data$participantNum), data$participantNum),] # Create dataframe with one line per unique participant 
sprintf("Number of participants: %.f", nrow(t.first))
sprintf("Number of Men: %.f. Number of Women: %.f.", sum(t.first$Sex == 'M') , sum(t.first$Sex == 'F')) 
sprintf("Age, Mean: %.2f, SD: %.2f.", mean(t.first$Age) , sd(t.first$Age))

######## Analysis ########
cbPalette <- c("#56B4E9", "#E69F00") # Define Colorblind proof plotting colors
####### Speech features #######
###### Speech features: F0 ######
formula <- 'F0semitoneFrom27.5Hz_sma3nz_amean ~ condition + Sex + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=dataModel, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=dataModel, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III') #d0.2

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ condition, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

# Effect size
z_to_d(
  z = c(2.756),
  df_error = 1,
  n = 71
)

F0_plot <- audio_pretty_plot(emm0.1, "F0")
F0_plot <- F0_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='**', size=7)
ggsave(F0_plot, file=paste0(plotPrefix, "Figure_F0.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: Jitter ######
formula <- 'jitterLocal_sma3nz_amean ~ condition + (1|participantNum)' # Declare formula | Sex dit not contribute significantly

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=dataModel, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=dataModel, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1) # Other models did not converge
tabel <- cbind(AIC(d0.1))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ condition, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

Jitter_plot <- audio_pretty_plot(emm0.1, "Jitter")
Jitter_plot <- Jitter_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='', size=7)
ggsave(Jitter_plot, file=paste0(plotPrefix, "Figure_Jitter.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: Shimmer ######
formula <- 'shimmerLocaldB_sma3nz_amean ~ condition + Sex + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=dataModel, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=dataModel, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ condition, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

# Effect size
z_to_d(
  z = c(2.881),
  df_error = 1,
  n = 71
)

Shimmer_plot <- audio_pretty_plot(emm0.1, "Shimmer")
Shimmer_plot <- Shimmer_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='**', size=7)
ggsave(Shimmer_plot, file=paste0(plotPrefix, "Figure_Shimmer.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: HNR ######
formula <- 'HNRdBACF_sma3nz_amean ~ condition + Sex + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=dataModel, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=dataModel, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ condition, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

# Effect size
z_to_d(
  z = c(2.858),
  df_error = 1,
  n = 71
)

HNR_plot <- audio_pretty_plot(emm0.1, "HNR")
HNR_plot <- HNR_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='**', size=7)
ggsave(HNR_plot, file=paste0(plotPrefix, "Figure_HNR.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: Voiced per sec ######
formula <- 'VoicedSegmentsPerSec ~ condition + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=dataModel, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=dataModel, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ condition, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

VoicedSeg_plot <- audio_pretty_plot(emm0.1, "Voiced segments per sec")
VoicedSeg_plot <- VoicedSeg_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='', size=7)
ggsave(VoicedSeg_plot, file=paste0(plotPrefix, "Figure_VoicedPerSec.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Speech features: Mean voiced segment length ######
formula <- 'MeanVoicedSegmentLengthSec ~ condition + Sex + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=dataModel, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=dataModel, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ condition, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value)

SegLength_plot <- audio_pretty_plot(emm0.1, "Mean voiced segment length")
SegLength_plot <- SegLength_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='', size=7)
ggsave(SegLength_plot, file=paste0(plotPrefix, "Figure_MeanVoicedLength.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

####### Self-reports #######
###### Self-reports: Valence ######
formula <- 'valence ~ condition + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=dataModel, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=dataModel, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ condition, adjust ="fdr", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts

# Effect size
param_tab <- parameters::model_parameters(chosenModel[[1]])
t_to_d(param_tab$t[2], param_tab$df_error[2])

valence_plot <- audio_pretty_plot(emm0.1, "Valence")
valence_plot <- valence_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 4, label='***', size=7)
ggsave(valence_plot, file=paste0(plotPrefix, "Figure_Valence.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### Self-reports: Arousal ######
formula <- 'arousal ~ condition + (1|participantNum)' # Declare formula

dataModel = data # Ensure correct data is taken
rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=dataModel, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=dataModel, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.3)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ condition, adjust ="fdr", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts

# Effect size
z_to_d(
  z = c(2.116),
  df_error = 1,
  n = 71
)

arousal_plot <- audio_pretty_plot(emm0.1, "Arousal")
arousal_plot <- arousal_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + (max(emm0.1$emmean) - min(emm0.1$emmean)) / 2, label='*', size=7)
ggsave(arousal_plot, file=paste0(plotPrefix, "Figure_Arousal.jpeg"), width = 2000, height = 1500, dpi = 300, units = "px")

###### p adjust #####
names = c('F0', 'Jitter', 'Shimmer', 'HNR', 'VoicedperSec', 'MeanVoicedSegLength')
ps = list()
ps[names] = p.adjust(pvalues, method = "fdr", length(pvalues)) # Create list containing fdr corrected pvalues
