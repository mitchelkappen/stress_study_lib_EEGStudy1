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

nAGQ = 1 # Set to 1 for eventual analysis

BASEPATH <- "D:/Data/EEG_Study_1/aligned_data/features/"
plotPrefix <- paste0(BASEPATH, "figures/")
plotPrefix <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/figures/")

data <-
  as.data.frame(read.csv(paste0(BASEPATH, "dataComplete.csv"))) # This contains everything except for the ones that were too bad -- see Excel drive

data <- data[1:(length(data)-8)] # Remove the physiology measures here, are treated in other script

data$participantNum <- as.factor(data$participantNum)

agesex <-
  as.data.frame(read.csv(paste0(BASEPATH, "SexAge.csv")))

data <- merge(data, agesex, by = c("participantNum"))
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

# Get relevant data
dataBackup = data # Backup data so we can go back to this whenever
data = data[data$block ==  'control' | data$block == 'stress',] # Get only block 1 and 2 data

data = data[data$HNRdBACF_sma3nz_amean > 0, ] # Kick out all the lines with negative HNR, that means the signal is too noisy to trust
######## Analysis ########
####### Speech features #######
###### Speech features: F0 ######
formula <- 'F0semitoneFrom27.5Hz_sma3nz_amean ~ block +(1|participantNum)' # FirstMenstrual had zero effect so was removed from the model | Order showed no effect and was removed from model

dataModel = data

rm(d0.1, d0.2, d0.3) # Just to be sure you're not comparing former models for this comparison

d0.1 <- lmer(formula,data=dataModel)
d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

modelNames = c(d0.1,d0.2,d0.3)

# Model Selection
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))

chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')
plot(effect("block", chosenModel[[1]], se = TRUE))

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ block, adjust ="fdr", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts

max_y<-max(data$F0semitoneFrom27.5Hz_sma3nz_amean) # max PSS is 33
emm0.2 <- data.frame('block'=emm0.1$block, 'F0'= emm0.1$emmean) #dataframe with all the emmeans

F0_plot <- audio_pretty_plot(emm0.1, "F0")
F0_plot + annotate('text', x=1.5, y=mean(emm0.1$emmean) + mean(emm0.1$emmean)/750, label='**', size=7)

ggplot()+
  ggtitle('F0 ~ block')+ #title 
  geom_flat_violin(data= data, aes(x= block, y= F0semitoneFrom27.5Hz_sma3nz_amean, fill=block),position = position_nudge(x =.3, y = 0), adjust = 1.5, alpha = .5, colour = NA)+ # flat violin distribution, .3 points to the right. alpha=.5 so see-through
  geom_boxplot(data= data, aes(x=block, y=F0semitoneFrom27.5Hz_sma3nz_amean, fill=block), outlier.shape=NA, alpha=.5, width=.3, colour='black')+ #boxplot, see through, no outline, 
  geom_point(data= emm0.2, aes(x = block, y = F0, fill=block), position= position_dodge(0.3), size=4) #points representing the emmeans

ggplot(data, aes(x = block, y = F0semitoneFrom27.5Hz_sma3nz_amean)) +
  geom_flat_violin(aes(fill=block),position = position_nudge(x =.2, y = 0), alpha=.5, adjust = 1.5, colour = NA)+
  # geom_point(aes(colour=PMS, fill=PMS),position=position_jitter(width=.15), alpha=.5, size=.25)+
  geom_boxplot(aes(x = block, y = F0semitoneFrom27.5Hz_sma3nz_amean, fill = block),outlier.shape= NA, alpha = .45, width = .1, colour = "black")+
  geom_point(data= emm0.2, aes(x = block, y = F0, fill=block),outlier.shape= NA, width = .5, size=4)+
  ggtitle('DASS_Anxiety~PMS')+
  geom_segment(aes(x = 1, y=max_y, xend= 2, yend=max_y), size= 1)+
  annotate('text', x=1.5, y=max_y+0.3, label='***', size=7)+
  geom_segment(aes(x = 2, y=max_y+1, xend= 3, yend=max_y+1), size= 1)+
  annotate('text', x=2.5, y=max_y+ 1.3, label='**', size=7)+
  geom_segment(aes(x = 1, y=max_y+2, xend= 3, yend=max_y+2), size= 1)+
  annotate('text', x=2, y=max_y+2.3, label='***', size=7)+
  scale_fill_manual(values = c("blue", 'red', 'purple'),
                    name='',labels=c('noPMS \n n=128 ', 'PMS \n n=74', 'PMDD \n n=35'))+
  guides(fill = guide_legend(reverse=TRUE))+
  theme(
    legend.key.size=unit(1.3, 'cm'),
    legend.text=element_text(size=13),
    plot.title = element_text(size=rel(2)),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major.y = element_line( size=.1, color="#dedede" ),
    axis.text.x=element_text(size=rel(1.5)),
    axis.title.y=element_text(size=rel(1.4)),
    axis.title.x = element_blank()) 

par(mfcol = c(1, 1))
plotAROUSAL <- pirateplot(
  formula = 'F0semitoneFrom27.5Hz_sma3nz_amean ~ block',
  data = data,
  theme = 1,
  pal = "info",
  main = 'pirate',
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

library(sjPlot)
plot_model(
  chosenModel[[1]], 
  colors = "Accent", 
  show.values = TRUE,
  value.offset = .4,
  value.size = 4,
  dot.size = 3,
  line.size = 1.5,
  vline.color = "blue",
  width = 1.5
)

plot_model(chosenModel[[1]], show.values = TRUE, value.offset = .3)
################
# formulas = c('F0semitoneFrom27.5Hz_sma3nz_amean ~ time', 'jitterLocal_sma3nz_amean ~ time', 'shimmerLocaldB_sma3nz_amean ~ time', 'HNRdBACF_sma3nz_amean ~ time', 'VoicedSegmentsPerSec ~ time', 'MeanVoicedSegmentLengthSec ~ time')
formulas = c('F0semitoneFrom27.5Hz_sma3nz_amean ~ block', 'F0BaselineCorrected ~ block', 'jitterLocal_sma3nz_amean ~ block', 'shimmerLocaldB_sma3nz_amean ~ block', 'HNRdBACF_sma3nz_amean ~ block', 'VoicedSegmentsPerSec ~ block', 'MeanVoicedSegmentLengthSec ~ block', 'valence ~ block', 'arousal ~ block')
# formulas = c('valence ~ block', 'arousal ~ block')
# formulas = c('valence ~ fileNum', 'arousal ~ fileNum')

plotTitles = c('F0', 'F0 - baselinecorrected', 'Jitter', 'Shimmer', 'HNR', 'VoicedPerSec', 'MeanVoicedLength', 'Valence', 'Arousal')
# plotTitles = c('Valence', 'Arousal')

# for(i in 1:length(formulas)) {
for(i in 1) {
  formula <- paste0(formulas[i], ' + (1|Sex) + (1|participantNum)')
  formula <- paste0(formulas[i], ' + (1|participantNum)')
  formula <- paste0(formulas[i], ' + Sex + Age + (1|participantNum)')
  
  # Model
  d0.1 <- lmer(formula,data=data)

  d0.4 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

  d0.7 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

  # modelNames = c(d0.1,d0.4,d0.7)
  modelNames = c(d0.1)
  
  # Model Selection
  # tabel <- cbind(AIC(d0.1), AIC(d0.4), AIC(d0.7))
  tabel <- cbind(AIC(d0.1))
  
  chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC
  
  print(Anova(chosenModel[[1]], type = 'III')) # Run Anova, double square brackets because of list properties
  plot(effect("block", chosenModel[[1]])) #just to check
  print("Stats control vs stress:")
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
  
  # dev.off()
}
