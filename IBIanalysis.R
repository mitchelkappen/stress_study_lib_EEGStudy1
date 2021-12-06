rm(list = ls()) # Clear environment
cat("\014") # Clear console
dev.off() # Clear plot window

library(yarrr)
library(lme4)
library(emmeans)
library(pander)

library(reshape)
library(pander)
library(dplyr)
library(arrow)
library(car)

nAGQ = 0

############################# CARDIO DATA ######################
# Set and Get directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location

IBIdata <-
  as.data.frame(read_parquet("df_tot_merged.parquet"))

# Compute dataframe with relevant variables
data <- data.frame(IBIdata[,c("user", "answered_correctly", "answered_in_time", "Running[Trial]",  "Trial", "Procedure[Block]")], select(IBIdata,contains("IBI_pos")))
groupingVars <- c("pptNum", "answered_correctly", "answered_in_time", "subBlock", "Trial", "Block") # Give easier to use names
names(data)[1:6] <- groupingVars

data$pptNum <- as.factor(data$pptNum)
data$answered_correctly <- as.factor(data$answered_correctly)
data$answered_in_time <- as.factor(data$answered_in_time)
data$subBlock <- as.factor(data$subBlock)
data$Trial <- as.factor(data$Trial)
data$Block <- as.factor(data$Block)

# Something weird going on here still.. 
data <- data[is.na(data$IBI_pos.5) == FALSE, ] # Take out all NA's for IBI's || Using IBI5 because that is where (in addition to generally missing data), some participants start missing data: probably has to do with 6seconds
data <- data[is.na(data$IBI_pos.6) == FALSE, ] # Take out all NA's for IBI's || Using IBI5 because that is where (in addition to generally missing data), some participants start missing data: probably has to do with 6seconds
data <- data[is.na(data$IBI_pos.7) == FALSE, ] # Take out all NA's for IBI's || Using IBI5 because that is where (in addition to generally missing data), some participants start missing data: probably has to do with 6seconds
data <- data[is.na(data$IBI_pos8) == FALSE, ] # Take out all NA's for IBI's || Using IBI5 because that is where (in addition to generally missing data), some participants start missing data: probably has to do with 6seconds

sprintf("Length of all data is: %.0f, and remaining size after removing NA is: %.0f", nrow(IBIdata), nrow(data))
dataBackup <- data #backup the data

# data <- data[,1:(length(dataBackup)-2)] # There are NA's for IBI7 and IBI8 for some participants. Probably has to do with the 6 second stuff. Look into, erase for now

data <- melt(dataBackup, id.vars = groupingVars) # Get it to long format

names(data)[names(data) == "variable"] <- "IBIno"
names(data)[names(data) == "value"] <- "IBIdelta_ms"
data$IBIno <- as.factor(data$IBIno)
levels(data$IBIno) = c("-7","-6","-5","-4","-3","-2", "1", "0", "1", "2", "3", "4", "5", "6", "7", "8")


# Formula time
formulas = c('IBIdelta_ms ~ Block * IBIno', 'arousal ~ fileNum')

plotTitles = c('Control vs Stress', 'Arousal')

for(i in 1) {
  formula <- paste0(formulas[i], ' + (1|pptNum)')
  # Model
  d0.1 <- lmer(formula,data=data)

  modelNames = c(d0.1)
  
  # Model Selection
  tabel <- cbind(AIC(d0.1))
  chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC
  
  Anova(chosenModel[[1]]) # Run Anova, double square brackets because of list properties
  print("Stats baseline vs control vs stress:")
  print(emmeans(chosenModel[[1]], pairwise ~ Block , adjust ="fdr", type="response")) # This is the right one for self-reports
  
  # Plotting
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
    
    xlab = "Sequential IBI",
    ylab = "Delta IBI (ms)"
  )
}


ggplot(data, aes(Time, Value)) + 
  geom_smooth(stat = 'summary', fun.data = mean_cl_quantile)

ggplot(data=data, aes(x=IBIno, y=IBIdelta_ms, group=Block)) +
  geom_line()+
  geom_point()


plottest1 <- ggplot(data, aes(x=IBIno, y=IBIdelta_ms, fill=Block))


plottest1 <- ggplot(data, aes(x=IBIno, y=IBIdelta_ms, fill=Block)) +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.125, aes(ymin=emmean-SE, ymax=emmean+SE)) +
  coord_cartesian(ylim = c(0, 75)) +
  theme_bw(base_size = 14) +
  theme(legend.position="bottom") +
  labs(y = "Counterfactual thinking") +
  labs(x = "Self-critical rumination") +
  ggtitle("A")  +
  geom_signif(annotations = c(formatC("*", digits=3)), y_position = c(55), xmin=c(0.775), xmax=c(1.225), textsize = 6, vjust = 0.5, tip_length = c(0.025, 0.025)) +
  geom_signif(annotations = c(formatC("*", digits=3)), y_position = c(73), xmin=c(2.775), xmax=c(3.225), textsize = 6, vjust = 0.5, tip_length = c(0.025, 0.025))

# Load document where functions are stored
source("functions.R")

#### LINE PLOTTING TIME ####
# First, make a new dataframe with the means per group

library(ggplot2)
# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

# data_summary <- summarySE(data, measurevar="IBIdelta_ms", groupvars=c("Block"))
data_summary <- summarySE(data, measurevar="IBIdelta_ms", groupvars=c("Block","IBIno"))

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

# Standard error of the mean
ggplot(data_summary, aes(x=IBIno, y=IBIdelta_ms, colour=Block, group = Block)) + 
  geom_errorbar(aes(ymin=IBIdelta_ms-se, ymax=IBIdelta_ms+se), width=.1, position=pd) +
  geom_line() +
  geom_point()

# 95% confidence intervals
ggplot(data_summary, aes(x=IBIno, y=IBIdelta_ms, colour=Block, group = Block)) + 
  geom_errorbar(aes(ymin=IBIdelta_ms-ci, ymax=IBIdelta_ms+ci), width=.1, position=pd) +
  geom_line() +
  geom_point()






