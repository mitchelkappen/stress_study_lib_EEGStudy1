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
data$fileNum <- ordered(data$fileNum)

# Audio Sample descriptives
t.first <- data[match(unique(data$participantNum), data$participantNum),] # Create dataframe with one line per unique participant 
sprintf("Number of participants: %.f", nrow(t.first))
sprintf("Number of Men: %.f. Number of Women: %.f.", sum(t.first$Sex == 'M') , sum(t.first$Sex == 'F')) 
sprintf("Age, Mean: %.2f, SD: %.2f.", mean(t.first$Age) , sd(t.first$Age))

######## Analysis ########
sum(is.na(data))
backupdata <- data
data <- data[, c("F0semitoneFrom27.5Hz_sma3nz_amean", 
                        "jitterLocal_sma3nz_amean", 
                        "shimmerLocaldB_sma3nz_amean", 
                        "HNRdBACF_sma3nz_amean", 
                        "VoicedSegmentsPerSec", 
                        "MeanVoicedSegmentLengthSec", 
                        "valence", "arousal", "Sex", "condition")]

condition_1_data <- data[data$condition == 'Control', 1:6]
condition_2_data <- data[data$condition == 'Negative', 1:6]

# Scale data
condition_1_scaled <- scale(condition_1_data)
condition_2_scaled <- scale(condition_2_data)

######## PCA ########
# Condition Neutral
# Run PCA
pca_condition_1 <- prcomp(condition_1_scaled, center = TRUE, scale. = TRUE)

# View summary and scree plot
summary(pca_condition_1)
plot(pca_condition_1, type = "l")

# Biplot
biplot(pca_condition_1)

# Condition Negative
# Run PCA
pca_condition_2 <- prcomp(condition_2_scaled, center = TRUE, scale. = TRUE)

# View summary and scree plot
summary(pca_condition_2)
plot(pca_condition_2, type = "l")

# Biplot
biplot(pca_condition_2)
