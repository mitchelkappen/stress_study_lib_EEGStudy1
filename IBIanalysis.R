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
library(ggplot2)
library(effects)
library(ggsignif)
library(gridExtra) #gridarrange

nAGQ = 0
IBIlength = "big"
############################# CARDIO DATA ######################
# Set and Get directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location

IBIdata <-
  as.data.frame(read_parquet("df_tot_merged_v2_ibi_pos_-2.parquet"))

# Compute dataframe with relevant variables
data <- data.frame(IBIdata[,c("user", "answered_correctly", "answered_in_time", "Running[Trial]",  "Trial", "Procedure[Block]", "sex", "answered_correctly")], select(IBIdata,contains("IBI_pos")))
groupingVars <- c("pptNum", "answered_correctly", "answered_in_time", "subBlock", "Trial", "Block", "Sex", "Correct") # Give easier to use names
names(data)[1:8] <- groupingVars

# Factorize relevant variables and clean up data
data$pptNum <- as.factor(data$pptNum)
data$answered_correctly <- as.factor(data$answered_correctly)
data$answered_in_time <- as.factor(data$answered_in_time)
data$subBlock <- as.factor(data$subBlock)
data$Trial <- as.factor(data$Trial)
data$Block <- as.factor(data$Block)
data$Sex <- as.factor(data$Sex)
data$Correct <- as.factor(data$Correct)
data$answered_in_time <- as.factor(data$answered_in_time)
data$subBlock <- as.factor(data$subBlock)

data <- data[is.na(data$IBI_pos.2) == FALSE, ] # Take out all NA's for IBI's 
data <- data[data$answered_in_time == TRUE, ] # Take out the timed-out trials

sprintf("Length of all data is: %.0f, and remaining size after removing NA is: %.0f", nrow(IBIdata), nrow(data))
dataBackup <- data #backup the data

# data <- data[,1:(length(dataBackup)-2)] # There are NA's for IBI7 and IBI8 for some participants. Probably has to do with the 6 second stuff. Look into, erase for now

data <- melt(dataBackup, id.vars = groupingVars) # Get it to long format
names(data)[names(data) == "variable"] <- "IBIno"
names(data)[names(data) == "value"] <- "IBIdelta_ms"

if(IBIlength == "small"){
  print('hi')
  # Create a plotting variable from IBI-4 for clarity in the viz
  plotdata = data[!(data$IBIno=="IBI_pos.7" | data$IBIno=="IBI_pos.6" | data$IBIno=="IBI_pos.5"), ]
  # But create a stats dataframe where the irrelevant datapoint will not be considered
  data = data[!(data$IBIno=="IBI_pos.7" | data$IBIno=="IBI_pos.6" | data$IBIno=="IBI_pos.5" | data$IBIno=="IBI_pos.4" | data$IBIno=="IBI_pos.3" | data$IBIno=="IBI_pos.2" | data$IBIno=="IBI_pos.1"), ]
} else if (IBIlength == "big"){
  data = data[!(data$IBIno=="IBI_pos.7" | data$IBIno=="IBI_pos.6" | data$IBIno=="IBI_pos.5" | data$IBIno=="IBI_pos.4"), ]
}

levels(data$IBIno) = c("-7","-6","-5","-4","-3","-2", "-1", "0", "1", "2", "3", "4", "5", "6", "7", "8")
data$IBIno <- as.factor(data$IBIno)
data$IBIno <- as.ordered(data$IBIno) 

############ STATISTICS #############

data$subBlock2 = data$subBlock
levels(data$subBlock2) = c("1","2","3","1","2","3")

# Full formula
formula <- IBIdelta_ms ~ Block * IBIno * subBlock2 + (1|pptNum)

# Load document where functions are stored
source("functions.R")

# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
# Make a data summary based on grouping vars

################## PLOTTING ##########
# Plot Settings
# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

### Plot 1 - Control vs Stress
d0.1 <- lmer(formula,data=data) # Fit the lmer
emmeans0.1 <- emmeans(d0.1, pairwise ~ Block | IBIno, adjust ="fdr", type = "response") # Compute a variable containing all emmeans/contrasts
emm0.1 <- summary(emmeans0.1)$emmeans

Anova(d0.1)
plot(effect("Block:IBIno:subBlock2", d0.1)) #just to check
plot(effect("IBIno:subBlock2", d0.1)) #just to check
plot(effect("Block:IBIno", d0.1)) #just to check

print("We see a main effect for control vs stress block, but not an interaction effect of block x IBIno. So we continue with post-hoc testing")
print("Also an effect for subblock is present, but not for subblock x Block.")

if(IBIlength == "small"){
  xplotPosition = 4.1
} else if(IBIlength == "big"){
  xplotPosition = 7.1
}
## LINEPLOT
ggplot(emm0.1, aes(x=IBIno, y=emmean, color=Block)) +
  geom_point(size = 1) + 
  geom_line(aes(group = Block),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  ggtitle("Feedback: Group x IBI")+
  labs(y = "Delta IBI (s)")+
  annotate(geom="text", x=xplotPosition, y=-37.5, label="**", color="#000000")+ #IBI3
  annotate(geom="text", x=xplotPosition + 1, y=-34.5, label="**", color="#000000")+ #IBI4
  annotate(geom="text", x=xplotPosition + 2, y=-26, label="***", color="#000000")+ #IBI5
  annotate(geom="text", x=xplotPosition + 3, y=-17, label="***", color="#000000")+ #IBI6
  annotate(geom="text", x=xplotPosition + 4, y=-12.5, label="**", color="#000000") #IBI7

# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

#### SUBBLOCK PLOTS #### Block:IBIno:subBlock2
emmeans0.2 <- emmeans(d0.1, pairwise ~ Block | IBIno | subBlock2, adjust ="fdr", type = "response") # Compute a variable containing all emmeans/contrasts
emm0.2 <- summary(emmeans0.2)$emmeans
ggplot(emm0.2, aes(x=IBIno, y=emmean, group = 3, colour = Block, linetype = subBlock2, shape = subBlock2)) +
  geom_point(size = 3) + 
  geom_line(aes(group = interaction(Block, subBlock2)),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  ggtitle("Feedback: Group x IBI")+
  labs(y = "Delta IBI (s)")
  # annotate(geom="text", x=xplotPosition, y=-37.5, label="**", color="#000000")+ #IBI3
  # annotate(geom="text", x=xplotPosition + 1, y=-34.5, label="**", color="#000000")+ #IBI4
  # annotate(geom="text", x=xplotPosition + 2, y=-26, label="***", color="#000000")+ #IBI5
  # annotate(geom="text", x=xplotPosition + 3, y=-17, label="***", color="#000000")+ #IBI6
  # annotate(geom="text", x=xplotPosition + 4, y=-12.5, label="**", color="#000000") #IBI7

#### SUBBLOCK PLOTS #### Block:IBIno:subBlock2 ## 1x3 plots
par(mfrow = c(1, 3))

# SubBlock 1:
block1Data = emm0.2[emm0.2$subBlock2 == 1, ]
b1 <- ggplot(block1Data, aes(x=IBIno, y=emmean, group = 1, colour = Block)) +
  geom_point(size = 3) + 
  geom_line(aes(group = Block),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  ggtitle("Feedback: Group x IBI at subBlock 1")+
  labs(y = "Delta IBI (s)") +
  ylim(-53, 14) +
  annotate(geom="text", x=xplotPosition, y=-30, label="**", color="#000000")+ #IBI3
  annotate(geom="text", x=xplotPosition + 1, y=-23.5, label="*", color="#000000")+ #IBI4
  annotate(geom="text", x=xplotPosition + 2, y=-15.5, label="**", color="#000000")+ #IBI5
  annotate(geom="text", x=xplotPosition + 3, y=-9, label="**", color="#000000")+ #IBI6
  annotate(geom="text", x=xplotPosition + 4, y=-8, label="*", color="#000000") #IBI7

# SubBlock 2:
block2Data = emm0.2[emm0.2$subBlock2 == 2, ]
b2 <- ggplot(block2Data, aes(x=IBIno, y=emmean, group = 1, colour = Block)) +
  geom_point(size = 3) + 
  geom_line(aes(group = Block),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  ggtitle("Feedback: Group x IBI at subBlock 2")+
  labs(y = "Delta IBI (s)") +
  ylim(-53, 14) +
  annotate(geom="text", x=xplotPosition + 2, y=-28, label="*", color="#000000")+ #IBI5
  annotate(geom="text", x=xplotPosition + 3, y=-19, label="*", color="#000000")+ #IBI6
  annotate(geom="text", x=xplotPosition + 4, y=-13, label="*", color="#000000") #IBI7

# SubBlock 3:
block3Data = emm0.2[emm0.2$subBlock2 == 3, ]
b3 <- ggplot(block3Data, aes(x=IBIno, y=emmean, group = 1, colour = Block)) +
  geom_point(size = 3) + 
  geom_line(aes(group = Block),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  ggtitle("Feedback: Group x IBI at subBlock 3")+
  labs(y = "Delta IBI (s)") +
  ylim(-53, 14) +
  annotate(geom="text", x=xplotPosition + 3, y=-23.5, label="*", color="#000000") #IBI6



grid.arrange(b1, b2, b3, ncol=3, nrow =1)

















# Correct/incorrect