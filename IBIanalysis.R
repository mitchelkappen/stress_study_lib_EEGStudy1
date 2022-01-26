##### Set environment #####
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
library(gridExtra)
library(viridis)

IBIlength = "big"

# Set and Get directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location
BASEPATH <- "D:/Data/EEG_Study_1/aligned_data/features/"
plotPrefix <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/figures/")

##### Loading data ##### 

IBIdata <-
  as.data.frame(read_parquet("df_tot_merged_v2_ibi_pos_-2.parquet"))

agesex <-
  as.data.frame(read.csv(paste0(BASEPATH, "SexAge.csv")))
agesex$user <- agesex$participantNum # Rename pptnum column for succesful merge
agesex <- subset(agesex, select = -c(participantNum, Sex)) # Drop irrelevent columns to not cloud the dataframe

IBIdata <- merge(IBIdata, agesex, by = c("user"))
##### Data cleanup #####
# Compute dataframe with relevant variables
data <- data.frame(IBIdata[,c("user", "answered_correctly", "answered_in_time", "Running[Trial]",  "Trial", "Procedure[Block]", "sex", "Age", "answered_correctly")], select(IBIdata,contains("IBI_pos")))
groupingVars <- c("pptNum", "answered_correctly", "answered_in_time", "subBlock", "Trial", "Block", "Sex", "Age", "Correct") # Give easier to use names
names(data)[1:9] <- groupingVars

# Factorize relevant variables and clean up data
data$pptNum <- as.factor(data$pptNum)
data$answered_correctly <- as.factor(data$answered_correctly)
data$answered_in_time <- as.factor(data$answered_in_time)
data$subBlock <- as.factor(data$subBlock)
data$Trial <- as.factor(data$Trial)
data$Block[data$Block == "Controle"] = "Neutral"
data$Block[data$Block == "Stress"] = "Negative"
data$Block <- factor(data$Block, levels = c("Neutral", "Negative"))
data$Sex <- as.factor(data$Sex)
data$Correct <- as.factor(data$Correct)
data$answered_in_time <- as.factor(data$answered_in_time)
data$subBlock <- as.factor(data$subBlock)

data <- data[is.na(data$IBI_pos.2) == FALSE, ] # Take out all NA's for IBI's 
data <- data[data$answered_in_time == TRUE, ] # Take out the timed-out trials

sprintf("Length of all data is: %.0f, and remaining size after removing NA is: %.0f", nrow(IBIdata), nrow(data))
dataBackup <- data #backup the data

data <- melt(dataBackup, id.vars = groupingVars) # Get it to long format
names(data)[names(data) == "variable"] <- "IBIno"
names(data)[names(data) == "value"] <- "IBIdelta_ms"

if(IBIlength == "small"){
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

######## Analysis ########
##### Statistics ####
data$subBlock2 = data$subBlock # Data backup
levels(data$subBlock2) = c("1","2","3","1","2","3") # Set data levels
# Load document where functions are stored
source("functions.R")

# Full formula
formula <- IBIdelta_ms ~ Block * IBIno + (1|pptNum)

d0.1 <- lmer(formula,data=data) # Fit the lmer

Anova(d0.1, type = 'III')
plot(effect("Block:IBIno", d0.1)) #just to check
# emmeans(d0.1, pairwise ~ Block | IBIno, adjust ="fdr", type = "response")
emmeans0.1 <- emmeans(d0.1, pairwise ~ Block | IBIno, adjust ="fdr", type = "response") # Compute a variable containing all emmeans/contrasts
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts

#### Visualisation ####
pd <- position_dodge(0.05) # To prevent errorbars overlapping, use position_dodge to move them horizontally - move them .05 to the left and right

print("Significant interaction effect block and IBIno ")

if(IBIlength == "small"){
  xplotPosition = 4.1
} else if(IBIlength == "big"){
  xplotPosition = 7.1
}

cbPalette <- c("#56B4E9", "#E69F00")

## LINEPLOT
IBI_plot <- ggplot(emm0.1, aes(x=IBIno, y=emmean, color=Block)) +
  geom_point(size = 2) + 
  geom_line(aes(linetype = Block, group = Block),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_colour_manual(values=cbPalette)+
  scale_linetype_manual(values=c("dashed", "solid")) +
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  labs(y = "Delta IBI (ms)", x = "IBI no.")+
  annotate(geom="text", x=xplotPosition - 1, y=mean(emm0.1$emmean[11:12]), label="*", color="#000000")+ #IBI2
  annotate(geom="text", x=xplotPosition, y=mean(emm0.1$emmean[13:14]), label="**", color="#000000")+ #IBI3
  annotate(geom="text", x=xplotPosition + 1, y=mean(emm0.1$emmean[15:16]), label="**", color="#000000")+ #IBI4
  annotate(geom="text", x=xplotPosition + 2, y=mean(emm0.1$emmean[17:18]), label="***", color="#000000")+ #IBI5
  annotate(geom="text", x=xplotPosition + 3, y=mean(emm0.1$emmean[19:20]), label="***", color="#000000")+ #IBI6
  annotate(geom="text", x=xplotPosition + 4, y=mean(emm0.1$emmean[21:22]), label="***", color="#000000")+ #IBI7
  theme(axis.text.x = element_text(size = 16))+ # X Axis ticks
  theme(axis.text.y = element_text(size = 10))+ # Y axis ticks
  theme(axis.title = element_text(size = 16))+ # Axis titles
  theme(legend.text = element_text(size = 16))+ # Legend text
  theme(legend.title = element_text(size = 14))+ # Legend title
  plot_theme_apa()+
  scale_x_discrete(labels=c("-3", "-2(r)", "-1", "0", "1", "2", "3", "4", "5", "6", "7", "8"))+
  theme(
    axis.text.x=element_text(size=rel(3)),
    axis.text.y=element_text(size=rel(2)),
    axis.title.y=element_text(size=rel(1)),
    axis.title.x = element_text(size=rel(1)),
    # legend.position = "bottom",
    legend.position = c(.8,.85),
    legend.title = element_blank()
  )
IBI_plot
ggsave(IBI_plot, file=paste0(plotPrefix, "Figure_IBI.jpeg"), width = 3000, height = 1500, dpi = 300, units = "px")

##### Descriptives #####
# Get No. of participants
length(match(unique(data$pptNum), data$pptNum))
# Get mean and sd of Age
mean(data$Age[match(unique(data$pptNum), data$pptNum)])
sd(data$Age[match(unique(data$pptNum), data$pptNum)])
# Get No. of women and men
sum(data$Sex[match(unique(data$pptNum), data$pptNum)] == 'F')
sum(data$Sex[match(unique(data$pptNum), data$pptNum)] == 'M')