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
plotPrefix = "Plots/"

# Set and Get directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location

# BASE_PATH = "Z:/ghepmk_data/2020_EEGStudy1/"
# DATA_PATH = paste0(BASE_PATH,"aligned_data")
# FEATURES_PATH = paste0(DATA_PATH,"/features_gemaps/")

# Load the csv file

data <-
  as.data.frame(read_parquet("user_data_per_block.parquet"))

# dataSmall <- data[complete.cases(data$HR), ]
data <- data[ , -which(names(data) %in% c("start_ts","stop_ts"))]

data$block <- as.factor(data$block)
data$ptcpt_id <- as.factor(data$ptcpt_id)

formulas = c('HR ~ block', 'HRV ~ block', 'SCR_RATE ~ block', 'SCR_AMPL ~ block')

plotTitles = c('HR', 'HRV', 'SCR_RATE', 'SCR_AMPL')

data_temp <- data #to get back to

# for(i in 1:length(formulas)) {
for(i in 4) {
  data <- data[complete.cases(get(names(data)[i + 2], data)), ] # get only complete cases for this specific variable
  # data <- data[is.na(get(names(data)[i + 2], data)), i + 2] <- 0 # First, Turn NA into Zero
  data <- data[get(names(data)[i + 2], data) != 0, ] # Now remove all rows with zeroes
  
  formula <- paste0(formulas[i], ' + (1|ptcpt_id)')
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
  )
  # abline(lm(formulas[i], data=data), lwd=4, lty=2, col = "red")
  
  # Secondary points and axis:
  # mtext("fileNum",1,line=1,at=0.2,col="red")
  
  # mtext("Experiment phase",1,line=3,at=0.2,col="blue")
  # axis(side=1, at=c(1:9), line = 3, labels=c('prerest','control1','control2','control3','midrest','stress1','stress2','stress3','postrest' ))
  # axis(side=1, at=c(1:6), line = 3, labels=c('control1','control2','control3','stress1','stress2','stress3' ))

  # dev.off()
   data <- data_temp # make sure you have all the data again in next round
}
