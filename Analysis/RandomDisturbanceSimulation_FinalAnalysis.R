############## Final Analysis script ####
#Full working script for feature selection, Random Forest Regression with importance scaling and significance, cross-validation
#This was the final script/analysis used for Eco Apps and Ch2 of dissertation (6/6/2023)

#####packages
#install.packages('randomForest')
#install.packages('rfUtilities')
#install.packages('spatialEco')
#install.packages('ranger')
#install.packages('pdp')
#install.packages('spm')
#install.packages('ggpubr')
library('dplyr')
library('data.table')
library('rlang')
library('randomForest')
library('rfUtilities')
library('spatialEco') #collinearity
library('ranger') #running random forests - if librarying this it masks importance()!!!
library('pdp')
library('spm')
library('gridExtra') #multi panel figures
library('ggplot2') # graphics
library('ggpubr') #ggarrange

###Functions
source("./rf.ImpScale.R") #function script in same folder as this one

pfun <- function(object, newdata, i = 2) {
  predict(object, data = newdata)$predictions[,i]
}

rmse <- function(y,x) { sqrt(mean((y - x)^2)) }
set.seed(666)
#***********************************************

## Read in and format data - verified works correct for OGcorrected
Combined <- read.csv(file="C:/Users/bca17002/Desktop/Ch3_DisturbanceSimulation/Data_simExperiment/AllMasterDatasets/Corrected/SimulationExperiment_ExperimentRuns_MasterDataset_OGcorrected.csv")
#Combined <- Combined[1:3000, ] #take only NEGL
#Combined <- Combined[3001:6000, ] #take only Moss
#Combined <- Combined[6001:9000, ] #take only FentonWest
#Combined <- Combined[9001:12000, ] #take only FentonEast
Site <- Combined[ ,1]
Plot <- Combined[ ,2]
Severity <- Combined[ ,4]
predictorSet1 <- Combined[ ,8:80]
DeltaRugosity <- Combined[ ,106]
Combined <- cbind(Site, Plot, Severity, predictorSet1, DeltaRugosity) #uncomment if including site and plot
#Combined <- cbind(Plot, Severity, predictorSet1, DeltaRugosity) #if leaving out plot
sum(is.na(Combined)) # check for NAs

#mutate the Plot column so plots indicated by numbers (for RFR) - add "" around number if want to make factor - verified
#this converts to true numeric value (not factor) which verified is needed here
Combined <- Combined %>% 
  mutate(Plot = case_when(
    Plot == "W" & Site == "NorthEagleville" ~ 1,
    Plot == "WC" & Site == "NorthEagleville"  ~ 2,
    Plot == "EC" & Site == "NorthEagleville"  ~ 3,
    Plot == "S" & Site == "Moss"  ~ 4,
    Plot == "SC" & Site == "Moss" ~ 5,
    Plot == "NW" & Site == "Moss" ~ 6,
    Plot == "NW" & Site == "FentonWest" ~ 7,
    Plot == "SW" & Site == "FentonWest" ~ 8,
    Plot == "E" & Site == "FentonWest" ~ 9,
    Plot == "NE" & Site == "FentonEast" ~ 10,
    Plot == "NW" & Site == "FentonEast" ~ 11,
    Plot == "SW" & Site == "FentonEast" ~ 12
  ))

#mutate the Site column so sites indicated by numbers (for RFR) - verified ########changed so non-numeric 4/22/23
#this converts to true numeric value (not factor) which verified is needed here
Combined <- Combined %>% 
  mutate(Site = case_when(
    Site == "NorthEagleville" ~ "North Eagleville",
    Site == "Moss" ~ "Moss",
    Site == "FentonWest" ~ "Fenton West",
    Site == "FentonEast" ~ "Fenton East"
  ))

Combined$Site <- as.factor(Combined$Site) ############
Combined$Plot <-as.factor(Combined$Plot) #############
is.factor(Combined$Site) #checks if Site a factor #############
is.factor(Combined$Plot) #checks if Plot a factor ############
sum(is.na(Combined)) # check for NAs

## Feature selection/elimination
#dat <- Combined[ ,1:76] #subset of only predictors to be evaluated for colinearity/multicol
dat <- Combined[ ,3:76] # if excluding site and plot (continuous predictor columns only)
#dat <- Combined[ ,2:76] # if excluding site

all.vars <- colnames(dat) #character string indicating names of all predictors to be evaluated for colinearity/multicol

# remove collinearity - test for collinearity # p actually reflects pearsons r (if nonlinear = FALSE) - automatically removes vars with higher mean pearson r
cl.vars <- spatialEco::collinear(dat[,all.vars], p = 0.85, 
                                 nonlinear = FALSE, p.value = 0.05)
if(length(cl.vars) > 0) {
  cat("Dropping colinear variables; ", cl.vars, "\n")  
  dat2 <- dat[,-which(names(dat) %in% cl.vars)]
  vars <- all.vars[-which(all.vars %in% cl.vars)]
} else {
  cat("No collinear variable found", "\n")
}
# multicolinearity - dat2 is new dataframe with collinear removed - tests for multicollinearity
#p originally 0.005 but changed as had 13 preds in final model
( mc <- multi.collinear(dat2[,vars], p=0.05,  perm = TRUE, n = 1000) ) #vignette (rfUtils) suggests p=0.05 for large number variables
mc.vars <- mc[which( mc$frequency > 5 ),]$variables 
if(length(mc.vars) > 0) { 
  dat2 <- dat2[,-which(names(dat2) %in% mc.vars)]
} else {
  cat("No multicollinear variable found", "\n")
}

#only need to run this line if multicollinearity identified
vars <- colnames(dat2) #character string indicating names of remaining predictors after colinear and multicol removed

#only need to run this line if re-adding Site and/or Plot after col/multicol (needed if factors)
vars <- c("Site", vars) #c("Site", "Plot", vars) #############

#if want to verify no NAs introduced
sum(is.na(dat2))
#Combined <- na.omit(Combined) #no NAs in this data so non-issue

# set up/run initial Random Forest Model (all preds except dropped colinear/multicollinear) as final step in feature selection/elimination
b = 2500 # number of bootstraps

#Split data by severities (optional for Sankey diagram test)
#Combined <- subset(Combined, Severity <= 0.420) #take only severities for slope1 (seg regression)
#Combined <- subset(Combined, between(Severity, 0.420, 0.571)) #take only severities for slope2 (seg regression)
#Combined <- subset(Combined, between(Severity, 0.571, 0.747)) #take only severities for slope3 (seg regression)
#Combined <- subset(Combined, Severity >= 0.747) #take only severities for slope4 (seg regression)

#Find optimal m for mtry - based on recommendations in Liaw & Wiener 2002 optimize by testing default, half, and double and pick lowest OOB
#line below does this 
optimalM1 <- tuneRF(x = Combined[,vars], y = Combined[,"DeltaRugosity"], ntreeTry = b, improve = 0.05, trace = TRUE, plot = TRUE, doBest = FALSE)
#m1 <- min(optimalM1$OOBError) 
#optimalM1 <- as.data.frame(arrange(optimalM1, OOBError))
m1 <- 12 #UPDATE FOR EACH RESPONSE VARIABLE - mtry=12 for DeltaRugosity (all obs), use mtry=24 for mean.height & mean.std (took forever in tuneRF)

# Evaluate importance significance and remove anything insignificant or with importance <= 0.1 (arbitrary)
# imp.reg appears to run an initial RF regression with all predictors remaining after removing colinear and multicolinear
( imp.reg <- ranger(x=Combined[,vars], y=Combined[,"DeltaRugosity"],
                    importance = "permutation", mtry = m1, replace = TRUE,
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=100)) ) #calls external function and scales importance scores/computes p-values by permutation (not by SE) considering non-collinear/multicollinear predictors (n=50-100 recommended for stable results by Altmann et al 2010) - relative importance among predictors after scaling is identical to pre-scaling as coded (seemed ok as both done via permutation) - took nearly an hour with full data/final rf params
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.05,]$parameter ) # removes anything with negative value and keeps everything positive that's significant (originally at 0.1 but testing 0.05)
#imp.vars <- c("Site", "speciesCount", "propRemStemsCDorD", "propBAremSuppressed", "propBAremD", "meanHeightCrownBase", "clarkEvans", "deltaPropBACD", "deltaSdDBH") #uncomment to force specific impvars - OTHERWISE USE ABOVE LINE INSTEAD

## Train the final Random Forest Model
# Split data into train and test data
#train.idx <- sample(nrow(Combined), 4/5 * nrow(Combined))
#Combined.train <- Combined[train.idx, ]
Combined.train <- Combined #if want to train using all observations (12,000) - SEEMS MOST APPROPRIATE FOR MAKING STATISTICAL INFERENCE (DISSERTATION)
#Combined.test <- Combined[-train.idx, ]
#verify train and test sets unique
#reCombined <- rbind.data.frame(Combined.train, Combined.test)
#nrow(unique(reCombined)) # should be 3000 if 1 site and 12000 if all sites

#Find optimal m for mtry - based on recommendations in Liaw & Wiener 2002 optimize by testing default, half, and double and pick lowest OOB
#line below does this
optimalM2 <- tuneRF(x = Combined.train[,imp.vars], y = Combined.train[,"DeltaRugosity"], ntreeTry = b, improve = 0.05, trace = TRUE, plot = TRUE, doBest = FALSE)
#m2 <- min(optimalM2$OOBError) 
#optimalM2 <- as.data.frame(arrange(optimalM2, OOBError))
m2 <- 3 #UPDATE FOR EACH RESPONSE VARIABLE - mtry=3 for DeltaRugosity (all obs)

# Fit (train) model - trains the final model based on training data subset
#INPUT ALL DATA OR TRAINING HERE?
( deltaRug.fit <- ranger(x=Combined.train[,imp.vars], y=Combined.train[,"DeltaRugosity"],
                         importance = "permutation", num.trees=b, mtry = m2, replace = TRUE,
                         keep.inbag=TRUE) )
cat("R-square of model:", deltaRug.fit$r.squared, "\n") #model R2 (coefficient of determination; amount variance explained based on oob data) - REPORT
cat("Prediction RMSE:", rmse(Combined.train[,"DeltaRugosity"], deltaRug.fit$predictions), "\n") #prediction root mean square error (predicted(based on oob data) vs observed) - REPORT
cat("Prediction correlation:", cor(Combined.train[,"DeltaRugosity"], deltaRug.fit$predictions), "\n") #pearson's r of model predicted(based on oob) vs observed - REPORT
cat("Prediction MSE:", deltaRug.fit$prediction.error, "\n") #(bca) prediction mean square error (overall oob prediction error) (predicted(based on oob data) vs observed) - REPORT

## Scale importance for final model and plot it/export scaled importance and p-values each predictor as csv
( imp <- na.omit(rf.ImpScale(deltaRug.fit, scaling="p", n=100)) ) #calls external function and scales importance scores/computes p-values by permutation (not by SE) considering final predictors only (n=50-100 recommended for stable results by Altmann et al 2010) - relative importance among predictors after scaling is identical to pre-scaling as coded (seemed ok as both done via permutation)
p <- imp[which(imp$parameter %in% imp.vars),]   
p <- p[order(p$importance),] 
deltaRugLabs <- c("Species Richness of Removals", "Change in Standard Deviation DBH", "Prop. Plot Basal Area Removed Dominant", "Mean Height-to-Crown of Removals", "Prop. Removed Stems CD or D", "Change in Prop. Plot Basal Area Codominant",  "Clark-Evans Index Removed Stems (X-Y)", "Site", "Prop. Plot Basal Area Removed Suppressed")
#vaiLabs <- c("Change in Minimum DBH", "Prop. Plot Basal Area Removed Dominant", "Species Richness of Removals", "Change in Prop. Plot Basal Area Intermediate", "Fisher Alpha DBH of Removals", "Change in Prop. Plot Basal Area Codominant", "Change in Prop. Plot Basal Area Overtopped", "Site", "Mean Height-to-Crown of Removals")
#meanHeightLabs <- c("Prop. Plot Basal Area Removed Suppressed", "Conifer Fraction of Removals", "Clark-Evans Index Removed Stems (X-Y)", "Shannon Diversity Spp of Removals", "Change in Prop. Plot Basal Area Intermediate", "Change in Fisher Alpha DBH", "Change in Range DBH", "Mean Height-to-Crown of Removals", "Change in Prop. Plot Basal Area Overtopped", "Site", "Prop. Plot Basal Area Removed Dominant")
#meanStdLabs <- c("Prop. Plot Basal Area Removed Dominant", "Change in Conifer Fraction", "Prop. Plot Basal Area Removed Snag", "Mean Height-to-Crown of Removals", "Site", "Shannon Diversity Spp of Removals", "Change in Prop. Plot Basal Area Overtopped", "Change in Prop. Plot Basal Area Intermediate", "Prop. Plot Basal Area Removed Suppressed", "Conifer Fraction of Removals", "Change in Prop. Plot Basal Area Codominant")
#porosityLabs <- c("Site", "Mean Height-to-Crown of Removals", "Clark-Evans Index Removed Stems (X-Y)", "Prop. Plot Basal Area Removed Dominant", "Species Richness of Removals", "Change in Fisher Alpha Mean Crown Diameter", "Conifer Fraction of Removals", "Change in Contagion Index Stems (X-Y)", "Change in Fisher Alpha DBH", "Fisher Alpha DBH of Removals", "Prop. Plot Basal Area Removed Suppressed", "Change in Prop. Plot Basal Area Overtopped", "Change in Prop. Plot Basal Area Intermediate", "Change in Prop. Plot Basal Area Codominant")
#check to see if order correct
p2 <- cbind(p, deltaRugLabs)
#plot it
x11()
dotchart(p2$importance, labels = p2$deltaRugLabs, lcolor="darkgray",
         main="Scaled Variable Importance To \u394 Rugosity", pch=19, pt.cex = 1.7)  

#export a matching csv if need to make multipanel version
fwrite(p2, file = "C:/Users/bca17002/Desktop/Ch2_FinalFiguresAndStatsOutputRFR/DeltaRugosity_scaledImportanceFinalPredictorsOnly.csv")

## Output raw partial dependency (as yhat) if needed (or alternative way of plotting PDP plots (see below) if plot = TRUE)
#works same whether specify type = "regression" or not; train = SHOULD MATCH WHATEVER USED TO TRAIN MODEL INDICATED HERE
#for partial - chull = FALSE by default so should match ggplot pdp plots
#took nearly an hour to plot with full data/final rf params
pd=list()
x11()
for(i in 1:length(imp.vars)){
  pd[[i]] <- partial(deltaRug.fit, pred.var=imp.vars[i], 
                     train=Combined.train[,imp.vars],
                     main=paste0("Partial dependency of ", imp.vars[i]),
                     plot = FALSE, ice=FALSE) #change to plot = TRUE if want to plot - otherwise outputs raw yhat relationships to csv
  pathToOutput <- "C:/Users/bca17002/Desktop/Ch2_FinalFiguresAndStatsOutputRFR/"
  fwrite(pd[[i]], file = paste0(pathToOutput, "DeltaRugosity_", i, ".csv"))
}

## Partial dependency plots (1 way) with autoplot (uses ggplot2) 
#probably justified to use chull=FALSE for 1 way plots as long as all data used to train RFR models (no extrapolation)
gg1 <- deltaRug.fit %>%
  partial(pred.var = "propBAremSuppressed", train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Prop. Plot Basal Area Removed Suppressed", ylab = "Influence (\u0177)") + geom_line(lwd = 1.7) + ylim(0, 2.8) + theme_classic(base_size = 16)

gg2 <- deltaRug.fit %>%
  partial(pred.var = "Site", train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Site", ylab = "Influence (\u0177)") + geom_line(lwd = 1.7) + ylim(0, 2.8) + theme_classic(base_size = 16)

gg3 <- deltaRug.fit %>%
  partial(pred.var = "clarkEvans", train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Clark-Evans Index Removed Stems (X-Y)", ylab = "Influence (\u0177)") + geom_line(lwd = 1.7) + ylim(0, 2.8) + theme_classic(base_size = 16)

gg4 <- deltaRug.fit %>%
  partial(pred.var = "deltaPropBACD", train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "\u394 Prop. Plot Basal Area Codominant", ylab = "Influence (\u0177)") + geom_line(lwd = 1.7) + ylim(0, 2.8) + theme_classic(base_size = 16)

gg5 <- deltaRug.fit %>%
  partial(pred.var = "propRemStemsCDorD", train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Prop. Removed Stems CD or D", ylab = "Influence (\u0177)") + geom_line(lwd = 1.7) + ylim(0, 2.8) + theme_classic(base_size = 16)

gg6 <- deltaRug.fit %>%
  partial(pred.var = "meanHeightCrownBase", train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Mean Height-to-Crown of Removals", ylab = "Influence (\u0177)") + geom_line(lwd = 1.7) + ylim(0, 2.8) + theme_classic(base_size = 16)

gg7 <- deltaRug.fit %>%
  partial(pred.var = "propBAremD", train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Prop. Plot Basal Area Removed Dominant", ylab = "Influence (\u0177)") + geom_line(lwd = 1.7) + ylim(0, 2.8) + theme_classic(base_size = 16)

gg8 <- deltaRug.fit %>%
  partial(pred.var = "deltaSdDBH", train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "\u394 Standard Deviation DBH", ylab = "Influence (\u0177)") + geom_line(lwd = 1.7) + ylim(0, 2.8) + theme_classic(base_size = 16)

gg9 <- deltaRug.fit %>%
  partial(pred.var = "speciesCount", train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Species Richness of Removals", ylab = "Influence (\u0177)") + geom_line(lwd = 1.7) + ylim(0, 2.8) + theme_classic(base_size = 16)

#ggarrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, ncol = 3, nrow = 3, align = "hv") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

multiPdp1 <- ggarrange(gg1, gg2, gg4, gg3, gg5, gg6, gg8, gg7, gg9, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), font.label = list(size = 24), ncol = 3, nrow = 3, align = "hv") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
annotate_figure(multiPdp1, fig.lab = "", fig.lab.pos = "top", fig.lab.face = "bold")

## Two-way partial dependency (with site) - combined effects of site and each important predictor (GOOD AS IS)
#must use chull = FALSE to get all 4 sites
#probably more justified to use chull=TRUE for 2 way plot - reduces risk of extrapolating (full range of each predictor likely not realized in each site as sites differ)
#option and direction from viridis (ignored)
#if need to re-plot can try adding "+ theme(panel.background = element_rect(fill = "white", color = "black"), plot.margin = margin(0.15, 0.15, 0.15, 0.15, "cm"))" to each plot to correct margins for g1 (the 0.025 value is slightly cut off in current plot)
g1 <- deltaRug.fit %>%
  partial(pred.var = c("propBAremSuppressed", "Site"), train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Prop. Plot Basal Area Removed Suppressed", ylab = "Influence (\u0177)", contour = FALSE, direction = -1, option = "B", legend.title = expression(hat(y))) + theme_classic(base_size = 20)

g2 <- deltaRug.fit %>%
  partial(pred.var = c("clarkEvans", "Site"), train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Clark-Evans Index Removed Stems (X-Y)", ylab = "Influence (\u0177)", contour = FALSE, direction = -1, option = "B", legend.title = expression(hat(y))) + theme_classic(base_size = 20)

g3 <- deltaRug.fit %>%
  partial(pred.var = c("deltaPropBACD", "Site"), train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "\u394 Prop. Plot Basal Area Codominant", ylab = "Influence (\u0177)", contour = FALSE, direction = -1, option = "B", legend.title = expression(hat(y))) + theme_classic(base_size = 20)

g4 <- deltaRug.fit %>%
  partial(pred.var = c("propRemStemsCDorD", "Site"), train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Prop. Removed Stems CD or D", ylab = "Influence (\u0177)", contour = FALSE, direction = -1, option = "B", legend.title = expression(hat(y))) + theme_classic(base_size = 20)

g5 <- deltaRug.fit %>%
  partial(pred.var = c("meanHeightCrownBase", "Site"), train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Mean Height-to-Crown of Removals", ylab = "Influence (\u0177)", contour = FALSE, direction = -1, option = "B", legend.title = expression(hat(y))) + theme_classic(base_size = 20)

g6 <- deltaRug.fit %>%
  partial(pred.var = c("propBAremD", "Site"), train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Prop. Plot Basal Area Removed Dominant", ylab = "Influence (\u0177)", contour = FALSE, direction = -1, option = "B", legend.title = expression(hat(y))) + theme_classic(base_size = 20)

g7 <- deltaRug.fit %>%
  partial(pred.var = c("deltaSdDBH", "Site"), train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "\u394 Standard Deviation DBH", ylab = "Influence (\u0177)", contour = FALSE, direction = -1, option = "B", legend.title = expression(hat(y))) + theme_classic(base_size = 20)

g8 <- deltaRug.fit %>%
  partial(pred.var = c("speciesCount", "Site"), train=Combined.train[,imp.vars], chull = FALSE, progress = TRUE) %>%
  autoplot(xlab = "Species Richness of Removals", ylab = "Influence (\u0177)", contour = FALSE, direction = -1, option = "B", legend.title = expression(hat(y))) + theme_classic(base_size = 20)

#grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, ncol = 3, nrow = 3) #use ggarrange if need to adjust plot margins

multiPdp2 <- ggarrange(g1, g3, g2, g4, g5, g7, g6, g8, labels = c("A", "B", "C", "D", "E", "F", "G", "H"), font.label = list(size = 24), ncol = 2, nrow = 4, align = "hv") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
annotate_figure(multiPdp2, fig.lab = "", fig.lab.pos = "top", fig.lab.face = "bold")

## Output raw partial dependency (as yhat) if needed (or alternative way of plotting 2-way PDP plots (see below) if plot = TRUE)
## Two-way partial dependency (with site) - combined effects of site and each important predictor
#works same whether specify type = "regression" or not; train = SHOULD MATCH WHATEVER USED TO TRAIN MODEL INDICATED HERE
#for partial - chull = FALSE by default so should match ggplot pdp plots
#took nearly an hour with full data/final rf params
imp.vars <- imp.vars[imp.vars != "Site"]
pd=list()
x11()
for(i in 1:length(imp.vars)){
  pd[[i]] <- partial(deltaRug.fit, pred.var=c("Site", imp.vars[i]), 
                     train=Combined.train[,c("Site", imp.vars)],
                     main=paste0("Partial dependency of ", imp.vars[i]),
                     plot = FALSE, ice=FALSE)
  pathToOutput <- "C:/Users/bca17002/Desktop/Ch2_FinalFiguresAndStatsOutputRFR/"
  fwrite(pd[[i]], file = paste0(pathToOutput, "DeltaRugosityTwoWayPD_", i, ".csv"))
}

##Cross-validation (n-fold here/possibly "nested")
#recommended to iterate over cross-validation multiple times to get the stabilized predictive accuracy (ie average or median) across multiple cv runs
#see https://cran.r-project.org/web/packages/spm/vignettes/spm.html
#     -all statistics based on differences between the predicted values for and the observed values of validation samples for cross-validation
#took ~20 mins at n = 10 with full training data and final rf params
n <- 10 # number of iterations, 60 to 100 is recommended
VEcv <- NULL
MSEcv <- NULL
RMSEcv <- NULL
Cohenf2cvMaybe <- NULL

#x11()
#par(mar=c(2,2,2,2))
par(mfrow=c(2,2))

for (i in 1:n) {
  rgcv1 <- rgcv(trainx = Combined.train[,imp.vars], trainy = Combined.train[,"DeltaRugosity"], cv.fold = 10, num.trees = b, mtry = m2, min.node.size = 5, predacc = "ALL") #min.node.size =5 is default for ranger and spm
  VEcv [i] <- rgcv1$vecv #variance explained (~R2) by predictive models based on cross-validation
  MSEcv [i] <- rgcv1$mse #prediction mse based on predicted vs observed values of validation samples from cross-validation
  RMSEcv [i] <- rgcv1$rmse #prediction rmse based on predicted vs observed values of validation samples from cross-validation
  Cohenf2cvMaybe [i] <- (rgcv1$vecv/100)/(1 - (rgcv1$vecv/100)) # Quasi-Cohen's f2 (global?) for predictive models based on (vecv and) cross-validation
}
plot(VEcv ~ c(1:n), xlab = "Iteration for RF Cross-Validation", ylab = "VEcv (%)")
points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
abline(h = median(VEcv), col = 'blue', lwd = 2)
median(VEcv) # REPORT

plot(MSEcv ~ c(1:n), xlab = "Iteration for RF Cross-Validation", ylab = "Prediction MSE")
points(cumsum(MSEcv) / c(1:n) ~ c(1:n), col = 2)
abline(h = median(MSEcv), col = 'blue', lwd = 2)
median(MSEcv) # REPORT

plot(RMSEcv ~ c(1:n), xlab = "Iteration for RF Cross-Validation", ylab = "Prediction RMSE")
points(cumsum(RMSEcv) / c(1:n) ~ c(1:n), col = 2)
abline(h = median(RMSEcv), col = 'blue', lwd = 2)
median(RMSEcv) # REPORT

plot(Cohenf2cvMaybe ~ c(1:n), xlab = "Iteration for RF Cross-Validation", ylab = "Cohen's f2 (from VEcv)")
points(cumsum(Cohenf2cvMaybe) / c(1:n) ~ c(1:n), col = 2)
abline(h = median(Cohenf2cvMaybe), col = 'blue', lwd = 2)
median(Cohenf2cvMaybe) # REPORT

########## End of analysis


