##### Figure making script - random disturbance simulations (Dissertation and Eco Apps versions)
#packages
library('ggplot2') # graphics
library('ggpubr')
library('ggridges')
library('ggpmisc')
library('car') # diagnostics; leveneTest - factors should be fully crossed, robust to slight non-normality
library('grid')
library('gridExtra')
library("stats") # diagnostics/summaries
library('dplyr')
library('data.table')
library('rlang')
library('viridis')

### Fig 2 - dotchart version (everything working correct including accurate/consistent names as of 4/6/23)
ImpDataDeltaRugosity <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/DeltaRugosity_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")

#dotChart version
x11()
#list labels lowest importance to highest
labs <- c("Species Richness of Removals", "Prop. Plot Basal Area Removed Dominant", "Δ Standard Deviation DBH", "Mean Height-to-Crown of Removals", "Prop. Removed Stems CD or D", "Clark-Evans Index Removed Stems (X-Y)", "Δ Prop. Plot Basal Area Codominant",  "Site", "Prop. Plot Basal Area Removed Suppressed")
#check to see if order correct
ImpDataDeltaRugosity <- cbind(ImpDataDeltaRugosity, labs)
#plot it
dotchart(ImpDataDeltaRugosity$Importance, labels = ImpDataDeltaRugosity$labs, lcolor = "darkgray",
         main="Scaled Variable Importance to Δ Rugosity", pch=19, pt.cex = 1.7)


### Segmented regressions - Fig 5, 6, 7 
#STREAMLINED VERSIONS APPEAR TO BE REASONABLE/CORRECT METHODOLOGY FROM VIGNETTE AND TESTING
#MAIN CONCERNS - PLAUSIBLY SERIAL CORRELATION 
#install.packages('segmented')
library('segmented')
library('car') #leveneTest
library('EnvStats')

##read in and format data - verified works correct for OGcorrected
Combined <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/MasterDatasets/Corrected/SimulationExperiment_ExperimentRuns_MasterDataset_OGcorrected.csv")

#Combined <- Combined[1:6000, ] #take first 2 sites for test only
Site <- Combined[ ,1]
Plot <- Combined[ ,2]
Severity <- Combined[ ,4]
predictorSet1 <- Combined[ ,8:80]
DeltaRugosity <- Combined[ ,106]
Combined <- cbind(Site, Plot, Severity, predictorSet1, DeltaRugosity) #uncomment if including site and plot
#Combined <- cbind(Severity, predictorSet1, DeltaRugosity) #if leaving out site and plot
sum(is.na(Combined)) # check for NAs
#Combined$Severity <- Combined$Severity * 100 #UNCLEAR IF NECESSARY

##severity (streamlined version) - appears to work correctly (no transform)
#fit lm
dRugMod <- lm(DeltaRugosity~Severity, data = Combined, na.action = na.exclude)
#test parametric assumptions - not ideal but generally ok
serialCorrelationTest(dRugMod, conf.level = 0.90, na.action = na.exclude) #TEST - appears serially corr
plot(cooks.distance(dRugMod))
par(mfrow = c(2, 2))
plot(dRugMod)
#transform if needed, refit lm if transformed
#minDeltaRugosity <- min(Combined$DeltaRugosity)
#Combined$DeltaRugosity <- Combined$DeltaRugosity + (abs(minDeltaRugosity)+1)
#dRugMod <- lm(log(DeltaRugosity)~Severity, data = Combined)
#determine optimal number breakpoints (max of 5 here)to specify for segmented lm and fit segmented lm with optimal number breakpoints
  #refit settings controlled by seg.control (defaults used as specified here - see them with "call")
dRugModSeg <- selgmented(dRugMod, alpha = 0.05, type = "aic", bonferroni = TRUE, Kmax = 5, refit = TRUE, return.fit = TRUE) #identify the number breakpoints in segmented relationship, then refit/repeat to get further optimized estimates and return the optimized model
#dRugModSeg <- selgmented(dRugMod, alpha = 0.05, type = "davies", bonferroni = TRUE, refit = TRUE, return.fit = TRUE) #sequentially hypothesis test for the number breakpoints (max of 2) in segmented relationship, then refit/repeat to get further optimized estimates and return the optimized model
#dRugModSeg <- segmented(dRugMod, npsi = 3) #(functionally equal alternative to above if specify return.fit = FALSE) fit segmented regression for single predictor and estimate 3 breakpoints
dRugModSeg[["call"]] #see what the model parameters were (includes breakpoint estimates that can verify against summary)
#plot result and get associated statistics to report
plot(DeltaRugosity~Severity, data = Combined, pch = 20, xlab = "Severity (Prop. Basal Area Removed)", ylab = "Δ Rugosity", col = "darkgray")
abline(a = 0, b = 0, col = "black", lty = "dashed", lwd = 3)
plot(dRugModSeg, col = "black", lwd = 4, conf.level = 0.95, add = TRUE)
ciSeg <- confint.segmented(dRugModSeg, method="delta", var.diff = TRUE) #use this as ci for estimated BREAKPOINT(S) (estimated breakpoint value also seems to match that of summary) - could potentially use to assess diffs among sites (in place of IQR)?
#ciLower <- ciSeg[ , 2]
#ciUpper <- ciSeg[ , 3]
#abline(v = ciLower, col = "black", lwd = 1)
#abline(v = ciUpper, col = "black", lwd = 1)
points.segmented(dRugModSeg, col = "black", pch = 1)
lines(dRugModSeg, bottom = FALSE, conf.level = 0.95, k = 2.5,  col = "black") #supposedly slightly less accurate than confint but display options better
slope(dRugModSeg, conf.level = 0.95, var.diff = TRUE) #use this as point estimate/ci/etc for SLOPES from segmented (ie either side of breakpoint) -appears ci based on quantiles
summary.segmented(dRugModSeg, var.diff=TRUE) #use this for value of estimated BREAKPOINT(S) (var.diff appears to allow different variance structures by segment)
#can use if one or more break point (lm), test for non-constant regression param (ie change in slope) - possibly more like ANOVA if multiple breakpoints?
#davies.test(dRugMod, seg.Z = ~Severity, values = c(0.420, 0.571, 0.747)) #if more than one change point davies.test preferred (here appears to test the specified breakpoint values, which overwrite "k=" if included as indicated by n.points =3 in output)
#davies.test(dRugMod, seg.Z = ~Severity, type = "lrt", values = c(0.420, 0.571, 0.747)) #identical result to above line
mtext(text = expression(paste("Adj. Multiple R" ^ 2, " = 0.08")), adj = 0.03, line = -1.5, cex = 0.9)
#mtext(text = expression(paste(italic("P"), "< 0.01")), adj = 0, line = -2, cex = 0.9)

sev80orAbove <- subset(Combined, Severity > 0.80) #check number obs above 80% sev


##vertical
#par(mfrow = c(1, 2)) #run if want to plot bacd and meanHeightCrown
#deltaPropBACD - appears to work correctly (no transform)
#fit lm
dRugMod <- lm(DeltaRugosity~deltaPropBACD, data = Combined)
#test parametric assumptions - generally acceptable but not ideal
serialCorrelationTest(dRugMod, conf.level = 0.90, na.action = na.exclude) #TEST - appears serially corr
plot(cooks.distance(dRugMod))
par(mfrow = c(2, 2))
plot(dRugMod)
#transform if needed, refit lm if transformed
#minDeltaRugosity <- min(Combined$DeltaRugosity)
#Combined$DeltaRugosity <- Combined$DeltaRugosity + (abs(minDeltaRugosity)+1)
#dRugMod <- lm(log(DeltaRugosity)~deltaPropBACD, data = Combined)
#determine optimal number breakpoints (max of 5 here)to specify for segmented lm (experimental function) and fit segmented lm with optimal number breakpoints
  #refit settings controlled by seg.control (defaults used as specified here - see them with "call")
dRugModSeg <- selgmented(dRugMod, alpha = 0.05, type = "aic", bonferroni = TRUE, Kmax = 5, refit = TRUE, return.fit = TRUE) #identify the number breakpoints in segmented relationship, then refit/repeat to get further optimized estimates and return the optimized model
#dRugModSeg <- segmented(dRugMod, npsi = 4) #(functionally equal alternative to above if specify return.fit = FALSE) fit segmented regression for single predictor and estimate 4 breakpoints
dRugModSeg[["call"]] #see what the model parameters were (includes breakpoint estimates that can verify against summary)
#plot result and get associated statistics to report
plot(DeltaRugosity~deltaPropBACD, data = Combined, pch = 20, xlab = "Δ Prop. Plot Basal Area Codominant", ylab = "Δ Rugosity", col = "darkgray")
abline(a = 0, b = 0, col = "black", lty = "dashed", lwd = 3)
plot(dRugModSeg, col = "black", lwd = 4, conf.level = 0.95, add = TRUE)
ciSeg <- confint.segmented(dRugModSeg, method="delta", var.diff = TRUE) #use this as ci for estimated BREAKPOINT(S) (estimated breakpoint value also seems to match that of summary) - could potentially use to assess diffs among sites (in place of IQR)?
#ciLower <- ciSeg[ , 2]
#ciUpper <- ciSeg[ , 3]
#abline(v = ciLower, col = "black", lwd = 1)
#abline(v = ciUpper, col = "black", lwd = 1)
points.segmented(dRugModSeg, col = "black", pch = 1)
lines(dRugModSeg, bottom = FALSE, conf.level = 0.95, k = 2.5, col = "black") #supposedly slightly less accurate than confint but display options better
slope(dRugModSeg, conf.level = 0.95, var.diff = TRUE) #use this as point estimate/ci/etc for SLOPES from segmented (ie either side of breakpoint) -appears ci based on quantiles
summary.segmented(dRugModSeg, var.diff=TRUE) #use this for value of estimated BREAKPOINT(S) (var.diff appears to allow different variance structures by segment)
#can use if one or more break point (lm), test for non-constant regression param (ie change in slope) - possibly more like ANOVA if multiple breakpoints?
#davies.test(dRugMod, seg.Z = ~deltaPropBACD, values = c(-0.567, -0.515, -0.350, -0.135)) #if more than one change point davies.test preferred (here appears to test the specified breakpoint values, which overwrite "k=" if included as indicated by n.points =3 in output)
#davies.test(dRugMod, seg.Z = ~Severity, type = "lrt", values = c(0.420, 0.571, 0.747)) #identical result to above line
#title(outer=FALSE, adj=0.05, main="A", cex.main=2, col="black", font=2, line=-2) #run if want to plot bacd and meanHeightCrown
mtext(text = expression(paste("Adj. Multiple R" ^ 2, " = 0.10")), adj = 0.95, line = -13.5, cex = 0.9)
#mtext(text = expression(paste(italic("P"), "> 0.05")), adj = 0, line = -2, cex = 0.9)


##horizontal
#clarkEvans - appears to work correctly (no transform)
#fit lm
dRugMod <- lm(DeltaRugosity~clarkEvans, data = Combined)
#test parametric assumptions - generally acceptable but not ideal
serialCorrelationTest(dRugMod, conf.level = 0.90, na.action = na.exclude) #TEST - appears serially corr
plot(cooks.distance(dRugMod))
par(mfrow = c(2, 2))
plot(dRugMod)
#transform if needed, refit lm if transformed
#minDeltaRugosity <- min(Combined$DeltaRugosity)
#Combined$DeltaRugosity <- Combined$DeltaRugosity + (abs(minDeltaRugosity)+1)
#dRugMod <- lm(log(DeltaRugosity)~clarkEvans, data = Combined)
#determine optimal number breakpoints (max of 5 here)to specify for segmented lm (experimental function) and fit segmented lm with optimal number breakpoints
  #refit settings controlled by seg.control (defaults used as specified here - see them with "call")
dRugModSeg <- selgmented(dRugMod, alpha = 0.05, type = "aic", Kmax = 5, bonferroni = TRUE, refit = TRUE, return.fit = TRUE) #identify the number breakpoints in segmented relationship, then refit/repeat to get further optimized estimates and return the optimized model
#dRugModSeg <- segmented(dRugMod, npsi = 3) #(functionally equal alternative to above if specify return.fit = FALSE) fit segmented regression for single predictor and estimate 3 breakpoints
dRugModSeg[["call"]] #see what the model parameters were (includes breakpoint estimates that can verify against summary)
#plot result and get associated statistics to report
plot(DeltaRugosity~clarkEvans, data = Combined, pch = 20, xlab = "Clark-Evans Index Removed Stems (X-Y)", ylab = "Δ Rugosity", col = "darkgray")
abline(a = 0, b = 0, col = "black", lty = "dashed", lwd = 3)
plot(dRugModSeg, col = "black", lwd = 4, conf.level = 0.95, add = TRUE)
ciSeg <- confint.segmented(dRugModSeg, method="delta", var.diff = TRUE) #use this as ci for estimated BREAKPOINT(S) (estimated breakpoint value also seems to match that of summary) - could potentially use to assess diffs among sites (in place of IQR)?
#ciLower <- ciSeg[ , 2]
#ciUpper <- ciSeg[ , 3]
#abline(v = ciLower, col = "black", lwd = 1)
#abline(v = ciUpper, col = "black", lwd = 1)
points.segmented(dRugModSeg, col = "black", pch = 1)
lines(dRugModSeg, bottom = FALSE, conf.level = 0.95, k = 3.5, col = "black") #supposedly slightly less accurate than confint but display options better
slope(dRugModSeg, conf.level = 0.95, var.diff = TRUE) #use this as point estimate/ci/etc for SLOPES from segmented (ie either side of breakpoint) -appears ci based on quantiles
summary.segmented(dRugModSeg, var.diff=TRUE) #use this for value of estimated BREAKPOINT(S) (var.diff appears to allow different variance structures by segment)
#can use if one or more break point (lm), test for non-constant regression param (ie change in slope) - possibly more like ANOVA if multiple breakpoints?
#davies.test(dRugMod, seg.Z = ~clarkEvans, values = c(0.861, 0.938, 1.038)) #if more than one change point davies.test preferred (here appears to test the specified breakpoint values, which overwrite "k=" if included as indicated by n.points =3 in output)
#davies.test(dRugMod, seg.Z = ~Severity, type = "lrt", values = c(0.420, 0.571, 0.747)) #identical result to above line
mtext(text = expression(paste("Adj. Multiple R" ^ 2, " = 0.25")), adj = 0.05, line = -1.5, cex = 0.9)
#mtext(text = expression(paste(italic("P"), "> 0.05")), adj = 0, line = -2, cex = 0.9)


###create correlation matrix among all non-forestr predictors and selected forestr responses
#UPDATED/PROOFED 4/19/23 - works great!
#install.packages("corrplot")
#install.packages("Hmisc")
library('corrplot') # draw correlogram from cor() or rcorr()
library('Hmisc') # correlations with significance
library('tidyr')

##read in the OGcorrected dataset
Combined <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/MasterDatasets/Corrected/SimulationExperiment_ExperimentRuns_MasterDataset_OGcorrected.csv")

#rename response variables so more intuitive
colnames(Combined)[101] <- "DeltaMeanVAI"
colnames(Combined)[94] <- "DeltaMeanLeafHeight"
colnames(Combined)[105] <- "DeltaMeanSDLeafHeight"
colnames(Combined)[103] <- "DeltaPorosity"
colnames(Combined)[106] <- "DeltaRugosity"

##read in data with final predictors each response/create list of all unique predictors (with updated names)
ImpDataDeltaRugosity <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/DeltaRugosity_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")
labs <- c("Species Richness of Removals", "Prop. Plot Basal Area Removed Dominant", "Δ Standard Deviation DBH", "Mean Height-to-Crown of Removals", "Prop. Removed Stems CD or D", "Clark-Evans Index Removed Stems (X-Y)", "Δ Prop. Plot Basal Area Codominant",  "Site", "Prop. Plot Basal Area Removed Suppressed")
ImpDataDeltaRugosity <- cbind(ImpDataDeltaRugosity, labs)

ImpDataDeltaMeanVAI <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaMeanVAI/DeltaMeanVAI_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")
labs <- c("Mean Height-to-Crown of Removals", "Δ Prop. Plot Basal Area Codominant", "Δ Prop. Plot Basal Area Overtopped", "Site")
ImpDataDeltaMeanVAI <- cbind(ImpDataDeltaMeanVAI, labs)

ImpDataDeltaMeanLeafHeight <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaMeanHeight/DeltaMeanLeafHeight_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")
labs <- c("Conifer Fraction of Removals", "Clark-Evans Index Removed Stems (X-Y)", "Δ Prop. Plot Basal Area Intermediate", "Δ Range DBH", "Shannon Diversity Spp of Removals", "Δ Fisher Alpha DBH", "Mean Height-to-Crown of Removals", "Δ Prop. Plot Basal Area Overtopped", "Site", "Prop. Plot Basal Area Removed Dominant")
ImpDataDeltaMeanLeafHeight <- cbind(ImpDataDeltaMeanLeafHeight, labs)

ImpDataDeltaMeanSDLeafHeight <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaMeanStd/DeltaMeanSDLeafHeight_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")
labs <- c("Prop. Plot Basal Area Removed Dominant", "Δ Conifer Fraction", "Prop. Plot Basal Area Removed Snag", "Mean Height-to-Crown of Removals", "Shannon Diversity Spp of Removals", "Δ Prop. Plot Basal Area Overtopped", "Site","Prop. Plot Basal Area Removed Suppressed", "Δ Prop. Plot Basal Area Intermediate", "Conifer Fraction of Removals", "Δ Prop. Plot Basal Area Codominant")
ImpDataDeltaMeanSDLeafHeight <- cbind(ImpDataDeltaMeanSDLeafHeight, labs)

ImpDataDeltaPorosity <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaPorosity/DeltaPorosity_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")
labs <- c("Site", "Clark-Evans Index Removed Stems (X-Y)", "Mean Height-to-Crown of Removals", "Prop. Plot Basal Area Removed Dominant", "Conifer Fraction of Removals", "Species Richness of Removals", "Δ Fisher Alpha Mean Crown Diameter", "Δ Contagion Index Stems (X-Y)", "Δ Fisher Alpha DBH", "Fisher Alpha DBH of Removals", "Prop. Plot Basal Area Removed Suppressed", "Δ Prop. Plot Basal Area Overtopped", "Δ Prop. Plot Basal Area Intermediate", "Δ Prop. Plot Basal Area Codominant")
ImpDataDeltaPorosity <- cbind(ImpDataDeltaPorosity, labs)

finalPredsDF <- rbind.data.frame(ImpDataDeltaRugosity, ImpDataDeltaMeanVAI, ImpDataDeltaMeanLeafHeight, ImpDataDeltaMeanSDLeafHeight, ImpDataDeltaPorosity)
uniqueFinalPreds <- as.character(unique(finalPredsDF$Parameter))
nrow(subset(finalPredsDF, Parameter %in% uniqueFinalPreds)) #checks if all final predictors were identified

#remove "Site" and "Plot" from uniqueFinalPredictors if included as only non-numeric can be used for correlations
uniqueFinalPreds <- uniqueFinalPreds[-which(uniqueFinalPreds %in% c("Site", "Plot"))]
uniqueFinalPreds #check to verify worked

#add the updated names for response variables to list of predictor names
#given want to test for correlations between predictors and responses too
uniqueCorrTestVars <-c(uniqueFinalPreds, "DeltaMeanVAI", "DeltaMeanLeafHeight", "DeltaMeanSDLeafHeight", "DeltaPorosity", "DeltaRugosity")
uniqueCorrTestVars #check to verify worked
#uniqueCorrTestVars <- uniqueFinalPreds #uncomment if don't want to include response variables

sum(is.na(Combined[, uniqueCorrTestVars])) # check for NAs

##compute correlations and plot correlogram full dataset with only vars of interest (no sig levels) - WORKS but ratio goofy if not using jpeg() or similar
res <- cor(Combined[ ,uniqueCorrTestVars], method = "spearman", use = "complete.obs")

#(optionally) rename variable names if needed
#verified works with data used here (6/8/23)
varNames <- as.character(unique(finalPredsDF$labs))
nrow(subset(finalPredsDF, labs %in% varNames)) #checks if all final predictors were identified
varNames <- varNames[-which(varNames %in% c("Site", "Plot"))] #remove site and/or plot as non-numeric
varNames #check to verify worked
#comment the below line out if did not include response vars above
varNames <-c(varNames, "Δ MeanVAI", "Δ MeanLeafHeight", "Δ MeanSDLeafHeight", "Δ Porosity", "Δ Rugosity") # add response vars

colnames(res) <- varNames
rownames(res) <- varNames

#run these lines in order to plot it - outputs to last folder exported to
filetag <- "corrplot_resultTest.jpg"

jpeg(filetag, height = 1500, width = 1500, quality = 100)

corrplot(res, type = "upper", order = "hclust", cl.cex = 2, tl.cex = 1.5, tl.col = "black", tl.srt = 45)
dev.off()
#test to change aspect ratio (line below)- did not help...
#corrplot(res, method = "square", type = "full", order = "hclust", tl.cex = 0.5, tl.col = "black", win.asp = 0.5)

#OR (to use below must clear environment and re-run everything from above down to "res")

##compute correlations and plot correlogram full dataset with only vars of interest (with significance) - works correctly!
res2 <- rcorr(as.matrix(Combined[ ,uniqueCorrTestVars]), type = "spearman") #from Hmisc - creates NAs in $P giving unequal #obs P vs r

#isolate p-values matrix and replace NA values
P <- res2$P
is.matrix(P) # check if matrix
ncol(P)
nrow(P)
sum(is.na(P)) # check for NAs
P[is.na(P)] <- 0 #can replace NAs with 0 as only NA if cor = 1 (otherwise error)
sum(is.na(P)) # check for NAs

#isolate correlations matrix
r <- res2$r
is.matrix(r) # check if matrix
ncol(r)
nrow(r)
sum(is.na(r)) # check for NAs

#get min and max absolute spearman rho for each predictor for respective response
#DeltaRugosity - works correctly
allCorrs <- as.data.frame(r)
rugCorrs <- subset(allCorrs, row.names(allCorrs) %in% ImpDataDeltaRugosity$Parameter)
min(abs(rugCorrs$DeltaRugosity))
max(abs(rugCorrs$DeltaRugosity))

#(optionally) rename variable names if needed and (optionally) add the updated names for response variables to list of renamed predictors
#verified works with data used here (6/8/23)
varNames <- as.character(unique(finalPredsDF$labs))
nrow(subset(finalPredsDF, labs %in% varNames)) #checks if all final predictors were identified
varNames <- varNames[-which(varNames %in% c("Site", "Plot"))] #remove site and/or plot as non-numeric
varNames #check to verify worked
varNames <-c(varNames, "Δ MeanVAI", "Δ MeanLeafHeight", "Δ MeanSDLeafHeight", "Δ Porosity", "Δ Rugosity") # add response vars

colnames(P) <- varNames
rownames(P) <- varNames

colnames(r) <- varNames
rownames(r)<- varNames

#plot it
filetag <- "corrplot_resultTestSig.jpg"

jpeg(filetag, height = 1500, width = 1500, quality = 100)

corrplot(corr = r, type = "upper", order = "original", cl.cex = 2, tl.cex = 1.5, tl.col = "black", tl.srt = 45, p.mat = P, sig.level = 0.05, insig = "p-value", na.label = "?", is.corr = TRUE)
dev.off()


### Composition summary figure for appendix (color scales etc best I could do)
#call CSCME field data as easier than laoding all plots sep
CSCME <- read.csv(file="C:/Users/balve/Desktop/SideProjects/CSCME/Data/Field_Data/CSCME_PreCutData_Summer2019_FinalWithTreatmentNotes.csv")

#mutate the Site column so sites have numbers too - verified
CSCME <- CSCME %>% 
  mutate(Site = case_when(
    Site == "NorthEagleville" ~ "(1) North Eagleville",
    Site == "Moss" ~ "(2) Moss",
    Site == "FentonWest" ~ "(3) Fenton West",
    Site == "FentonEast" ~ "(4) Fenton East"))

CSCME <- CSCME %>% 
  mutate(CommonName = case_when(
    CommonName == "Bigtooth Aspen" ~ "Uncommon/Other",
    CommonName == "Flowering Dogwood" ~ "Uncommon/Other",
    CommonName == "Black Cherry" ~ "Uncommon/Other",
    CommonName == "American Elm" ~ "Uncommon/Other",
    CommonName == "Norway Maple" ~ "Uncommon/Other",
    CommonName == "Unidentified" ~ "Uncommon/Other",
    CommonName == "Sassafras" ~ "Uncommon/Other",
    CommonName == "Witch Hazel" ~ "Uncommon/Other",
    CommonName == "Common Hawthorne" ~ "Uncommon/Other",
    CommonName == "Sugar Maple" ~ "Acer saccharum",
    CommonName == "Black Birch" ~ "Betula lenta",
    CommonName == "Eastern Hemlock" ~ "Tsuga canadensis",
    CommonName == "Yellow Birch" ~ "Betula alleghaniensis",
    CommonName == "Shagbark Hickory" ~ "Carya ovata",
    CommonName == "Mockernut Hickory" ~ "Carya tomentosa",
    CommonName == "Bitternut Hickory" ~ "Carya cordiformis",
    CommonName == "Pignut Hickory" ~ "Carya glabra",
    CommonName == "Red Maple" ~ "Acer rubrum",
    CommonName == "American Hornbeam" ~ "Carpinus caroliniana",
    CommonName == "American Beech" ~ "Fagus grandifolia",
    CommonName == "Scarlet Oak" ~ "Quercus coccinea",
    CommonName == "Black Oak" ~ "Quercus velutina",
    CommonName == "White Oak" ~ "Quercus alba",
    CommonName == "Eastern Hophornbeam" ~ "Ostrya virginiana",
    CommonName == "White Ash" ~ "Fraxinus americana",
    CommonName == "Red Oak" ~ "Quercus rubra"
  ))

colnames(CSCME)[5] <- "Species"
unique(CSCME$Species) # should be 18

ggplot(CSCME, aes(x = Site, y = BasalArea_m2, fill = Species)) +
  geom_col(position = "fill") + labs(y = "Relative Basal Area") + scale_fill_viridis(discrete = TRUE, option = "H") +
  guides(fill = guide_legend(label.theme = element_text(face = "italic"))) + theme_classic(base_size = 20)


### Calculate stand characteristic summaries for appendix table
#call CSCME field data as easier than laoding all plots sep
CSCME <- read.csv(file="C:/Users/balve/Desktop/SideProjects/CSCME/Data/Field_Data/CSCME_PreCutData_Summer2019_FinalWithTreatmentNotes.csv")

#subset by site
NorthEagleville <- CSCME %>% filter(Site == "NorthEagleville")
Moss <- CSCME %>% filter(Site == "Moss")
FentonWest <- CSCME %>% filter(Site == "FentonWest")
FentonEast <- CSCME %>% filter(Site == "FentonEast")

#stand structure
stemsPerHa <- nrow(FentonEast)/0.75 #stems per ha
baPerHa <- sum(FentonEast$BasalArea_m2)/0.75 #basal area per ha (m2)
meanDBH <- mean(FentonEast$DBH_cm)

#composition
sppRichness <- n_distinct(FentonEast$CommonName) #spp richness
n_distinct(FentonEast$Species_Code) #qc check spp richness
coniferBA <- FentonEast %>% filter(CommonName == "Eastern Hemlock")
coniferBA <- sum(coniferBA$BasalArea_m2)
coniferFraction <- coniferBA/sum(FentonEast$BasalArea_m2)
####

### Final updated Figure S3 2-way (each final predictor by site) pdp figure for Eco Apps pub (ch2) - works correctly (8/10/23)/matches dissertation figure trends
#COULD LIKELY USE THIS SYNTAX/METHOD TO REVISE THE 1-WAY VERSION TOO IF NEEDED

#propBAremSuppressed - pd_3
pdppPropBAremSupp <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/RawPartialDependencyOutput_TwoWayWithSite/DeltaRugosityTwoWayPD_3.csv")

ggPropBAremSupp <- ggplot(pdppPropBAremSupp, aes(propBAremSuppressed, y=yhat, col=Site))+
  geom_line(lwd=1.2)+xlab("Prop. Plot Basal Area Removed Suppressed")+
  scale_y_continuous(limits = c(0, 3.5)) + theme_classic()+ylab("Influence (ŷ)")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "none",
        legend.title=element_text(size=14),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=14))

#deltaPropBACD - pd_7
pdppDeltaBACD <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/RawPartialDependencyOutput_TwoWayWithSite/DeltaRugosityTwoWayPD_7.csv")

ggDeltaBACD <- ggplot(pdppDeltaBACD, aes(deltaPropBACD, y=yhat, col=Site))+
  geom_line(lwd=1.2)+xlab("Δ Prop. Plot Basal Area Codominant")+
  scale_y_continuous(limits = c(0, 3.5)) + theme_classic()+ylab("Influence (ŷ)")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "none",
        legend.title=element_text(size=14),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=14))

#clarkEvans - pd_6
pdppClarkEvans <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/RawPartialDependencyOutput_TwoWayWithSite/DeltaRugosityTwoWayPD_6.csv")

ggClarkEvans <- ggplot(pdppClarkEvans, aes(clarkEvans, y=yhat, col=Site))+
  geom_line(lwd=1.2)+xlab("Clark-Evans Index Removed Stems (X-Y)")+
  scale_y_continuous(limits = c(0, 3.5)) + theme_classic()+ylab("Influence (ŷ)")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "none",
        legend.title=element_text(size=14),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=14))

#propRemStemsCDorD - pd_2
pdppPropRemStemsCDorD <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/RawPartialDependencyOutput_TwoWayWithSite/DeltaRugosityTwoWayPD_2.csv")

ggPropRemStemsCDorD <- ggplot(pdppPropRemStemsCDorD, aes(propRemStemsCDorD, y=yhat, col=Site))+
  geom_line(lwd=1.2)+xlab("Prop. Removed Stems CD or D")+
  scale_y_continuous(limits = c(0, 3.5)) + theme_classic()+ylab("Influence (ŷ)")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "none",
        legend.title=element_text(size=14),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=14))

#meanHeightCrownBase - pd_5
pdppMeanHtCrown <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/RawPartialDependencyOutput_TwoWayWithSite/DeltaRugosityTwoWayPD_5.csv")

ggMeanHtCrown <- ggplot(pdppMeanHtCrown, aes(meanHeightCrownBase, y=yhat, col=Site))+
  geom_line(lwd=1.2)+xlab("Mean Height-To-Crown of Removals (m)")+
  scale_y_continuous(limits = c(0, 3.5)) + theme_classic()+ylab("Influence (ŷ)")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "none",
        legend.title=element_text(size=14),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=14))

#deltaSdDBH -pd_8
pdppDeltaSdDBH <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/RawPartialDependencyOutput_TwoWayWithSite/DeltaRugosityTwoWayPD_8.csv")

ggDeltaSdDBH <- ggplot(pdppDeltaSdDBH, aes(deltaSdDBH, y=yhat, col=Site))+
  geom_line(lwd=1.2)+xlab("Δ Standard Deviation DBH")+
  scale_y_continuous(limits = c(0, 3.5)) + theme_classic()+ylab("Influence (ŷ)")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "none",
        legend.title=element_text(size=14),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=14))

#propBAremD - pd_4
pdppPropBAremD <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/RawPartialDependencyOutput_TwoWayWithSite/DeltaRugosityTwoWayPD_4.csv")

ggPropBAremD <- ggplot(pdppPropBAremD, aes(propBAremD, y=yhat, col=Site))+
  geom_line(lwd=1.2)+xlab("Prop. Plot Basal Area Removed Dominant")+
  scale_y_continuous(limits = c(0, 3.5)) + theme_classic()+ylab("Influence (ŷ)")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "none",
        legend.title=element_text(size=14),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=14))

#speciesCount - pd_1
pdppSpeciesCount <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/RawPartialDependencyOutput_TwoWayWithSite/DeltaRugosityTwoWayPD_1.csv")

ggSpeciesCount <- ggplot(pdppSpeciesCount, aes(speciesCount, y=yhat, col=Site))+
  geom_line(lwd=1.2)+xlab("Species Richness of Removals")+
  scale_y_continuous(limits = c(0, 3.5))+theme_classic()+ylab("Influence (ŷ)")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "none",
        legend.title=element_text(size=14),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=14))

#plot multipanel
multiPlot <- ggarrange(ggPropBAremSupp, ggDeltaBACD, ggClarkEvans, ggPropRemStemsCDorD, ggMeanHtCrown, ggDeltaSdDBH, ggPropBAremD, ggSpeciesCount, ncol=2, nrow=4, labels = c("A", "B", "C", "D", "E", "F", "G", "H", font.label = list(size = 26)), legend = "bottom", common.legend = TRUE)
annotate_figure(multiPlot, fig.lab="", fig.lab.pos="top", top = textGrob("Partial Dependence of Δ Rugosity on Top Predictors", gp=gpar(fontsize=26)))


###Figure 3 - ggdotchart version - verified to match trends in dissertation version fig (8/11/23)
#COULD USE THIS SYNTAX/METHOD TO REVISE THE DELTARUGOSITY DOTCHART TOO IF NEEDED
#read in non-rugosity responses, add better labels and check to see if match order parameters (verified working correctly as of 8/10/23)
ImpDataDeltaMeanVAI <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaMeanVAI/DeltaMeanVAI_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")
labs <- c("Mean Height-to-Crown of Removals*", "Δ Prop. Plot Basal Area Codominant*", "Δ Prop. Plot Basal Area Overtopped", "Site*")
ImpDataDeltaMeanVAI <- cbind(ImpDataDeltaMeanVAI, labs)

ImpDataDeltaMeanLeafHeight <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaMeanHeight/DeltaMeanLeafHeight_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")
labs <- c("Conifer Fraction of Removals", "Clark-Evans Index Removed Stems (X-Y)*", "Δ Prop. Plot Basal Area Intermediate", "Δ Range DBH", "Shannon Diversity Spp of Removals", "Δ Fisher Alpha DBH", "Mean Height-to-Crown of Removals*", "Δ Prop. Plot Basal Area Overtopped", "Site*", "Prop. Plot Basal Area Removed Dominant*")
ImpDataDeltaMeanLeafHeight <- cbind(ImpDataDeltaMeanLeafHeight, labs)

ImpDataDeltaMeanSDLeafHeight <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaMeanStd/DeltaMeanSDLeafHeight_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")
labs <- c("Prop. Plot Basal Area Removed Dominant*", "Δ Conifer Fraction", "Prop. Plot Basal Area Removed Snag", "Mean Height-to-Crown of Removals*", "Shannon Diversity Spp of Removals", "Δ Prop. Plot Basal Area Overtopped", "Site*","Prop. Plot Basal Area Removed Suppressed*", "Δ Prop. Plot Basal Area Intermediate", "Conifer Fraction of Removals", "Δ Prop. Plot Basal Area Codominant*")
ImpDataDeltaMeanSDLeafHeight <- cbind(ImpDataDeltaMeanSDLeafHeight, labs)

ImpDataDeltaPorosity <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaPorosity/DeltaPorosity_scaledImportanceFinalPredictorsOnly_UpdatedNames.csv")
labs <- c("Site*", "Clark-Evans Index Removed Stems (X-Y)*", "Mean Height-to-Crown of Removals*", "Prop. Plot Basal Area Removed Dominant*", "Conifer Fraction of Removals", "Species Richness of Removals*", "Δ Fisher Alpha Mean Crown Diameter", "Δ Contagion Index Stems (X-Y)", "Δ Fisher Alpha DBH", "Fisher Alpha DBH of Removals", "Prop. Plot Basal Area Removed Suppressed*", "Δ Prop. Plot Basal Area Overtopped", "Δ Prop. Plot Basal Area Intermediate", "Δ Prop. Plot Basal Area Codominant*")
ImpDataDeltaPorosity <- cbind(ImpDataDeltaPorosity, labs)

#create the individual panel for each response
#DeltaMeanVAI
DeltaMeanVAI <- ggdotchart(ImpDataDeltaMeanVAI, x = "labs", y = "Importance", xlab = FALSE, ylab = FALSE, 
                      dot.size = 4, add = "segment", sorting = "descending", rotate = TRUE, title = "")
DeltaMeanVAI <- ggpar(DeltaMeanVAI, ylim = c(0, 1), yticks.by = 0.2, orientation = "horizontal")

DeltaMeanLeafHeight <- ggdotchart(ImpDataDeltaMeanLeafHeight, x = "labs", y = "Importance", xlab = FALSE, ylab = FALSE, 
                      dot.size = 4, add = "segment", sorting = "descending", rotate = TRUE, title = "")
DeltaMeanLeafHeight <- ggpar(DeltaMeanLeafHeight, ylim = c(0, 1), yticks.by = 0.2, orientation = "horizontal")

DeltaMeanSDLeafHeight <- ggdotchart(ImpDataDeltaMeanSDLeafHeight, x = "labs", y = "Importance", xlab = FALSE, ylab = FALSE, 
                                  dot.size = 4, add = "segment", sorting = "descending", rotate = TRUE, title = "")
DeltaMeanSDLeafHeight <- ggpar(DeltaMeanSDLeafHeight, ylim = c(0, 1), yticks.by = 0.2, orientation = "horizontal")

DeltaPorosity <- ggdotchart(ImpDataDeltaPorosity, x = "labs", y = "Importance", xlab = FALSE, ylab = FALSE, 
                           dot.size = 4, add = "segment", sorting = "descending", rotate = TRUE, title = "")
DeltaPorosity <- ggpar(DeltaPorosity, ylim = c(0, 1), yticks.by = 0.2, orientation = "horizontal")

#plot multipanel
multiDotChartPlot <- ggarrange(DeltaMeanVAI, DeltaMeanLeafHeight, DeltaMeanSDLeafHeight, DeltaPorosity, ncol=2, nrow=2, align = "hv", labels = c("A", "B", "C", "D"), font.label = list(size = 24))
annotate_figure(multiDotChartPlot, bottom = text_grob("Scaled Variable Importance", size = 16, vjust = 0.2, hjust = -0.145))


###Figure 8 - Sankey equivalent - predictor importance to deltaRugosity - relative density at diff sev ranges - works correctly
ImpDataDeltaRugSankey <- read.csv(file="C:/Users/balve/Desktop/Dissertation/Chapter2_RandomDisturbanceSimulations/FinalResults_RFR/ImportanceFromFinalRFR/DeltaRugosity/BySeverityRange/DeltaRugosity_scaledImportanceFinalPredictorsOnly_AllSeverities_UpdatedNames.csv")

#mutate the Predictor column so "delta" is "Δ" and severity column so indicates % or proportion
ImpDataDeltaRugSankey <- ImpDataDeltaRugSankey %>% 
  mutate(Predictor = case_when(
    Parameter == "deltaPropBACD" ~ "Δ Prop. Plot Basal Area Codominant",
    Parameter == "speciesCount" ~ "Species Richness of Removals",
    Parameter == "meanHeightCrownBase" ~ "Mean Height-to-Crown of Removals",
    Parameter == "propBAremD" ~ "Prop. Plot Basal Area Removed Dominant",
    Parameter == "deltaSdDBH" ~ "Δ Standard Deviation DBH",
    Parameter == "clarkEvans" ~ "Clark-Evans Index Removed Stems (X-Y)",
    Parameter == "propRemStemsCDorD" ~ "Prop. Removed Stems CD or D",
    Parameter == "propBAremSuppressed" ~ "Prop. Plot Basal Area Removed Suppressed",
    Parameter == "Site" ~ "Site"))

ImpDataDeltaRugSankey <- ImpDataDeltaRugSankey %>% 
  mutate(Severity = case_when(
    Severity == "20-42" ~ "0.20 - 0.42",
    Severity == "42-57" ~ "0.42 - 0.57",
    Severity == "57-75" ~ "0.57 - 0.75",
    Severity == "75-85" ~ "0.75 - 0.85"))

level_order <- c("Prop. Plot Basal Area Removed Suppressed", "Δ Prop. Plot Basal Area Codominant", "Prop. Removed Stems CD or D", "Mean Height-to-Crown of Removals", "Prop. Plot Basal Area Removed Dominant", "Site", "Clark-Evans Index Removed Stems (X-Y)", "Δ Standard Deviation DBH", "Species Richness of Removals")

ggplot(ImpDataDeltaRugSankey, aes(x = Severity, y = Importance, fill = factor(Predictor, level = level_order))) +
  geom_col(position = "fill") + labs(x = "Severity (Prop. Basal Area Removed)", y = "Relative Scaled Importance to Δ Rugosity", fill = "Predictor") + scale_fill_viridis(discrete = TRUE, option = "H") +
  theme_classic(base_size = 20)


