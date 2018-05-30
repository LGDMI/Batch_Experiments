#### Statistics - Batch Experiments ####

### This script will perform som statistical tests for the cell count and the VFA-data of the Batch-experiments
### It will display means as well as standard deviations
### Remark: for the power analysis and the ANOVA-tests, we should actually assume that we know that there are no differences between the donors
### in their response to the Emulsifiers. We know for a fact that they DO differ, but we ignore this fact for a moment and take their variance.

# install packages & load libraries #
install.packages("ggplot2")
install.packages('dplyr')
install.packages('tidyr')
install.packages('Hmisc')
install.packages('pwr')


library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(pwr)
library(readxl)

set.seed(2)


#-----------------------------------------------------------------------------------------------------------------#

###CELLCOUNTS
###Statistics will be performed on Live/dead ratio's and/or liver percentages

### read-in data from raw data excel ###
#data.cells <- read.csv2('SGPI_All_Donors.csv', header= TRUE)
data.cells <- readxl::read_xlsx("SGPI_All_Donors.xlsx",
                                sheet = "Overview_Cellcounts")
#data.cells<- na.omit(data.cells)


levene.test(y, group, location=c("median", "mean", "trim.mean"), trim.alpha=0.25,
            bootstrap = FALSE, num.bootstrap=1000, kruskal.test=FALSE, 
            correction.method=c("none","correction.factor","zero.removal","zero.correction"))



#-----------------------------------------------------------------------------------------------------------------#

### Effect size calculation ###
### Per ANOVA you want to do, the effect size needs to be calculated in order 
### to calc the required sample size.

#Mean.Ratio = mean(data.cells[["Ratio"]])
#Var.Ratio = var(data.cells[["Ratio"]])



### CMC
  Grouped.data  <- group_by(data.cells,Emulsifier,EM_Concentration,Timepoint)
  Summary.data <- summarise(Grouped.data,mean.ratio= mean(Ratio), sd.ratio = sd(Ratio))

#TO
  S.data.CMC.T0 <- filter(Summary.data,Timepoint =="0", Emulsifier == 'CMC')
  #Mean.R.CMC= mean(S.data.CMC.T0[["mean.ratio.CMC"]])

  Grouped.CMC.T0 <- filter(Grouped.data, Timepoint == "0", Emulsifier == 'CMC')
  Var.R.CMC <- var(Grouped.CMC.T0[["Ratio"]])              
  Sd.R.CMC <- sd(Grouped.CMC.T0[["Ratio"]])                 
  Mean.R.CMC <- mean(Grouped.CMC.T0[["Ratio"]])             

  p <- 0.25
  Eff.Size.CMC.T0 <- as.numeric(sqrt(p*(((S.data.CMC.T0[1,4]-Mean.R.CMC)^2)+((S.data.CMC.T0[2,4]-Mean.R.CMC)^2)+((S.data.CMC.T0[3,4]-Mean.R.CMC)^2)+((S.data.CMC.T0[4,4]-Mean.R.CMC)^2))/Var.R.CMC))
  Sample_Size_CMC_T0 <- pwr.anova.test(k = 4, n = , f = Eff.Size.CMC.T0, sig.level = 0.05, power = 0.8)

  
  #Cohen's D, calculated versus control
  D_0.005_CMC_T0 <- as.numeric((S.data.CMC.T0[2,4]-S.data.CMC.T0[1,4])/Var.R.CMC)
  D_0.05_CMC_T0 <- as.numeric((S.data.CMC.T0[3,4]-S.data.CMC.T0[1,4])/Var.R.CMC)
  D_0.5_CMC_T0 <- as.numeric((S.data.CMC.T0[4,4]-S.data.CMC.T0[1,4])/Var.R.CMC)

  SS_0.005_CMC_T0 <- pwr.t.test(n = , d = D_0.005_CMC_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_CMC_T0 <- pwr.t.test(n = , d = D_0.05_CMC_T0, sig.level = 0.05, power = 0.8, type = "two.sample")    ## => doesn't work 
  SS_0.5_CMC_T0 <- pwr.t.test(n = , d = D_0.5_CMC_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  
#T1
  S.data.CMC.T1 <- filter(Summary.data,Timepoint =="24", Emulsifier == 'CMC')
  #Mean.R.CMC= mean(S.data.CMC.T1[["mean.ratio.CMC"]])
  
  Grouped.CMC.T1 <- filter(Grouped.data, Timepoint == "24", Emulsifier == 'CMC')
  Var.R.CMC <- var(Grouped.CMC.T1[["Ratio"]])
  Sd.R.CMC <- sd(Grouped.CMC.T1[["Ratio"]])
  Mean.R.CMC <- mean(Grouped.CMC.T1[["Ratio"]])
  
  p <- 0.25
  Eff.Size.CMC.T1 <- as.numeric (sqrt(p*(((S.data.CMC.T1[1,4]-Mean.R.CMC)^2)+((S.data.CMC.T1[2,4]-Mean.R.CMC)^2)+((S.data.CMC.T1[3,4]-Mean.R.CMC)^2)+((S.data.CMC.T1[4,4]-Mean.R.CMC)^2))/Var.R.CMC))
  Sample_Size_CMC_T1 <- pwr.anova.test(k = 4, n = , f = Eff.Size.CMC.T1, sig.level = 0.05, power = 0.8)
  
  #Cohen's D, calculated versus control
  D_0.005_CMC_T1 <- as.numeric((S.data.CMC.T1[2,4]-S.data.CMC.T1[1,4])/Var.R.CMC)
  D_0.05_CMC_T1 <- as.numeric((S.data.CMC.T1[3,4]-S.data.CMC.T1[1,4])/Var.R.CMC)
  D_0.5_CMC_T1 <- as.numeric((S.data.CMC.T1[4,4]-S.data.CMC.T1[1,4])/Var.R.CMC)
  
  SS_0.005 <- pwr.t.test(n = , d = D_0.005_CMC_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05 <- pwr.t.test(n = , d = D_0.05_CMC_T1, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5 <- pwr.t.test(n = , d = D_0.5_CMC_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  
#T2
  S.data.CMC.T2 <- filter(Summary.data,Timepoint =="48", Emulsifier == 'CMC')
  #Mean.R.CMC= mean(S.data.CMC.T2[["mean.ratio.CMC"]])
  
  Grouped.CMC.T2 <- filter(Grouped.data, Timepoint == "48", Emulsifier == 'CMC')
  Var.R.CMC <- var(Grouped.CMC.T2[["Ratio"]])
  Sd.R.CMC <- sd(Grouped.CMC.T2[["Ratio"]])
  Mean.R.CMC <- mean(Grouped.CMC.T2[["Ratio"]])
  
  p <- 0.25
  Eff.Size.CMC.T2 <- as.numeric(sqrt(p*(((S.data.CMC.T2[1,4]-Mean.R.CMC)^2)+((S.data.CMC.T2[2,4]-Mean.R.CMC)^2)+((S.data.CMC.T2[3,4]-Mean.R.CMC)^2)+((S.data.CMC.T2[4,4]-Mean.R.CMC)^2))/Var.R.CMC))
  Sample_Size_CMC_T2 <- pwr.anova.test(k = 4, n = , f = Eff.Size.CMC.T2, sig.level = 0.05, power = 0.8)
  
  
  #Cohen's D, calculated versus control
  D_0.005_CMC_T2 <- as.numeric((S.data.CMC.T2[2,4]-S.data.CMC.T2[1,4])/Var.R.CMC)
  D_0.05_CMC_T2 <- as.numeric((S.data.CMC.T2[3,4]-S.data.CMC.T2[1,4])/Var.R.CMC)
  D_0.5_CMC_T2 <- as.numeric((S.data.CMC.T2[4,4]-S.data.CMC.T2[1,4])/Var.R.CMC)
  
  SS_0.005 <- pwr.t.test(n = , d = D_0.005_CMC_T2, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05 <- pwr.t.test(n = , d = D_0.05_CMC_T2, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5 <- pwr.t.test(n = , d = D_0.5_CMC_T2, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  
  
### Tween80
  
  
#TO
  S.data.P80.T0 <- filter(Summary.data,Timepoint =="0", Emulsifier == 'P80')
  #Mean.R.P80= mean(S.data.P80.T0[["mean.ratio.P80"]])
  
  Grouped.P80.T0 <- filter(Grouped.data, Timepoint == "0", Emulsifier == 'P80')
  Var.R.P80 <- var(Grouped.P80.T0[["Ratio"]])
  Sd.R.P80 <- sd(Grouped.P80.T0[["Ratio"]])
  Mean.R.P80 <- mean(Grouped.P80.T0[["Ratio"]])
  
  p <- 0.25
  Eff.Size.P80.T0 <- as.numeric (sqrt(p*(((S.data.P80.T0[1,4]-Mean.R.P80)^2)+((S.data.P80.T0[2,4]-Mean.R.P80)^2)+((S.data.P80.T0[3,4]-Mean.R.P80)^2)+((S.data.P80.T0[4,4]-Mean.R.P80)^2))/Var.R.P80))
  Sample_Size_P80_T0 <- pwr.anova.test(k = 4, n = , f = Eff.Size.P80.T0, sig.level = 0.05, power = 0.8)
  
  #Cohen's D, calculated versus control
  D_0.005_P80_T0 <- as.numeric(S.data.P80.T0[2,4]-S.data.P80.T0[1,4])/Var.R.P80
  D_0.05_P80_T0 <- as.numeric(S.data.P80.T0[3,4]-S.data.P80.T0[1,4])/Var.R.P80
  D_0.5_P80_T0 <- as.numeric(S.data.P80.T0[4,4]-S.data.P80.T0[1,4])/Var.R.P80
  
  SS_0.005 <- pwr.t.test(n = , d = D_0.005, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05 <- pwr.t.test(n = , d = D_0.05, sig.level = 0.05, power = 0.8, type = "two.sample")    
  ## => doesn't work but since the ES is larger here than for the other concentrations, the sample size will be efen smaller
  SS_0.5 <- pwr.t.test(n = , d = D_0.5, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  
#T1
  S.data.P80.T0 <- filter(Summary.data,Timepoint =="24", Emulsifier == 'P80')
  #Mean.R.P80= mean(S.data.P80.T0[["mean.ratio.P80"]])
  
  Grouped.P80.T0 <- filter(Grouped.data, Timepoint == "24", Emulsifier == 'P80')
  Var.R.P80 <- var(Grouped.P80.T0[["Ratio"]])
  Sd.R.P80 <- sd(Grouped.P80.T0[["Ratio"]])
  Mean.R.P80 <- mean(Grouped.P80.T0[["Ratio"]])
  
  p <- 0.25
  Eff.Size.P80.T0 <- as.numeric (sqrt(p*(((S.data.P80.T0[1,4]-Mean.R.P80)^2)+((S.data.P80.T0[2,4]-Mean.R.P80)^2)+((S.data.P80.T0[3,4]-Mean.R.P80)^2)+((S.data.P80.T0[4,4]-Mean.R.P80)^2))/Var.R.P80))
  Sample_Size_P80_T0 <- pwr.anova.test(k = 4, n = , f = Eff.Size.P80.T0, sig.level = 0.05, power = 0.8)
  
  #Cohen's D, calculated versus control
  D_0.005_P80_T0 <- as.numeric(S.data.P80.T0[2,4]-S.data.P80.T0[1,4])/Var.R.P80
  D_0.05_P80_T0 <- as.numeric(S.data.P80.T0[3,4]-S.data.P80.T0[1,4])/Var.R.P80
  D_0.5_P80_T0 <- as.numeric(S.data.P80.T0[4,4]-S.data.P80.T0[1,4])/Var.R.P80
  
  SS_0.005 <- pwr.t.test(n = , d = D_0.005, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05 <- pwr.t.test(n = , d = D_0.05, sig.level = 0.05, power = 0.8, type = "two.sample")    
  ## => doesn't work but since the ES is larger here than for the other concentrations, the sample size will be efen smaller
  SS_0.5 <- pwr.t.test(n = , d = D_0.5, sig.level = 0.05, power = 0.8, type = "two.sample")
  
#T2
  
# Soy Lecithin
  
# T0
  
# T1
  
# T2  
  
# Sophorolipids
  
# T0
  
# T1
  
# T2 
  
# Rhamnolipids
  
# T0
  
# T1
  
# T2 
  
  
  
  
  
  rm(list = ls())