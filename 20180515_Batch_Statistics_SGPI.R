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


levene.test(data.cells, group, location=c("median", "mean", "trim.mean"), trim.alpha=0.25,
            bootstrap = FALSE, num.bootstrap=1000, kruskal.test=FALSE, 
            correction.method=c("none","correction.factor","zero.removal","zero.correction"))



#-----------------------------------------------------------------------------------------------------------------#

### Effect size calculation ###
### Per ANOVA you want to do, the effect size needs to be calculated in order 
### to calc the required sample size.

#Mean.Ratio = mean(data.cells[["Ratio"]])
#Var.Ratio = var(data.cells[["Ratio"]])


### Calculate means and Stdev of Ratio => pull together all donors
  Grouped.data  <- group_by(data.cells,Emulsifier,EM_Concentration,Timepoint)
  Summary.data <- summarise(Grouped.data,mean.ratio= mean(Ratio), sd.ratio = sd(Ratio))

  
  ### create matrix in which sample sizes will be written out. 
  SS_Anova = matrix(nrow=15, ncol=4)
  colnames(SS_Anova) <- c("ANOVA", "T_0.005", "T_0.05", "T_0.5")
  row.names(SS_Anova) <- c("CMC_T0","CMC_T1","CMC_T2","P80_T0","P80_T1","P80_T2","Soy_T0","Soy_T1","Soy_T2","SL_T0","SL_T1","SL_T2","RL_T0","RL_T1","RL_T2")
  
### CMC
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
                          ### k = aantal klassen
                          ### f = calculated effect size
                          ### sig.levele = significance level you'd use in an ANOVA-test
                          ### power = DESIRED power of your ANOVA-test.

  SS_Anova[1,1] <-  as.numeric(Sample_Size_CMC_T0 [2])     # wegschrijven naar matrix 
  
  #Cohen's D, calculated versus control
  D_0.005_CMC_T0 <- abs(as.numeric((S.data.CMC.T0[2,4]-S.data.CMC.T0[1,4])/Var.R.CMC))
  D_0.05_CMC_T0 <- abs(as.numeric((S.data.CMC.T0[3,4]-S.data.CMC.T0[1,4])/Var.R.CMC))
  D_0.5_CMC_T0 <- abs(as.numeric((S.data.CMC.T0[4,4]-S.data.CMC.T0[1,4])/Var.R.CMC))

  SS_0.005_CMC_T0 <- pwr.t.test(n = , d = D_0.005_CMC_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_CMC_T0 <- pwr.t.test(n = , d = D_0.05_CMC_T0, sig.level = 0.05, power = 0.8, type = "two.sample")    ## => doesn't work 
  SS_0.5_CMC_T0 <- pwr.t.test(n = , d = D_0.5_CMC_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[1,2] <-  as.numeric(SS_0.005_CMC_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[1,3] <-  as.numeric(SS_0.05_CMC_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[1,4] <-  as.numeric(SS_0.5_CMC_T0 [1])     # wegschrijven naar matrix 
  
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
  
  SS_Anova[2,1] <-  as.numeric(Sample_Size_CMC_T1 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_CMC_T1 <- abs(as.numeric((S.data.CMC.T1[2,4]-S.data.CMC.T1[1,4])/Var.R.CMC))
  D_0.05_CMC_T1 <- abs(as.numeric((S.data.CMC.T1[3,4]-S.data.CMC.T1[1,4])/Var.R.CMC))
  D_0.5_CMC_T1 <- abs(as.numeric((S.data.CMC.T1[4,4]-S.data.CMC.T1[1,4])/Var.R.CMC))
  
  SS_0.005_CMC_T1 <- pwr.t.test(n = , d = D_0.005_CMC_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_CMC_T1 <- pwr.t.test(n = , d = D_0.05_CMC_T1, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_CMC_T1 <- pwr.t.test(n = , d = D_0.5_CMC_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[2,2] <-  as.numeric(SS_0.005_CMC_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[2,3] <-  as.numeric(SS_0.05_CMC_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[2,4] <-  as.numeric(SS_0.5_CMC_T1 [1])   
  
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
  
  SS_Anova[3,1] <-  as.numeric(Sample_Size_CMC_T2 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_CMC_T2 <- abs(as.numeric((S.data.CMC.T2[2,4]-S.data.CMC.T2[1,4])/Var.R.CMC))
  D_0.05_CMC_T2 <- abs(as.numeric((S.data.CMC.T2[3,4]-S.data.CMC.T2[1,4])/Var.R.CMC))
  D_0.5_CMC_T2 <- abs(as.numeric((S.data.CMC.T2[4,4]-S.data.CMC.T2[1,4])/Var.R.CMC))
  
  SS_0.005_CMC_T2 <- pwr.t.test(n = , d = D_0.005_CMC_T2, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_CMC_T2 <- pwr.t.test(n = , d = D_0.05_CMC_T2, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_CMC_T2 <- pwr.t.test(n = , d = D_0.5_CMC_T2, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[3,2] <-  as.numeric(SS_0.005_CMC_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[3,3] <-  as.numeric(SS_0.05_CMC_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[3,4] <-  as.numeric(SS_0.5_CMC_T2 [1])     # wegschrijven naar matrix 
  
### TWEEN 80
  #TO
  S.data.P80.T0 <- filter(Summary.data,Timepoint =="0", Emulsifier == 'P80')
  #Mean.R.P80= mean(S.data.P80.T0[["mean.ratio.P80"]])
  
  Grouped.P80.T0 <- filter(Grouped.data, Timepoint == "0", Emulsifier == 'P80')
  Var.R.P80 <- var(Grouped.P80.T0[["Ratio"]])              
  Sd.R.P80 <- sd(Grouped.P80.T0[["Ratio"]])                 
  Mean.R.P80 <- mean(Grouped.P80.T0[["Ratio"]])             
  
  p <- 0.25
  Eff.Size.P80.T0 <- as.numeric(sqrt(p*(((S.data.P80.T0[1,4]-Mean.R.P80)^2)+((S.data.P80.T0[2,4]-Mean.R.P80)^2)+((S.data.P80.T0[3,4]-Mean.R.P80)^2)+((S.data.P80.T0[4,4]-Mean.R.P80)^2))/Var.R.P80))
  Sample_Size_P80_T0 <- pwr.anova.test(k = 4, n = , f = Eff.Size.P80.T0, sig.level = 0.05, power = 0.8)
  ### k = aantal klassen
  ### f = calculated effect size
  ### sig.levele = significance level you'd use in an ANOVA-test
  ### power = DESIRED power of your ANOVA-test.
  
  SS_Anova[4,1] <-  as.numeric(Sample_Size_P80_T0 [2])     # wegschrijven naar matrix 
  
  #Cohen's D, calculated versus control
  D_0.005_P80_T0 <- abs(as.numeric((S.data.P80.T0[2,4]-S.data.P80.T0[1,4])/Var.R.P80))
  D_0.05_P80_T0 <- abs(as.numeric((S.data.P80.T0[3,4]-S.data.P80.T0[1,4])/Var.R.P80))
  D_0.5_P80_T0 <- abs(as.numeric((S.data.P80.T0[4,4]-S.data.P80.T0[1,4])/Var.R.P80))
  
  SS_0.005_P80_T0 <- pwr.t.test(n = , d = D_0.005_P80_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_P80_T0 <- pwr.t.test(n = , d = D_0.05_P80_T0, sig.level = 0.05, power = 0.8, type = "two.sample")    ## => doesn't work 
  SS_0.5_P80_T0 <- pwr.t.test(n = , d = D_0.5_P80_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[4,2] <-  as.numeric(SS_0.005_P80_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[4,3] <-  as.numeric(SS_0.05_P80_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[4,4] <-  as.numeric(SS_0.5_P80_T0 [1])     # wegschrijven naar matrix 
  
  
  #T1
  S.data.P80.T1 <- filter(Summary.data,Timepoint =="24", Emulsifier == 'P80')
  #Mean.R.P80= mean(S.data.P80.T1[["mean.ratio.P80"]])
  
  Grouped.P80.T1 <- filter(Grouped.data, Timepoint == "24", Emulsifier == 'P80')
  Var.R.P80 <- var(Grouped.P80.T1[["Ratio"]])
  Sd.R.P80 <- sd(Grouped.P80.T1[["Ratio"]])
  Mean.R.P80 <- mean(Grouped.P80.T1[["Ratio"]])
  
  p <- 0.25
  Eff.Size.P80.T1 <- as.numeric (sqrt(p*(((S.data.P80.T1[1,4]-Mean.R.P80)^2)+((S.data.P80.T1[2,4]-Mean.R.P80)^2)+((S.data.P80.T1[3,4]-Mean.R.P80)^2)+((S.data.P80.T1[4,4]-Mean.R.P80)^2))/Var.R.P80))
  Sample_Size_P80_T1 <- pwr.anova.test(k = 4, n = , f = Eff.Size.P80.T1, sig.level = 0.05, power = 0.8)
  
  SS_Anova[5,1] <-  as.numeric(Sample_Size_P80_T1 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_P80_T1 <- abs(as.numeric((S.data.P80.T1[2,4]-S.data.P80.T1[1,4])/Var.R.P80))
  D_0.05_P80_T1 <- abs(as.numeric((S.data.P80.T1[3,4]-S.data.P80.T1[1,4])/Var.R.P80))
  D_0.5_P80_T1 <- abs(as.numeric((S.data.P80.T1[4,4]-S.data.P80.T1[1,4])/Var.R.P80))
  
  SS_0.005_P80_T1 <- pwr.t.test(n = , d = D_0.005_P80_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_P80_T1 <- pwr.t.test(n = , d = D_0.05_P80_T1, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_P80_T1 <- pwr.t.test(n = , d = D_0.5_P80_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[5,2] <-  as.numeric(SS_0.005_P80_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[5,3] <-  as.numeric(SS_0.05_P80_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[5,4] <-  as.numeric(SS_0.5_P80_T1 [1])     # wegschrijven naar matrix 
  
  #T2
  S.data.P80.T2 <- filter(Summary.data,Timepoint =="48", Emulsifier == 'P80')
  #Mean.R.P80= mean(S.data.P80.T2[["mean.ratio.P80"]])
  
  Grouped.P80.T2 <- filter(Grouped.data, Timepoint == "48", Emulsifier == 'P80')
  Var.R.P80 <- var(Grouped.P80.T2[["Ratio"]])
  Sd.R.P80 <- sd(Grouped.P80.T2[["Ratio"]])
  Mean.R.P80 <- mean(Grouped.P80.T2[["Ratio"]])
  
  p <- 0.25
  Eff.Size.P80.T2 <- as.numeric(sqrt(p*(((S.data.P80.T2[1,4]-Mean.R.P80)^2)+((S.data.P80.T2[2,4]-Mean.R.P80)^2)+((S.data.P80.T2[3,4]-Mean.R.P80)^2)+((S.data.P80.T2[4,4]-Mean.R.P80)^2))/Var.R.P80))
  Sample_Size_P80_T2 <- pwr.anova.test(k = 4, n = , f = Eff.Size.P80.T2, sig.level = 0.05, power = 0.8)
  
  SS_Anova[6,1] <-  as.numeric(Sample_Size_P80_T2 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_P80_T2 <- abs(as.numeric((S.data.P80.T2[2,4]-S.data.P80.T2[1,4])/Var.R.P80))
  D_0.05_P80_T2 <- abs(as.numeric((S.data.P80.T2[3,4]-S.data.P80.T2[1,4])/Var.R.P80))
  D_0.5_P80_T2 <- abs(as.numeric((S.data.P80.T2[4,4]-S.data.P80.T2[1,4])/Var.R.P80))
  
  SS_0.005_P80_T2 <- pwr.t.test(n = , d = D_0.005_P80_T2, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_P80_T2 <- pwr.t.test(n = , d = D_0.05_P80_T2, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_P80_T2 <- pwr.t.test(n = , d = D_0.5_P80_T2, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[6,2] <-  as.numeric(SS_0.005_P80_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[6,3] <-  as.numeric(SS_0.05_P80_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[6,4] <-  as.numeric(SS_0.5_P80_T2 [1])     # wegschrijven naar matrix 
  
### Soy Lecithin
 
  #TO
  S.data.Soy.T0 <- filter(Summary.data,Timepoint =="0", Emulsifier == 'SoyLec')
  #Mean.R.Soy= mean(S.data.Soy.T0[["mean.ratio.Soy"]])
  
  Grouped.Soy.T0 <- filter(Grouped.data, Timepoint == "0", Emulsifier == 'SoyLec')
  Var.R.Soy <- var(Grouped.Soy.T0[["Ratio"]])              
  Sd.R.Soy <- sd(Grouped.Soy.T0[["Ratio"]])                 
  Mean.R.Soy <- mean(Grouped.Soy.T0[["Ratio"]])             
  
  p <- 0.25
  Eff.Size.Soy.T0 <- as.numeric(sqrt(p*(((S.data.Soy.T0[1,4]-Mean.R.Soy)^2)+((S.data.Soy.T0[2,4]-Mean.R.Soy)^2)+((S.data.Soy.T0[3,4]-Mean.R.Soy)^2)+((S.data.Soy.T0[4,4]-Mean.R.Soy)^2))/Var.R.Soy))
  Sample_Size_Soy_T0 <- pwr.anova.test(k = 4, n = , f = Eff.Size.Soy.T0, sig.level = 0.05, power = 0.8)
  ### k = aantal klassen
  ### f = calculated effect size
  ### sig.levele = significance level you'd use in an ANOVA-test
  ### power = DESIRED power of your ANOVA-test.
  
  SS_Anova[7,1] <-  as.numeric(Sample_Size_Soy_T0 [2])     # wegschrijven naar matrix 
  
  #Cohen's D, calculated versus control
  D_0.005_Soy_T0 <- abs(as.numeric((S.data.Soy.T0[2,4]-S.data.Soy.T0[1,4])/Var.R.Soy))
  D_0.05_Soy_T0 <- abs(as.numeric((S.data.Soy.T0[3,4]-S.data.Soy.T0[1,4])/Var.R.Soy))
  D_0.5_Soy_T0 <- abs(as.numeric((S.data.Soy.T0[4,4]-S.data.Soy.T0[1,4])/Var.R.Soy))
  
  SS_0.005_Soy_T0 <- pwr.t.test(n = , d = D_0.005_Soy_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_Soy_T0 <- pwr.t.test(n = , d = D_0.05_Soy_T0, sig.level = 0.05, power = 0.8, type = "two.sample")    ## => doesn't work 
  SS_0.5_Soy_T0 <- pwr.t.test(n = , d = D_0.5_Soy_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[7,2] <-  as.numeric(SS_0.005_Soy_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[7,3] <-  as.numeric(SS_0.05_Soy_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[7,4] <-  as.numeric(SS_0.5_Soy_T0 [1])     # wegschrijven naar matrix 
  
  #T1
  S.data.Soy.T1 <- filter(Summary.data,Timepoint =="24", Emulsifier == 'SoyLec')
  #Mean.R.Soy= mean(S.data.Soy.T1[["mean.ratio.Soy"]])
  
  Grouped.Soy.T1 <- filter(Grouped.data, Timepoint == "24", Emulsifier == 'SoyLec')
  Var.R.Soy <- var(Grouped.Soy.T1[["Ratio"]])
  Sd.R.Soy <- sd(Grouped.Soy.T1[["Ratio"]])
  Mean.R.Soy <- mean(Grouped.Soy.T1[["Ratio"]])
  
  p <- 0.25
  Eff.Size.Soy.T1 <- as.numeric (sqrt(p*(((S.data.Soy.T1[1,4]-Mean.R.Soy)^2)+((S.data.Soy.T1[2,4]-Mean.R.Soy)^2)+((S.data.Soy.T1[3,4]-Mean.R.Soy)^2)+((S.data.Soy.T1[4,4]-Mean.R.Soy)^2))/Var.R.Soy))
  Sample_Size_Soy_T1 <- pwr.anova.test(k = 4, n = , f = Eff.Size.Soy.T1, sig.level = 0.05, power = 0.8)
  
  SS_Anova[8,1] <-  as.numeric(Sample_Size_Soy_T1 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_Soy_T1 <- abs(as.numeric((S.data.Soy.T1[2,4]-S.data.Soy.T1[1,4])/Var.R.Soy))
  D_0.05_Soy_T1 <- abs(as.numeric((S.data.Soy.T1[3,4]-S.data.Soy.T1[1,4])/Var.R.Soy))
  D_0.5_Soy_T1 <- abs(as.numeric((S.data.Soy.T1[4,4]-S.data.Soy.T1[1,4])/Var.R.Soy))
  
  SS_0.005_Soy_T1 <- pwr.t.test(n = , d = D_0.005_Soy_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_Soy_T1 <- pwr.t.test(n = , d = D_0.05_Soy_T1, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_Soy_T1 <- pwr.t.test(n = , d = D_0.5_Soy_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[8,2] <-  as.numeric(SS_0.005_Soy_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[8,3] <-  as.numeric(SS_0.05_Soy_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[8,4] <-  as.numeric(SS_0.5_Soy_T1 [1])     # wegschrijven naar matrix 
  
  #T2
  S.data.Soy.T2 <- filter(Summary.data,Timepoint =="48", Emulsifier == 'SoyLec')
  #Mean.R.Soy= mean(S.data.Soy.T2[["mean.ratio.Soy"]])
  
  Grouped.Soy.T2 <- filter(Grouped.data, Timepoint == "48", Emulsifier == 'SoyLec')
  Var.R.Soy <- var(Grouped.Soy.T2[["Ratio"]])
  Sd.R.Soy <- sd(Grouped.Soy.T2[["Ratio"]])
  Mean.R.Soy <- mean(Grouped.Soy.T2[["Ratio"]])
  
  p <- 0.25
  Eff.Size.Soy.T2 <- as.numeric(sqrt(p*(((S.data.Soy.T2[1,4]-Mean.R.Soy)^2)+((S.data.Soy.T2[2,4]-Mean.R.Soy)^2)+((S.data.Soy.T2[3,4]-Mean.R.Soy)^2)+((S.data.Soy.T2[4,4]-Mean.R.Soy)^2))/Var.R.Soy))
  Sample_Size_Soy_T2 <- pwr.anova.test(k = 4, n = , f = Eff.Size.Soy.T2, sig.level = 0.05, power = 0.8)
  
  SS_Anova[9,1] <-  as.numeric(Sample_Size_Soy_T2 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_Soy_T2 <- abs(as.numeric((S.data.Soy.T2[2,4]-S.data.Soy.T2[1,4])/Var.R.Soy))
  D_0.05_Soy_T2 <- abs(as.numeric((S.data.Soy.T2[3,4]-S.data.Soy.T2[1,4])/Var.R.Soy))
  D_0.5_Soy_T2 <- abs(as.numeric((S.data.Soy.T2[4,4]-S.data.Soy.T2[1,4])/Var.R.Soy))
  
  SS_0.005_Soy_T2 <- pwr.t.test(n = , d = D_0.005_Soy_T2, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_Soy_T2 <- pwr.t.test(n = , d = D_0.05_Soy_T2, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_Soy_T2 <- pwr.t.test(n = , d = D_0.5_Soy_T2, sig.level = 0.05, power = 0.8, type = "two.sample")  
  
  SS_Anova[9,2] <-  as.numeric(SS_0.005_Soy_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[9,3] <-  as.numeric(SS_0.05_Soy_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[9,4] <-  as.numeric(SS_0.5_Soy_T2 [1])     # wegschrijven naar matrix
  
# Sophorolipids
  
  #TO
  S.data.SL.T0 <- filter(Summary.data,Timepoint =="0", Emulsifier == 'SL')
  #Mean.R.SL= mean(S.data.SL.T0[["mean.ratio.SL"]])
  
  Grouped.SL.T0 <- filter(Grouped.data, Timepoint == "0", Emulsifier == 'SL')
  Var.R.SL <- var(Grouped.SL.T0[["Ratio"]])              
  Sd.R.SL <- sd(Grouped.SL.T0[["Ratio"]])                 
  Mean.R.SL <- mean(Grouped.SL.T0[["Ratio"]])             
  
  p <- 0.25
  Eff.Size.SL.T0 <- as.numeric(sqrt(p*(((S.data.SL.T0[1,4]-Mean.R.SL)^2)+((S.data.SL.T0[2,4]-Mean.R.SL)^2)+((S.data.SL.T0[3,4]-Mean.R.SL)^2)+((S.data.SL.T0[4,4]-Mean.R.SL)^2))/Var.R.SL))
  Sample_Size_SL_T0 <- pwr.anova.test(k = 4, n = , f = Eff.Size.SL.T0, sig.level = 0.05, power = 0.8)
  ### k = aantal klassen
  ### f = calculated effect size
  ### sig.levele = significance level you'd use in an ANOVA-test
  ### power = DESIRED power of your ANOVA-test.
  
  SS_Anova[10,1] <-  as.numeric(Sample_Size_SL_T0 [2])     # wegschrijven naar matrix 
  
  #Cohen's D, calculated versus control
  D_0.005_SL_T0 <- abs(as.numeric((S.data.SL.T0[2,4]-S.data.SL.T0[1,4])/Var.R.SL))
  D_0.05_SL_T0 <- abs(as.numeric((S.data.SL.T0[3,4]-S.data.SL.T0[1,4])/Var.R.SL))
  D_0.5_SL_T0 <- abs(as.numeric((S.data.SL.T0[4,4]-S.data.SL.T0[1,4])/Var.R.SL))
  
  SS_0.005_SL_T0 <- pwr.t.test(n = , d = D_0.005_SL_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_SL_T0 <- pwr.t.test(n = , d = D_0.05_SL_T0, sig.level = 0.05, power = 0.8, type = "two.sample")    ## => doesn't work 
  SS_0.5_SL_T0 <- pwr.t.test(n = , d = D_0.5_SL_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  
  SS_Anova[10,2] <-  as.numeric(SS_0.005_SL_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[10,3] <-  as.numeric(SS_0.05_SL_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[10,4] <-  as.numeric(SS_0.5_SL_T0 [1])     # wegschrijven naar matrix
  
  #T1
  S.data.SL.T1 <- filter(Summary.data,Timepoint =="24", Emulsifier == 'SL')
  #Mean.R.SL= mean(S.data.SL.T1[["mean.ratio.SL"]])
  
  Grouped.SL.T1 <- filter(Grouped.data, Timepoint == "24", Emulsifier == 'SL')
  Var.R.SL <- var(Grouped.SL.T1[["Ratio"]])
  Sd.R.SL <- sd(Grouped.SL.T1[["Ratio"]])
  Mean.R.SL <- mean(Grouped.SL.T1[["Ratio"]])
  
  p <- 0.25
  Eff.Size.SL.T1 <- as.numeric (sqrt(p*(((S.data.SL.T1[1,4]-Mean.R.SL)^2)+((S.data.SL.T1[2,4]-Mean.R.SL)^2)+((S.data.SL.T1[3,4]-Mean.R.SL)^2)+((S.data.SL.T1[4,4]-Mean.R.SL)^2))/Var.R.SL))
  Sample_Size_SL_T1 <- pwr.anova.test(k = 4, n = , f = Eff.Size.SL.T1, sig.level = 0.05, power = 0.8)
  
  SS_Anova[11,1] <-  as.numeric(Sample_Size_SL_T1 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_SL_T1 <- abs(as.numeric((S.data.SL.T1[2,4]-S.data.SL.T1[1,4])/Var.R.SL))
  D_0.05_SL_T1 <- abs(as.numeric((S.data.SL.T1[3,4]-S.data.SL.T1[1,4])/Var.R.SL))
  D_0.5_SL_T1 <- abs(as.numeric((S.data.SL.T1[4,4]-S.data.SL.T1[1,4])/Var.R.SL))
  
  SS_0.005_SL_T1 <- pwr.t.test(n = , d = D_0.005_SL_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_SL_T1 <- pwr.t.test(n = , d = D_0.05_SL_T1, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_SL_T1 <- pwr.t.test(n = , d = D_0.5_SL_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[11,2] <-  as.numeric(SS_0.005_SL_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[11,3] <-  as.numeric(SS_0.05_SL_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[11,4] <-  as.numeric(SS_0.5_SL_T1 [1])     # wegschrijven naar matrix
  
  #T2
  S.data.SL.T2 <- filter(Summary.data,Timepoint =="48", Emulsifier == 'SL')
  #Mean.R.SL= mean(S.data.SL.T2[["mean.ratio.SL"]])
  
  Grouped.SL.T2 <- filter(Grouped.data, Timepoint == "48", Emulsifier == 'SL')
  Var.R.SL <- var(Grouped.SL.T2[["Ratio"]])
  Sd.R.SL <- sd(Grouped.SL.T2[["Ratio"]])
  Mean.R.SL <- mean(Grouped.SL.T2[["Ratio"]])
  
  p <- 0.25
  Eff.Size.SL.T2 <- as.numeric(sqrt(p*(((S.data.SL.T2[1,4]-Mean.R.SL)^2)+((S.data.SL.T2[2,4]-Mean.R.SL)^2)+((S.data.SL.T2[3,4]-Mean.R.SL)^2)+((S.data.SL.T2[4,4]-Mean.R.SL)^2))/Var.R.SL))
  Sample_Size_SL_T2 <- pwr.anova.test(k = 4, n = , f = Eff.Size.SL.T2, sig.level = 0.05, power = 0.8)
  
  SS_Anova[12,1] <-  as.numeric(Sample_Size_SL_T2 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_SL_T2 <- abs(as.numeric((S.data.SL.T2[2,4]-S.data.SL.T2[1,4])/Var.R.SL))
  D_0.05_SL_T2 <- abs(as.numeric((S.data.SL.T2[3,4]-S.data.SL.T2[1,4])/Var.R.SL))
  D_0.5_SL_T2 <- abs(as.numeric((S.data.SL.T2[4,4]-S.data.SL.T2[1,4])/Var.R.SL))
  
  SS_0.005_SL_T2 <- pwr.t.test(n = , d = D_0.005_SL_T2, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_SL_T2 <- pwr.t.test(n = , d = D_0.05_SL_T2, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_SL_T2 <- pwr.t.test(n = , d = D_0.5_SL_T2, sig.level = 0.05, power = 0.8, type = "two.sample")  
  
  SS_Anova[12,2] <-  as.numeric(SS_0.005_SL_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[12,3] <-  as.numeric(SS_0.05_SL_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[12,4] <-  as.numeric(SS_0.5_SL_T2 [1])     # wegschrijven naar matrix
  
# Rhamnolipids
  
  #TO
  S.data.RL.T0 <- filter(Summary.data,Timepoint =="0", Emulsifier == 'RL')
  #Mean.R.RL= mean(S.data.RL.T0[["mean.ratio.RL"]])
  
  Grouped.RL.T0 <- filter(Grouped.data, Timepoint == "0", Emulsifier == 'RL')
  Var.R.RL <- var(Grouped.RL.T0[["Ratio"]])              
  Sd.R.RL <- sd(Grouped.RL.T0[["Ratio"]])                 
  Mean.R.RL <- mean(Grouped.RL.T0[["Ratio"]])             
  
  p <- 0.25
  Eff.Size.RL.T0 <- as.numeric(sqrt(p*(((S.data.RL.T0[1,4]-Mean.R.RL)^2)+((S.data.RL.T0[2,4]-Mean.R.RL)^2)+((S.data.RL.T0[3,4]-Mean.R.RL)^2)+((S.data.RL.T0[4,4]-Mean.R.RL)^2))/Var.R.RL))
  Sample_Size_RL_T0 <- pwr.anova.test(k = 4, n = , f = Eff.Size.RL.T0, sig.level = 0.05, power = 0.8)
  ### k = aantal klassen
  ### f = calculated effect size
  ### sig.levele = significance level you'd use in an ANOVA-test
  ### power = DESIRED power of your ANOVA-test.
  
  SS_Anova[13,1] <-  as.numeric(Sample_Size_RL_T0 [2])     # wegschrijven naar matrix 
  
  #Cohen's D, calculated versus control
  D_0.005_RL_T0 <- abs(as.numeric((S.data.RL.T0[2,4]-S.data.RL.T0[1,4])/Var.R.RL))
  D_0.05_RL_T0 <- abs(as.numeric((S.data.RL.T0[3,4]-S.data.RL.T0[1,4])/Var.R.RL))
  D_0.5_RL_T0 <- abs(as.numeric((S.data.RL.T0[4,4]-S.data.RL.T0[1,4])/Var.R.RL))
  
  SS_0.005_RL_T0 <- pwr.t.test(n = , d = D_0.005_RL_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_RL_T0 <- pwr.t.test(n = , d = D_0.05_RL_T0, sig.level = 0.05, power = 0.8, type = "two.sample")    ## => doesn't work 
  SS_0.5_RL_T0 <- pwr.t.test(n = , d = D_0.5_RL_T0, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[13,2] <-  as.numeric(SS_0.005_RL_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[13,3] <-  as.numeric(SS_0.05_RL_T0 [1])     # wegschrijven naar matrix 
  SS_Anova[13,4] <-  as.numeric(SS_0.5_RL_T0 [1])     # wegschrijven naar matrix
  
  #T1
  S.data.RL.T1 <- filter(Summary.data,Timepoint =="24", Emulsifier == 'RL')
  #Mean.R.RL= mean(S.data.RL.T1[["mean.ratio.RL"]])
  
  Grouped.RL.T1 <- filter(Grouped.data, Timepoint == "24", Emulsifier == 'RL')
  Var.R.RL <- var(Grouped.RL.T1[["Ratio"]])
  Sd.R.RL <- sd(Grouped.RL.T1[["Ratio"]])
  Mean.R.RL <- mean(Grouped.RL.T1[["Ratio"]])
  
  p <- 0.25
  Eff.Size.RL.T1 <- as.numeric (sqrt(p*(((S.data.RL.T1[1,4]-Mean.R.RL)^2)+((S.data.RL.T1[2,4]-Mean.R.RL)^2)+((S.data.RL.T1[3,4]-Mean.R.RL)^2)+((S.data.RL.T1[4,4]-Mean.R.RL)^2))/Var.R.RL))
  Sample_Size_RL_T1 <- pwr.anova.test(k = 4, n = , f = Eff.Size.RL.T1, sig.level = 0.05, power = 0.8)
  
  SS_Anova[14,1] <-  as.numeric(Sample_Size_RL_T1 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_RL_T1 <- abs(as.numeric((S.data.RL.T1[2,4]-S.data.RL.T1[1,4])/Var.R.RL))
  D_0.05_RL_T1 <- abs(as.numeric((S.data.RL.T1[3,4]-S.data.RL.T1[1,4])/Var.R.RL))
  D_0.5_RL_T1 <- abs(as.numeric((S.data.RL.T1[4,4]-S.data.RL.T1[1,4])/Var.R.RL))
  
  SS_0.005_RL_T1 <- pwr.t.test(n = , d = D_0.005_RL_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_RL_T1 <- pwr.t.test(n = , d = D_0.05_RL_T1, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_RL_T1 <- pwr.t.test(n = , d = D_0.5_RL_T1, sig.level = 0.05, power = 0.8, type = "two.sample")
  
  SS_Anova[14,2] <-  as.numeric(SS_0.005_RL_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[14,3] <-  as.numeric(SS_0.05_RL_T1 [1])     # wegschrijven naar matrix 
  SS_Anova[14,4] <-  as.numeric(SS_0.5_RL_T1 [1])     # wegschrijven naar matrix
  
  #T2
  S.data.RL.T2 <- filter(Summary.data,Timepoint =="48", Emulsifier == 'RL')
  #Mean.R.RL= mean(S.data.RL.T2[["mean.ratio.RL"]])
  
  Grouped.RL.T2 <- filter(Grouped.data, Timepoint == "48", Emulsifier == 'RL')
  Var.R.RL <- var(Grouped.RL.T2[["Ratio"]])
  Sd.R.RL <- sd(Grouped.RL.T2[["Ratio"]])
  Mean.R.RL <- mean(Grouped.RL.T2[["Ratio"]])
  
  p <- 0.25
  Eff.Size.RL.T2 <- as.numeric(sqrt(p*(((S.data.RL.T2[1,4]-Mean.R.RL)^2)+((S.data.RL.T2[2,4]-Mean.R.RL)^2)+((S.data.RL.T2[3,4]-Mean.R.RL)^2)+((S.data.RL.T2[4,4]-Mean.R.RL)^2))/Var.R.RL))
  Sample_Size_RL_T2 <- pwr.anova.test(k = 4, n = , f = Eff.Size.RL.T2, sig.level = 0.05, power = 0.8)
  
  SS_Anova[15,1] <-  as.numeric(Sample_Size_RL_T2 [2])  # wegschrijven naar matrix
  
  #Cohen's D, calculated versus control
  D_0.005_RL_T2 <- abs(as.numeric((S.data.RL.T2[2,4]-S.data.RL.T2[1,4])/Var.R.RL))
  D_0.05_RL_T2 <- abs(as.numeric((S.data.RL.T2[3,4]-S.data.RL.T2[1,4])/Var.R.RL))
  D_0.5_RL_T2 <- abs(as.numeric((S.data.RL.T2[4,4]-S.data.RL.T2[1,4])/Var.R.RL))
  
  SS_0.005_RL_T2 <- pwr.t.test(n = , d = D_0.005_RL_T2, sig.level = 0.05, power = 0.8, type = "two.sample")
  SS_0.05_RL_T2 <- pwr.t.test(n = , d = D_0.05_RL_T2, sig.level = 0.05, power = 0.8, type = "two.sample")    
  SS_0.5_RL_T2 <- pwr.t.test(n = , d = D_0.5_RL_T2, sig.level = 0.05, power = 0.8, type = "two.sample") 
  
  SS_Anova[15,2] <-  as.numeric(SS_0.005_RL_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[15,3] <-  as.numeric(SS_0.05_RL_T2 [1])     # wegschrijven naar matrix 
  SS_Anova[15,4] <-  as.numeric(SS_0.5_RL_T2 [1])     # wegschrijven naar matrix
  
  library(xlsx)
  write.xlsx(SS_Anova, 'Batch_SampleSizes_SGPI_Ratio.xlsx', sheetName="Sheet1" )
  

  rm(list = ls())