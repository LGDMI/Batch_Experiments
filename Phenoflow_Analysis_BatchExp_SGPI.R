###Phenoflow analysis  Batch-experiments    SGPI-data
setwd("/media/projects2/LisaM/Batch_Experiments/Batch_Experiments")
# install packages & load libraries
#install.packages("ggplot2")
#install.packages('dplyr')
#install.packages('tidyr')

library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

#install_github("rprops/Phenoflow_package")

library("Phenoflow") # for fingerprinting
library("flowViz") # for plotting
library("ggplot2") # for plotting
library("flowAI") # for denoising

set.seed(777)


### DONOR 1 ###

path = "D1_-4_AllTimeppoints"
flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

#$frames
#$phenoData
#$colnames
#$class


#------------------------------------------------------------------------------------------------------------------#

### Denoise data => trek gates


# Select phenotypic features of interest and transform parameters
flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
remove(flowData)



### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)

# Intact Cells
    sqrcut1 <- matrix(c(9.5,9.5,14,14,6,8.2,13.2,10),ncol=2, nrow=4)
    colnames(sqrcut1) <- c("FL1-H","FL3-H")
    polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

# Damaged Cells
    sqrcut2 <- matrix(c(8.25,8.25,12,12,9,11.2,15,13),ncol=2, nrow=4)
    colnames(sqrcut2) <- c("FL1-H","FL3-H")
    polyGate2 <- polygonGate(.gate=sqrcut2, filterId = "Total Cells")

sqrcut3 <- matrix(c(9.5,8.25,16,16,5,12.5,18,10),ncol=2, nrow=4)
colnames(sqrcut3) <- c("FL1-H","FL3-H")
polyGate3 <- polygonGate(.gate=sqrcut3, filterId = "Total Cells")

###  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[31], filter=polyGate3,
       scales=list(y=list(limits=c(0,20)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)


### Isolate only the cellular information based on the polyGate3
flowData_transformed <- Subset(flowData_transformed, polyGate3)

### Extract metadata from sample names
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed),"-"), rbind)))
colnames(metadata) <- c("Emulsifier", "EM_Conc", "Dilution", "Timepoint", "Donor")


summary <- fsApply(x = flowData_transformed, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
maxval <- max(summary[,9])
mytrans <- function(x) x/maxval
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))

#------------------------------------------------------------------------------------------------------------------#


### FINGERPRINTING ###

### Randomly resample to the lowest sample size
flowData_transformed <- FCS_resample(flowData_transformed, replace=TRUE)
### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

#------------------------------------------------------------------------------------------------------------------#


### CALCULATE ALPHA DIVERSITY ###

### Calculate Diversity from normalized fingerprint 
Diversity.fbasis <- Diversity(fbasis,d=3,plot=FALSE, R=999)


p1 <- ggplot(data = Diversity.fbasis, aes(x = as.numeric(as.character(metadata$EM_Conc)), y = D2, color = metadata$Emulsifier))+
  geom_point(size = 8, alpha = 0.7)+
  geom_line()+
  scale_color_manual(values = c("#a65628", "red", 
                                "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
  theme_bw()+
  labs(color = "EMulsifier", y = "Phenotypic diversity (D2)", x = "EM-Conc")+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")
print(p1)

# Diversity assessment with cleaning
Diversity.clean <- Diversity_rf(flowData_transformed, param = param, R = 3, R.b = 3,
                                cleanFCS = TRUE)

p2 <- ggplot(data = Diversity.clean, aes(x = as.numeric(as.character(metadata$Emulsifier)), y = D2, color = metadata$EM_Conc))+
  geom_point(size = 8, alpha = 0.7)+
  geom_line()+
  scale_color_manual(values = c("#a65628", "red", 
                                "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
  theme_bw()+
  labs(color = "Reactor phase", y = "Phenotypic diversity (D2)", x = "Days")+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")
print(p2)

### Export diversity estimates to .csv file in the chosen directory
write.csv2(file="results.metrics.csv", Diversity.clean)

#------------------------------------------------------------------------------------------------------------------#

### Beta Diversity ANalysis ###

# Beta-diversity assessment of fingerprint
beta.div <- beta_div_fcm(fbasis, ord.type="PCoA")

# Plot ordination
plot_beta_fcm(beta.div, shape = metadata$EM_Conc, color = metadata$Emulsifier, labels="Emulsifier") + 
  theme_bw() +
  geom_point(size = 4, alpha = 0.5)

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#

### DONOR 6 ###

path = "D1_-4_AllTimeppoints"
flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

#$frames
#$phenoData
#$colnames
#$class


#------------------------------------------------------------------------------------------------------------------#

### Denoise data => trek gates


# Select phenotypic features of interest and transform parameters
flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
remove(flowData)



### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)

# Intact Cells
sqrcut1 <- matrix(c(9.5,9.5,14,14,6,8.2,13.2,10),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

# Damaged Cells
sqrcut2 <- matrix(c(8.25,8.25,12,12,9,11.2,15,13),ncol=2, nrow=4)
colnames(sqrcut2) <- c("FL1-H","FL3-H")
polyGate2 <- polygonGate(.gate=sqrcut2, filterId = "Total Cells")

sqrcut3 <- matrix(c(9.5,8.25,16,16,5,12.5,18,10),ncol=2, nrow=4)
colnames(sqrcut3) <- c("FL1-H","FL3-H")
polyGate3 <- polygonGate(.gate=sqrcut3, filterId = "Total Cells")

###  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[31], filter=polyGate3,
       scales=list(y=list(limits=c(0,20)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)


### Isolate only the cellular information based on the polyGate3
flowData_transformed <- Subset(flowData_transformed, polyGate3)

### Extract metadata from sample names
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed),"-"), rbind)))
#metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed),"-", rbind))))
colnames(metadata) <- c("Emulsifier", "EM_Conc", "Dilution", "Timepoint", "Donor")


summary <- fsApply(x = flowData_transformed, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
maxval <- max(summary[,9])
mytrans <- function(x) x/maxval
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))

#------------------------------------------------------------------------------------------------------------------#


### FINGERPRINTING ###

### Randomly resample to the lowest sample size
flowData_transformed <- FCS_resample(flowData_transformed, replace=TRUE)
### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

#------------------------------------------------------------------------------------------------------------------#


### CALCULATE ALPHA DIVERSITY ###

### Calculate Diversity from normalized fingerprint 
Diversity.fbasis <- Diversity(fbasis,d=3,plot=FALSE, R=999)


p1 <- ggplot(data = Diversity.fbasis, aes(x = as.numeric(as.character(metadata$EM_Conc)), y = D2, color = metadata$Emulsifier))+
  geom_point(size = 8, alpha = 0.7)+
  geom_line()+
  scale_color_manual(values = c("#a65628", "red", 
                                "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
  theme_bw()+
  labs(color = "EMulsifier", y = "Phenotypic diversity (D2)", x = "EM-Conc")+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")
print(p1)

# Diversity assessment with cleaning
Diversity.clean <- Diversity_rf(flowData_transformed, param = param, R = 3, R.b = 3,
                                cleanFCS = TRUE)

p2 <- ggplot(data = Diversity.clean, aes(x = as.numeric(as.character(metadata$Emulsifier)), y = D2, color = metadata$EM_Conc))+
  geom_point(size = 8, alpha = 0.7)+
  geom_line()+
  scale_color_manual(values = c("#a65628", "red", 
                                "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
  theme_bw()+
  labs(color = "Reactor phase", y = "Phenotypic diversity (D2)", x = "Days")+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")
print(p2)

### Export diversity estimates to .csv file in the chosen directory
write.csv2(file="results.metrics.csv", Diversity.clean)

#------------------------------------------------------------------------------------------------------------------#

### Beta Diversity ANalysis ###

# Beta-diversity assessment of fingerprint
beta.div <- beta_div_fcm(fbasis, ord.type="PCoA")

# Plot ordination
plot_beta_fcm(beta.div, shape = metadata$EM_Conc, color = metadata$Emulsifier, labels="Emulsifier") + 
  theme_bw() +
  geom_point(size = 4, alpha = 0.5)



#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#


### Extract Cellcounts ###

### Creating a rectangle gate for counting HNA and LNA cells
rGate_HNA <- rectangleGate("FL1-H"=c(asinh(20000), 20)/maxval,"FL3-H"=c(0,20)/maxval, 
                           filterId = "HNA bacteria")
### Normalize total cell gate
sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3)/maxval,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

### Check if rectangle gate is correct, if not, adjust rGate_HNA
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=rGate_HNA,
       scales=list(y=list(limits=c(0,1)),
                   x=list(limits=c(0.4,1))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, 
                                                          cex=2), smooth=FALSE)
### Extract the cell counts
a <- flowCore::filter(flowData_transformed, rGate_HNA) 
HNACount <- summary(a);HNACount <- toTable(HNACount)
s <- flowCore::filter(flowData_transformed, polyGate1)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)

### Extract the volume
vol <- as.numeric(flowCore::fsApply(flowData_transformed, FUN = function(x) x@description$`$VOL`))/1000

### Store the data
results_counts <- data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                             Total.cells = TotalCount$true/vol, HNA.cells = HNACount$true/vol)


### Exporting cell counts to .csv file to working directory
write.csv2(file="results.counts.csv", results_counts)



### Plot cell density
ggplot(data = results_counts, aes(x = as.numeric(as.character(metadata$day)), y = Total.cells, color = metadata$Reactor_phase))+
  geom_point(size = 8, alpha = 0.9)+
  scale_color_manual(values = c("#a65628", "red", 
                                "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
  theme_bw()+
  labs(color = "Reactor phase", y = "Total cell density (cells/?L)", x = "Days")  
