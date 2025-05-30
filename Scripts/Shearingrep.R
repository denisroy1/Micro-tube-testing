#### Shearingrep.R ####

# Script written to assesses the repeatability of the 
# shearing process for subset of individual MicroTubes 
# used for the MicroTube paper

# Written by dr May 2025

#### Preliminaries ####

## Clearing all instances
rm(list=ls())

# loading libraries
{
  library(lmerTest)
  library(lme4)
  library(car)
  library(ggplot2)
}

# Set the working directory
setwd("/path/to your/shearing/file/")

# Load the data. Looking for the "shearingrep.csv" file
sheardat<-read.csv(file.choose(), header=T, stringsAsFactors = T)

# Look at the first few lines of data
head(sheardat)

# Set individuals as a factor for easier use
sheardat$ind <- as.factor(sheardat$individual)

#### Repeatability Tests ####

# Use the lmerTest package and lmer to do a repeated measure ANOVA 
# using individual as the random factor.
reptest<-lmer(shearingpeak ~ 1 + (1|ind), data = sheardat)

# Get the summary stats showing a significant difference among 
# and within measurements
summary(reptest)

# Get the confidence intervals for the estimates.
repconf <- confint(reptest)

# Set the VarCorr results as a dataframe to access results
repmat <- as.data.frame(VarCorr(reptest))

# Calculate repeatability as per Nosil & Crespi 2006, Whitlock and Schluter 2020)
repbty <- (repmat$vcov[1])/(repmat$vcov[1] + repmat$vcov[2])
repbty

# Using the confint calculations above to calculate the 95% confidence in the 
# repeatibility estimate
repbl <- (repconf[1]^2)/(repconf[1]^2 + repconf[2]^2)
repbh <- (repconf[4]^2)/(repconf[4]^2 + repconf[5]^2)

# View the low 2.5% CI for repbty
repbl

# View the high 97% CI for repbty
repbh

#### Visualise Repeatability ####

# Plot showing the repeatability of shearing peaks (in bp) of individual 
# DNA extracted and sheared in re-used tubes (i.e., using the same DNA).
ggplot(sheardat, aes(y = shearingpeak, x = ind)) +  
  geom_point(size = 5, col = "firebrick4", alpha = 0.6, stroke = 1) + 
  geom_line(aes(group = ind), linewidth = 0.6, colour = "black") +
  geom_hline(yintercept = 150, linewidth = 0.6, color = "gray", linetype = "dashed") +
  geom_hline(yintercept = 300, linewidth = 0.6, color = "gray", linetype = "dashed") +
  labs(x = "microTube", y = "Shearing peak in BP") + 
  scale_y_continuous(limits = c(140, 300), breaks = seq(140, 320, by = 20)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 22),        
        axis.line = element_line(linewidth = 1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3, "cm"))

#### END ####
