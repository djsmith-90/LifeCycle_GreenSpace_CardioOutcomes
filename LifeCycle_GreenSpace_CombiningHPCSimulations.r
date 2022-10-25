### Combining simulation results for structured lifecourse models with confounders and interactions

### Created 23/8/2022 by Dan Major-Smith
### R version 4.0.4

### Aim: To test the power of this structured life course approach with interactions to select the 'true' interaction term under varying conditions. Will explore all combinations of the various parameters:
#   - Sample sizes: 1,000 vs 10,000
#   - Exposures: Binary (access to green space; yes/no) vs continuous (distance to green space), and centered vs uncentered. Also want to vary the correlation between exposures to see how collinearity impacts the power of the lasso to detect the true model
#   - Outcome: Binary (overweight/obese) vs continuous (BMI)

## For these simulations, will use the same set-up as in the example simulation script, with 'access to green space' as the exposure measured at three time points, cardiometabolic health as the outcome (BMI/obesity), and SEP as a confounder/interaction term. SEP causes access to green space and the outcome (lower BMI/obesity if higher SEP), while the interaction between SEP and the most recent green space time-point also causes the outcome (access to green space causes lower BMI/obesity, but only in those from lower SEP backgrounds).

## In addition to this scenario, we will also vary the strength of the interaction term, to explore how this impacts the power to detect the interaction, as well as varying the specific life course interaction (i.e., the main model will explore an interaction with the most recent critical period, while other simulations will examine interactions with accumulation and change, to see whether this impacts conclusions). Given the number of simulations to run, and that each simulation takes about 6 hours, processing time on a standard laptop may be a bit prohibitive, so will use this script to set-up and test the simulation/code, and then run the actual simulations using University of Bristol's High Performance Computing suite. This script takes the results of these HPC simulation scripts, combines them together, and analyses them


#######################################################################################################
#### Clear workspace and install/load packages

rm(list = ls())

#install.packages("tidyverse")
#install.packages("GGally")

library(tidyverse)
library(GGally)


## Working directory
setwd("C:\\Users\\ds16565\\OneDrive - University of Bristol\\MyFiles-Migrated\\Documents\\Projects\\Lifecycle\\SimResults\\HPC_Results")


### Read in each of the simulation results files and then combine together

## Critical period - No interaction
crit_noInt <- read_csv("simulationResults_crit_noInt.csv")
head(crit_noInt)
summary(crit_noInt)

crit_noInt <- crit_noInt %>%
  mutate(Model = "Crit", Interaction = "None") %>%
  relocate(c(Model, Interaction))

## Critical period - Very small interaction
crit_vSmallInt <- read_csv("simulationResults_crit_vSmallInt.csv")
head(crit_vSmallInt)
summary(crit_vSmallInt)

crit_vSmallInt <- crit_vSmallInt %>%
  mutate(Model = "Crit", Interaction = "vSmall") %>%
  relocate(c(Model, Interaction))

## Critical period - Small interaction
crit_smallInt <- read_csv("simulationResults_crit_smallInt.csv")
head(crit_smallInt)
summary(crit_smallInt)

crit_smallInt <- crit_smallInt %>%
  mutate(Model = "Crit", Interaction = "small") %>%
  relocate(c(Model, Interaction))

## Critical period - Moderate interaction
crit_moderateInt <- read_csv("simulationResults_crit_moderateInt.csv")
head(crit_moderateInt)
summary(crit_moderateInt)

crit_moderateInt <- crit_moderateInt %>%
  mutate(Model = "Crit", Interaction = "moderate") %>%
  relocate(c(Model, Interaction))

## Critical period - Large interaction
crit_largeInt <- read_csv("simulationResults_crit_largeInt.csv")
head(crit_largeInt)
summary(crit_largeInt)

crit_largeInt <- crit_largeInt %>%
  mutate(Model = "Crit", Interaction = "large") %>%
  relocate(c(Model, Interaction))

## Critical period - Very large interaction
crit_vLargeInt <- read_csv("simulationResults_crit_vLargeInt.csv")
head(crit_vLargeInt)
summary(crit_vLargeInt)

crit_vLargeInt <- crit_vLargeInt %>%
  mutate(Model = "Crit", Interaction = "vLarge") %>%
  relocate(c(Model, Interaction))


### Combine these into a single 'critical period' dataset
crit_res <- rbind(crit_noInt, crit_vSmallInt, crit_smallInt, 
                  crit_moderateInt, crit_largeInt, crit_vLargeInt)
head(crit_res)
summary(crit_res)

# Recode the string vars as factors
crit_res <- crit_res %>%
  mutate(Interaction = factor(Interaction, levels = c("None", "vSmall", "small", "moderate", 
                                          "large", "vLarge"))) %>%
  mutate(Exposure = factor(Exposure, levels = c("Binary", "Cont"))) %>%
  mutate(Centered = factor(Centered, levels = c("No", "Yes"))) %>%
  mutate(Collinear = factor(Collinear, levels = c("Low", "High"))) %>%
  mutate(Outcome = factor(Outcome, levels = c("Cont", "Binary")))

head(crit_res)
summary(crit_res)


## Accumulation - No interaction
accum_noInt <- read_csv("simulationResults_accum_noInt.csv")
head(accum_noInt)
summary(accum_noInt)

accum_noInt <- accum_noInt %>%
  mutate(Model = "Accum", Interaction = "None") %>%
  relocate(c(Model, Interaction))

## Accumulation - Very small interaction
accum_vSmallInt <- read_csv("simulationResults_accum_vSmallInt.csv")
head(accum_vSmallInt)
summary(accum_vSmallInt)

accum_vSmallInt <- accum_vSmallInt %>%
  mutate(Model = "Accum", Interaction = "vSmall") %>%
  relocate(c(Model, Interaction))

## Accumulation - Small interaction
accum_smallInt <- read_csv("simulationResults_accum_smallInt.csv")
head(accum_smallInt)
summary(accum_smallInt)

accum_smallInt <- accum_smallInt %>%
  mutate(Model = "Accum", Interaction = "small") %>%
  relocate(c(Model, Interaction))

## Accumulation - Moderate interaction
accum_moderateInt <- read_csv("simulationResults_accum_moderateInt.csv")
head(accum_moderateInt)
summary(accum_moderateInt)

accum_moderateInt <- accum_moderateInt %>%
  mutate(Model = "Accum", Interaction = "moderate") %>%
  relocate(c(Model, Interaction))

## Accumulation - Large interaction
accum_largeInt <- read_csv("simulationResults_accum_largeInt.csv")
head(accum_largeInt)
summary(accum_largeInt)

accum_largeInt <- accum_largeInt %>%
  mutate(Model = "Accum", Interaction = "large") %>%
  relocate(c(Model, Interaction))

## Accumulation - Very large interaction
accum_vLargeInt <- read_csv("simulationResults_accum_vLargeInt.csv")
head(accum_vLargeInt)
summary(accum_vLargeInt)

accum_vLargeInt <- accum_vLargeInt %>%
  mutate(Model = "Accum", Interaction = "vLarge") %>%
  relocate(c(Model, Interaction))


### Combine these into a single 'accumulation' dataset
accum_res <- rbind(accum_noInt, accum_vSmallInt, accum_smallInt, 
                   accum_moderateInt, accum_largeInt, accum_vLargeInt)
head(accum_res)
summary(accum_res)

# Recode the string vars as factors
accum_res <- accum_res %>%
  mutate(Interaction = factor(Interaction, levels = c("None", "vSmall", "small", "moderate", 
                                                      "large", "vLarge"))) %>%
  mutate(Exposure = factor(Exposure, levels = c("Binary", "Cont"))) %>%
  mutate(Centered = factor(Centered, levels = c("No", "Yes"))) %>%
  mutate(Collinear = factor(Collinear, levels = c("Low", "High"))) %>%
  mutate(Outcome = factor(Outcome, levels = c("Cont", "Binary")))

head(accum_res)
summary(accum_res)


## Change from time 2 to time 3 - No interaction
change_noInt <- read_csv("simulationResults_change_noInt.csv")
head(change_noInt)
summary(change_noInt)

change_noInt <- change_noInt %>%
  mutate(Model = "Change", Interaction = "None") %>%
  relocate(c(Model, Interaction))

## Change from time 2 to time 3 - Very small interaction
change_vSmallInt <- read_csv("simulationResults_change_vSmallInt.csv")
head(change_vSmallInt)
summary(change_vSmallInt)

change_vSmallInt <- change_vSmallInt %>%
  mutate(Model = "Change", Interaction = "vSmall") %>%
  relocate(c(Model, Interaction))

## Change from time 2 to time 3 - Small interaction
change_smallInt <- read_csv("simulationResults_change_smallInt.csv")
head(change_smallInt)
summary(change_smallInt)

change_smallInt <- change_smallInt %>%
  mutate(Model = "Change", Interaction = "small") %>%
  relocate(c(Model, Interaction))

## Change from time 2 to time 3 - Moderate interaction
change_moderateInt <- read_csv("simulationResults_change_moderateInt.csv")
head(change_moderateInt)
summary(change_moderateInt)

change_moderateInt <- change_moderateInt %>%
  mutate(Model = "Change", Interaction = "moderate") %>%
  relocate(c(Model, Interaction))

## Change from time 2 to time 3 - Large interaction
change_largeInt <- read_csv("simulationResults_change_largeInt.csv")
head(change_largeInt)
summary(change_largeInt)

change_largeInt <- change_largeInt %>%
  mutate(Model = "Change", Interaction = "large") %>%
  relocate(c(Model, Interaction))

## Change from time 2 to time 3 - Very large interaction
change_vLargeInt <- read_csv("simulationResults_change_vLargeInt.csv")
head(change_vLargeInt)
summary(change_vLargeInt)

change_vLargeInt <- change_vLargeInt %>%
  mutate(Model = "Change", Interaction = "vLarge") %>%
  relocate(c(Model, Interaction))


### Combine these into a single 'change' dataset
change_res <- rbind(change_noInt, change_vSmallInt, change_smallInt, 
                   change_moderateInt, change_largeInt, change_vLargeInt)
head(change_res)
summary(change_res)

# Recode the string vars as factors
change_res <- change_res %>%
  mutate(Interaction = factor(Interaction, levels = c("None", "vSmall", "small", "moderate", 
                                                      "large", "vLarge"))) %>%
  mutate(Exposure = factor(Exposure, levels = c("Binary", "Cont"))) %>%
  mutate(Centered = factor(Centered, levels = c("No", "Yes"))) %>%
  mutate(Collinear = factor(Collinear, levels = c("Low", "High"))) %>%
  mutate(Outcome = factor(Outcome, levels = c("Cont", "Binary")))

head(change_res)
summary(change_res)


### Finally, combine all the simulation results together
all_res <- rbind(crit_res, accum_res, change_res)
head(all_res)
summary(all_res)

# Recode the string vars as factors
all_res <- all_res %>%
  mutate(Model = factor(Model, levels = c("Crit", "Accum", "Change")))

head(all_res)
summary(all_res)



############################################################################################
##### Analyse/summarise some of the results

## Start with proportion getting correct interaction term for each interaction strength

# Critical period setting
crit_res %>%
  group_by(Interaction) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Accumulation setting
accum_res %>%
  group_by(Interaction) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# change from time 2 to 3 setting
change_res %>%
  group_by(Interaction) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# All in one table
(intcorrect_all <- all_res %>%
  group_by(Interaction, Model) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_all, "intCorrect_allModels.csv")



###############################
### Now group by different simulation parameters (in just the 'moderate' interaction scenario)

## By sample size

# Critical period setting
crit_res %>%
  filter(Interaction == "moderate") %>%
  group_by(sampleSize) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Accumulation setting
accum_res %>%
  filter(Interaction == "moderate") %>%
  group_by(sampleSize) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Change setting
change_res %>%
  filter(Interaction == "moderate") %>%
  group_by(sampleSize) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# All in one table
(intcorrect_moderate_sampleSize <- all_res %>%
  filter(Interaction == "moderate") %>%
  group_by(sampleSize, Model) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_moderate_sampleSize, "intCorrect_moderateInt_sampleSize.csv")


## By exposure

# Critical period setting
crit_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Accumulation setting
accum_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Change setting
change_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# All in one table
(intcorrect_moderate_exposure <- all_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure, Model) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_moderate_exposure, "intCorrect_moderateInt_exposure.csv")


## By whether exposures were centered

# Critical period setting
crit_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Centered) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Accumulation setting
accum_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Centered) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Change setting
change_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Centered) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# All in one table
(intcorrect_moderate_centered <- all_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Centered, Model) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_moderate_centered, "intCorrect_moderateInt_centered.csv")


## By whether exposures were centered and type of exposure

# Critical period setting
crit_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure, Centered) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Accumulation setting
accum_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure, Centered) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Change setting
change_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure, Centered) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# All in one table
(intcorrect_moderate_exposureCentered <- all_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure, Centered, Model) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_moderate_exposureCentered, "intCorrect_moderateInt_exposureCentered.csv")


## By collinearity

# Critical period setting
crit_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Collinear) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Accumulation setting
accum_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Collinear) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Change setting
change_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Collinear) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# All in one table
(intcorrect_moderate_collinearity <- all_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Collinear, Model) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_moderate_collinearity, "intCorrect_moderateInt_collinearity.csv")


## By collinearity and whether centered

# Critical period setting
crit_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Collinear, Centered) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Accumulation setting
accum_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Collinear, Centered) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Change setting
change_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Collinear, Centered) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# All in one table
(intcorrect_moderate_collinearityCentered <- all_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Collinear, Centered, Model) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_moderate_collinearityCentered, "intCorrect_moderateInt_collinearityCentered.csv")


## By outcome

# Critical period setting
crit_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Outcome) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Accumulation setting
accum_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Outcome) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Change setting
change_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Outcome) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# All in one table
(intcorrect_moderate_outcome <- all_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Outcome, Model) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_moderate_outcome, "intCorrect_moderateInt_outcome.csv")


## By outcome and exposure

## Critical period setting
crit_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure, Outcome) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

## Accumulation setting
accum_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure, Outcome) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

## Change setting
change_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure, Outcome) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

## All in one table
(intcorrect_moderate_exposureOutcome <- all_res %>%
  filter(Interaction == "moderate") %>%
  group_by(Exposure, Outcome, Model) %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_moderate_exposureOutcome, "intCorrect_moderateInt_exposureOutcome.csv")



### Make a scatterplot matrix of all the different simulation combinations for the 'moderate' setting

## Critical period setting
temp <- crit_res %>%
  filter(Interaction == "moderate") %>%
  select(AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

diagfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping) +
    geom_histogram(binwidth = 5) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 20))
}

(crit_intcorrect <- ggpairs(temp, upper = list(continuous = wrap(ggally_cor, stars = F)),
                            diag = list(continuous = wrap(diagfun))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit axis limits
for (i in 2:crit_intcorrect$nrow) {
  for (j in 1:(i-1)) {
    crit_intcorrect[i,j] <- crit_intcorrect[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      xlim(0, 100) + ylim(0, 100)
  }
}
crit_intcorrect

pdf("crit_scatterplot_moderateInt_intCorrect.pdf", height = 6, width = 8)
crit_intcorrect
dev.off()


## Accumulation setting
temp <- accum_res %>%
  filter(Interaction == "moderate") %>%
  select(AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

diagfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping) +
    geom_histogram(binwidth = 5) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 25))
}

(accum_intcorrect <- ggpairs(temp, upper = list(continuous = wrap(ggally_cor, stars = F)),
                            diag = list(continuous = wrap(diagfun))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit axis limits
for (i in 2:accum_intcorrect$nrow) {
  for (j in 1:(i-1)) {
    accum_intcorrect[i,j] <- accum_intcorrect[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      xlim(0, 100) + ylim(0, 100)
  }
}
accum_intcorrect

pdf("accum_scatterplot_moderateInt_intCorrect.pdf", height = 6, width = 8)
accum_intcorrect
dev.off()


## Change setting
temp <- change_res %>%
  filter(Interaction == "moderate") %>%
  select(AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

diagfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping) +
    geom_histogram(binwidth = 5) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 25))
}

(change_intcorrect <- ggpairs(temp, upper = list(continuous = wrap(ggally_cor, stars = F)),
                             diag = list(continuous = wrap(diagfun))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit axis limits
for (i in 2:change_intcorrect$nrow) {
  for (j in 1:(i-1)) {
    change_intcorrect[i,j] <- change_intcorrect[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      xlim(0, 100) + ylim(0, 100)
  }
}
change_intcorrect

pdf("change_scatterplot_moderateInt_intCorrect.pdf", height = 6, width = 8)
change_intcorrect
dev.off()


## All approaches combined (Can't get the overlapping histograms by group to work correctly, so am just using density on the diagonal as easier!)
temp <- all_res %>%
  filter(Interaction == "moderate") %>%
  select(Model, AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

(all_intcorrect <- ggpairs(temp, columns = 2:5, ggplot2::aes(colour = Model),
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
        lower = list(continuous = wrap("points", alpha = 0.5)),
        upper = list(continuous = wrap(ggally_cor, stars = F))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots
for (i in 2:all_intcorrect$nrow) {
  for (j in 1:(i-1)) {
    all_intcorrect[i,j] <- all_intcorrect[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5)
  }
}

# Edit axis range of diagonals as well
for (i in 1:all_intcorrect$nrow) {
  all_intcorrect[i, i] <- all_intcorrect[i, i] + xlim(0, 100)
}
all_intcorrect

pdf("allModels_scatterplot_moderateInt_intCorrect.pdf", height = 6, width = 8)
all_intcorrect
dev.off()


### Repeat the combined plot for each of the other interaction settings

## No interaction
temp <- all_res %>%
  filter(Interaction == "None") %>%
  select(Model, AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

(all_intcorrect <- ggpairs(temp, columns = 2:5, ggplot2::aes(colour = Model),
                           diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
                           lower = list(continuous = wrap("points", alpha = 0.5)),
                           upper = list(continuous = wrap(ggally_cor, stars = F))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit the axis range
for (i in 2:all_intcorrect$nrow) {
  for (j in 1:(i-1)) {
    all_intcorrect[i,j] <- all_intcorrect[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
      xlim(0, 100) + ylim(0, 100)
  }
}

# Edit axis range of diagonals as well
for (i in 1:all_intcorrect$nrow) {
  all_intcorrect[i, i] <- all_intcorrect[i, i] + xlim(0, 100)
}

all_intcorrect

pdf("allModels_scatterplot_noInt_intCorrect.pdf", height = 6, width = 8)
all_intcorrect
dev.off()


## Very small interaction
temp <- all_res %>%
  filter(Interaction == "vSmall") %>%
  select(Model, AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

(all_intcorrect <- ggpairs(temp, columns = 2:5, ggplot2::aes(colour = Model),
                           diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
                           lower = list(continuous = wrap("points", alpha = 0.5)),
                           upper = list(continuous = wrap(ggally_cor, stars = F))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit the axis range
for (i in 2:all_intcorrect$nrow) {
  for (j in 1:(i-1)) {
    all_intcorrect[i,j] <- all_intcorrect[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
      xlim(0, 100) + ylim(0, 100)
  }
}

# Edit axis range of diagonals as well
for (i in 1:all_intcorrect$nrow) {
  all_intcorrect[i, i] <- all_intcorrect[i, i] + xlim(0, 100)
}

all_intcorrect

pdf("allModels_scatterplot_vSmallInt_intCorrect.pdf", height = 6, width = 8)
all_intcorrect
dev.off()


## Small interaction
temp <- all_res %>%
  filter(Interaction == "small") %>%
  select(Model, AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

(all_intcorrect <- ggpairs(temp, columns = 2:5, ggplot2::aes(colour = Model),
                           diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
                           lower = list(continuous = wrap("points", alpha = 0.5)),
                           upper = list(continuous = wrap(ggally_cor, stars = F))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit the axis range
for (i in 2:all_intcorrect$nrow) {
  for (j in 1:(i-1)) {
    all_intcorrect[i,j] <- all_intcorrect[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
      xlim(0, 100) + ylim(0, 100)
  }
}

# Edit axis range of diagonals as well
for (i in 1:all_intcorrect$nrow) {
  all_intcorrect[i, i] <- all_intcorrect[i, i] + xlim(0, 100)
}

all_intcorrect

pdf("allModels_scatterplot_smallInt_intCorrect.pdf", height = 6, width = 8)
all_intcorrect
dev.off()


## Large interaction
temp <- all_res %>%
  filter(Interaction == "large") %>%
  select(Model, AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

(all_intcorrect <- ggpairs(temp, columns = 2:5, ggplot2::aes(colour = Model),
                           diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
                           lower = list(continuous = wrap("points", alpha = 0.5)),
                           upper = list(continuous = wrap(ggally_cor, stars = F))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit the axis range
for (i in 2:all_intcorrect$nrow) {
  for (j in 1:(i-1)) {
    all_intcorrect[i,j] <- all_intcorrect[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
      xlim(0, 100) + ylim(0, 100)
  }
}

# Edit axis range of diagonals as well
for (i in 1:all_intcorrect$nrow) {
  all_intcorrect[i, i] <- all_intcorrect[i, i] + xlim(0, 100)
}

all_intcorrect

pdf("allModels_scatterplot_largeInt_intCorrect.pdf", height = 6, width = 8)
all_intcorrect
dev.off()


## Very Large interaction
temp <- all_res %>%
  filter(Interaction == "vLarge") %>%
  select(Model, AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

(all_intcorrect <- ggpairs(temp, columns = 2:5, ggplot2::aes(colour = Model),
                           diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
                           lower = list(continuous = wrap("points", alpha = 0.5)),
                           upper = list(continuous = wrap(ggally_cor, stars = F))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit the axis range
for (i in 2:all_intcorrect$nrow) {
  for (j in 1:(i-1)) {
    all_intcorrect[i,j] <- all_intcorrect[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
      xlim(0, 100) + ylim(0, 100)
  }
}

# Edit axis range of diagonals as well
for (i in 1:all_intcorrect$nrow) {
  all_intcorrect[i, i] <- all_intcorrect[i, i] + xlim(0, 100)
}

all_intcorrect

pdf("allModels_scatterplot_vLargeInt_intCorrect.pdf", height = 6, width = 8)
all_intcorrect
dev.off()



############################
##### Also want to make an equivalent plot for the scenario where all of the 'change' hypotheses were removed, to see if this helped select the correct model (obv., this would only apply to the 'critical period' and 'accumulation' settings)

### Read in the data

## Critical period - Moderate interaction - No Change hypotheses
crit_moderateInt_noChange <- read_csv("simulationResults_crit_moderateInt_noChange.csv")
head(crit_moderateInt_noChange)
summary(crit_moderateInt_noChange)

crit_moderateInt_noChange <- crit_moderateInt_noChange %>%
  mutate(Model = "Crit", Interaction = "moderate") %>%
  relocate(c(Model, Interaction))

## Accumulation - Moderate interaction - No Change hypotheses
accum_moderateInt_noChange <- read_csv("simulationResults_accum_moderateInt_noChange.csv")
head(accum_moderateInt_noChange)
summary(accum_moderateInt_noChange)

accum_moderateInt_noChange <- accum_moderateInt_noChange %>%
  mutate(Model = "Accum", Interaction = "moderate") %>%
  relocate(c(Model, Interaction))


## Combine these into a single 'critical period' dataset
noChange_res <- rbind(crit_moderateInt_noChange, accum_moderateInt_noChange)
head(noChange_res)
summary(noChange_res)

# Recode the string vars as factors
noChange_res <- noChange_res %>%
  mutate(Model = factor(Model, levels = c("Crit", "Accum"))) %>%
  mutate(Interaction = as.factor(Interaction)) %>%
  mutate(Exposure = factor(Exposure, levels = c("Binary", "Cont"))) %>%
  mutate(Centered = factor(Centered, levels = c("No", "Yes"))) %>%
  mutate(Collinear = factor(Collinear, levels = c("Low", "High"))) %>%
  mutate(Outcome = factor(Outcome, levels = c("Cont", "Binary")))

head(noChange_res)
summary(noChange_res)


## Start with proportion getting correct interaction term for each interaction strength

# Critical period setting
crit_moderateInt_noChange %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# Accumulation setting
accum_moderateInt_noChange %>%
  summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
            min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
            mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
            min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
            mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
            min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
            mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
            min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect))

# All in one table
(intcorrect_mod_noChange <- noChange_res %>%
    group_by(Model) %>%
    summarise(mean_AIC = mean(AIC_intcorrect), median_AIC = median(AIC_intcorrect),
              min_AIC = min(AIC_intcorrect), max_AIC = max(AIC_intcorrect),
              mean_BIC = mean(BIC_intcorrect), median_BIC = median(BIC_intcorrect),
              min_BIC = min(BIC_intcorrect), max_BIC_int = max(BIC_intcorrect),
              mean_CV_1SE = mean(CV_1SE_intcorrect), median_CV_1SE = median(CV_1SE_intcorrect),
              min_CV_1SE = min(CV_1SE_intcorrect), max_CV_1SE = max(CV_1SE_intcorrect),
              mean_CV_minMSE = mean(CV_minMSE_intcorrect), median_CV_minMSE = median(CV_minMSE_intcorrect),
              min_CV_minMSE = min(CV_minMSE_intcorrect), max_CV_minMSE = max(CV_minMSE_intcorrect)))

# Save this table
write_csv(intcorrect_mod_noChange, "intCorrect_moderateInt_noChange.csv")


### Make a scatterplot matrix of all the different simulation combinations for the 'moderate' setting

## Critical period setting
temp <- crit_moderateInt_noChange %>%
  select(AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

diagfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping) +
    geom_histogram(binwidth = 5) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 20))
}

(crit_intcorrect_noChange <- ggpairs(temp, upper = list(continuous = wrap(ggally_cor, stars = F)),
                            diag = list(continuous = wrap(diagfun))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit axis limits
for (i in 2:crit_intcorrect_noChange$nrow) {
  for (j in 1:(i-1)) {
    crit_intcorrect_noChange[i,j] <- crit_intcorrect_noChange[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      xlim(0, 100) + ylim(0, 100)
  }
}
crit_intcorrect_noChange

pdf("crit_scatterplot_moderateInt_intCorrect_noChange.pdf", height = 6, width = 8)
crit_intcorrect_noChange
dev.off()


## Accumulation setting
temp <- accum_moderateInt_noChange %>%
  select(AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

diagfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping) +
    geom_histogram(binwidth = 5) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 25))
}

(accum_intcorrect_noChange <- ggpairs(temp, upper = list(continuous = wrap(ggally_cor, stars = F)),
                             diag = list(continuous = wrap(diagfun))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots and edit axis limits
for (i in 2:accum_intcorrect_noChange$nrow) {
  for (j in 1:(i-1)) {
    accum_intcorrect_noChange[i,j] <- accum_intcorrect_noChange[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      xlim(0, 100) + ylim(0, 100)
  }
}
accum_intcorrect_noChange

pdf("accum_scatterplot_moderateInt_intCorrect_noChange.pdf", height = 6, width = 8)
accum_intcorrect_noChange
dev.off()


## All approaches combined (Can't get the overlapping histograms by group to work correctly, so am just using density on the diagonal as easier!)
temp <- noChange_res %>%
  select(Model, AIC_intcorrect, BIC_intcorrect, CV_1SE_intcorrect, CV_minMSE_intcorrect) %>%
  rename(AIC = AIC_intcorrect, BIC = BIC_intcorrect, CV_1SE = CV_1SE_intcorrect, CV_minMSE = CV_minMSE_intcorrect)
head(temp)

(all_intcorrect_noChange <- ggpairs(temp, columns = 2:5, ggplot2::aes(colour = Model),
                           diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
                           lower = list(continuous = wrap("points", alpha = 0.5)),
                           upper = list(continuous = wrap(ggally_cor, stars = F))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# Add AB line to each of the lower plots
for (i in 2:all_intcorrect_noChange$nrow) {
  for (j in 1:(i-1)) {
    all_intcorrect_noChange[i,j] <- all_intcorrect_noChange[i,j] + 
      geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5)
  }
}

# Edit axis range of diagonals as well
for (i in 1:all_intcorrect_noChange$nrow) {
  all_intcorrect_noChange[i, i] <- all_intcorrect_noChange[i, i] + xlim(0, 100)
}
all_intcorrect_noChange

pdf("allModels_scatterplot_moderateInt_intCorrect_noChange.pdf", height = 6, width = 8)
all_intcorrect_noChange
dev.off()



################################
### Also look at whether the different methods detected the true model

# Critical period setting
crit_res %>%
  group_by(Interaction) %>%
  summarise(mean_AIC = mean(AIC_propcorrect), median_AIC = median(AIC_propcorrect),
            min_AIC = min(AIC_propcorrect), max_AIC = max(AIC_propcorrect),
            mean_BIC = mean(BIC_propcorrect), median_BIC = median(BIC_propcorrect),
            min_BIC = min(BIC_propcorrect), max_BIC_int = max(BIC_propcorrect),
            mean_CV_1SE = mean(CV_1SE_propcorrect), median_CV_1SE = median(CV_1SE_propcorrect),
            min_CV_1SE = min(CV_1SE_propcorrect), max_CV_1SE = max(CV_1SE_propcorrect),
            mean_CV_minMSE = mean(CV_minMSE_propcorrect), median_CV_minMSE = median(CV_minMSE_propcorrect),
            min_CV_minMSE = min(CV_minMSE_propcorrect), max_CV_minMSE = max(CV_minMSE_propcorrect))

# Accumulation setting
accum_res %>%
  group_by(Interaction) %>%
  summarise(mean_AIC = mean(AIC_propcorrect), median_AIC = median(AIC_propcorrect),
            min_AIC = min(AIC_propcorrect), max_AIC = max(AIC_propcorrect),
            mean_BIC = mean(BIC_propcorrect), median_BIC = median(BIC_propcorrect),
            min_BIC = min(BIC_propcorrect), max_BIC_int = max(BIC_propcorrect),
            mean_CV_1SE = mean(CV_1SE_propcorrect), median_CV_1SE = median(CV_1SE_propcorrect),
            min_CV_1SE = min(CV_1SE_propcorrect), max_CV_1SE = max(CV_1SE_propcorrect),
            mean_CV_minMSE = mean(CV_minMSE_propcorrect), median_CV_minMSE = median(CV_minMSE_propcorrect),
            min_CV_minMSE = min(CV_minMSE_propcorrect), max_CV_minMSE = max(CV_minMSE_propcorrect))

# Change setting
change_res %>%
  group_by(Interaction) %>%
  summarise(mean_AIC = mean(AIC_propcorrect), median_AIC = median(AIC_propcorrect),
            min_AIC = min(AIC_propcorrect), max_AIC = max(AIC_propcorrect),
            mean_BIC = mean(BIC_propcorrect), median_BIC = median(BIC_propcorrect),
            min_BIC = min(BIC_propcorrect), max_BIC_int = max(BIC_propcorrect),
            mean_CV_1SE = mean(CV_1SE_propcorrect), median_CV_1SE = median(CV_1SE_propcorrect),
            min_CV_1SE = min(CV_1SE_propcorrect), max_CV_1SE = max(CV_1SE_propcorrect),
            mean_CV_minMSE = mean(CV_minMSE_propcorrect), median_CV_minMSE = median(CV_minMSE_propcorrect),
            min_CV_minMSE = min(CV_minMSE_propcorrect), max_CV_minMSE = max(CV_minMSE_propcorrect))

# All in one table
(truemodel_all <- all_res %>%
  group_by(Interaction, Model) %>%
  summarise(mean_AIC = mean(AIC_propcorrect), median_AIC = median(AIC_propcorrect),
            min_AIC = min(AIC_propcorrect), max_AIC = max(AIC_propcorrect),
            mean_BIC = mean(BIC_propcorrect), median_BIC = median(BIC_propcorrect),
            min_BIC = min(BIC_propcorrect), max_BIC_int = max(BIC_propcorrect),
            mean_CV_1SE = mean(CV_1SE_propcorrect), median_CV_1SE = median(CV_1SE_propcorrect),
            min_CV_1SE = min(CV_1SE_propcorrect), max_CV_1SE = max(CV_1SE_propcorrect),
            mean_CV_minMSE = mean(CV_minMSE_propcorrect), median_CV_minMSE = median(CV_minMSE_propcorrect),
            min_CV_minMSE = min(CV_minMSE_propcorrect), max_CV_minMSE = max(CV_minMSE_propcorrect)))

# Save this table
write_csv(truemodel_all, "trueModel_allModels.csv")


## And whether the different methods detected the true model plus additional terms

# Critical period setting
crit_res %>%
  group_by(Interaction) %>%
  summarise(mean_AIC = mean(AIC_mainintextra), median_AIC = median(AIC_mainintextra),
            min_AIC = min(AIC_mainintextra), max_AIC = max(AIC_mainintextra),
            mean_BIC = mean(BIC_mainintextra), median_BIC = median(BIC_mainintextra),
            min_BIC = min(BIC_mainintextra), max_BIC_int = max(BIC_mainintextra),
            mean_CV_1SE = mean(CV_1SE_mainintextra), median_CV_1SE = median(CV_1SE_mainintextra),
            min_CV_1SE = min(CV_1SE_mainintextra), max_CV_1SE = max(CV_1SE_mainintextra),
            mean_CV_minMSE = mean(CV_minMSE_mainintextra), median_CV_minMSE = median(CV_minMSE_mainintextra),
            min_CV_minMSE = min(CV_minMSE_mainintextra), max_CV_minMSE = max(CV_minMSE_mainintextra))

# Accumulation setting
accum_res %>%
  group_by(Interaction) %>%
  summarise(mean_AIC = mean(AIC_mainintextra), median_AIC = median(AIC_mainintextra),
            min_AIC = min(AIC_mainintextra), max_AIC = max(AIC_mainintextra),
            mean_BIC = mean(BIC_mainintextra), median_BIC = median(BIC_mainintextra),
            min_BIC = min(BIC_mainintextra), max_BIC_int = max(BIC_mainintextra),
            mean_CV_1SE = mean(CV_1SE_mainintextra), median_CV_1SE = median(CV_1SE_mainintextra),
            min_CV_1SE = min(CV_1SE_mainintextra), max_CV_1SE = max(CV_1SE_mainintextra),
            mean_CV_minMSE = mean(CV_minMSE_mainintextra), median_CV_minMSE = median(CV_minMSE_mainintextra),
            min_CV_minMSE = min(CV_minMSE_mainintextra), max_CV_minMSE = max(CV_minMSE_mainintextra))

# Change setting
change_res %>%
  group_by(Interaction) %>%
  summarise(mean_AIC = mean(AIC_mainintextra), median_AIC = median(AIC_mainintextra),
            min_AIC = min(AIC_mainintextra), max_AIC = max(AIC_mainintextra),
            mean_BIC = mean(BIC_mainintextra), median_BIC = median(BIC_mainintextra),
            min_BIC = min(BIC_mainintextra), max_BIC_int = max(BIC_mainintextra),
            mean_CV_1SE = mean(CV_1SE_mainintextra), median_CV_1SE = median(CV_1SE_mainintextra),
            min_CV_1SE = min(CV_1SE_mainintextra), max_CV_1SE = max(CV_1SE_mainintextra),
            mean_CV_minMSE = mean(CV_minMSE_mainintextra), median_CV_minMSE = median(CV_minMSE_mainintextra),
            min_CV_minMSE = min(CV_minMSE_mainintextra), max_CV_minMSE = max(CV_minMSE_mainintextra))

# All in one table
(trueplusextra_all <- all_res %>%
  group_by(Interaction, Model) %>%
  summarise(mean_AIC = mean(AIC_mainintextra), median_AIC = median(AIC_mainintextra),
            min_AIC = min(AIC_mainintextra), max_AIC = max(AIC_mainintextra),
            mean_BIC = mean(BIC_mainintextra), median_BIC = median(BIC_mainintextra),
            min_BIC = min(BIC_mainintextra), max_BIC_int = max(BIC_mainintextra),
            mean_CV_1SE = mean(CV_1SE_mainintextra), median_CV_1SE = median(CV_1SE_mainintextra),
            min_CV_1SE = min(CV_1SE_mainintextra), max_CV_1SE = max(CV_1SE_mainintextra),
            mean_CV_minMSE = mean(CV_minMSE_mainintextra), median_CV_minMSE = median(CV_minMSE_mainintextra),
            min_CV_minMSE = min(CV_minMSE_mainintextra), max_CV_minMSE = max(CV_minMSE_mainintextra)))

# Save this table
write_csv(trueplusextra_all, "trueModelPlusExtra_allModels.csv")



