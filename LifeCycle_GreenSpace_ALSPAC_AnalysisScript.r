### LifeCycle Green Space Analysis - ALSPAC analysis script

## Created 27/1/2022 by Dan Major-Smith
## R version 4.0.4


### Note that originally planned to look at all combinations of green space exposures (4 - Access to nearby green space [binary], distance to nearest green space [cont], NDVI value [cont] and access to garden [binary]), outcomes (4 - BMI [cont], overweight [binary], systolic BP [cont] and diastolic BP [cont]) and SEP variables/interaction terms (3 - maternal education, area-level deprivation and household income [all binary]) in a multiverse type analyses. However, as this paper has become more methodological and less applied, have decided just to focus on 1 green space exposure (access to nearby green space), 1 SEP interaction term (maternal education) and the 4 cardiovascular outcomes. This simplifies the analysis and the presentation of this illustrative example (earlier versions of this script which looked at all possible combinations found little evidence for an association between any of the access to green space exposures and any of the cardiovascular outcomes, for any of the SEP interaction terms - so the shortened results presented here are representative of the wider analysis)


## Make sure workspace is clear
rm(list = ls())

# Install and load packages
#install.packages("tidyverse")
library(tidyverse)

#install.packages("glmnet")
library(glmnet)

#install.packages("dagitty")
library(dagitty)

## For the sankey plots (also know as 'alluvial' plots, will use the package 'ggalluvial' - see: https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html)
#install.packages("ggalluvial")
library(ggalluvial)

#install.packages("selectiveInference")
library(selectiveInference)


# Set working directory
setwd("C:\\Users\\ds16565\\OneDrive - University of Bristol\\MyFiles-Migrated\\Documents\\Projects\\Lifecycle\\Results")


# Set the seed, so all cross-validated lasso models give same results
set.seed(3930)


### Create DAG to show hypothesised causal associations between variables (haven't included green space by SEP interactions for all exposure time-points, as too messy to visualise; arrows from ethnicity and SEP to all green space variables have also been omitted)
dag <- dagitty('dag {
                SEP [pos = "1,1"]
                GreenPreg [pos = "0,2"]
                Green4 [pos = "0,3"]
                Green7 [pos = "0,4"]
                SEP_GreenPreg_int [pos = "1,2"]
                cardio [pos = "2,3"]
                ethnicity [pos = "1,0"]
                age [pos = "1.85,1"]
                sex [pos = "2.15,1"]
                
                SEP -> GreenPreg
                SEP -> cardio
                GreenPreg -> Green4
                Green4 -> Green7
                GreenPreg -> SEP_GreenPreg_int
                SEP -> SEP_GreenPreg_int
                SEP_GreenPreg_int -> cardio
                GreenPreg -> cardio
                Green4 -> cardio
                Green7 -> cardio
                ethnicity -> GreenPreg
                ethnicity -> SEP
                ethnicity -> cardio
                age -> cardio
                sex -> cardio
                }')
plot(dag)


## Will also make a simplified DAG, with all of the green space variables collapsed together
dag_simp <- dagitty('dag {
                SEP [pos = "1,1"]
                GreenSpace [pos = "0,2"]
                SEP_GreenSpace_int [pos = "1,1.5"]
                Cardio [pos = "2,2"]
                Ethnicity [pos = "1,0.5"]
                Age [pos = "1.85,1"]
                Sex [pos = "2.15,1"]
                
                SEP -> GreenSpace
                SEP -> Cardio
                GreenSpace -> SEP_GreenSpace_int
                SEP -> SEP_GreenSpace_int
                SEP_GreenSpace_int -> Cardio
                GreenSpace -> Cardio
                Ethnicity -> GreenSpace
                Ethnicity -> SEP
                Ethnicity -> Cardio
                Age -> Cardio
                Sex -> Cardio
                }')
plot(dag_simp)

# Save this DAG
pdf("GreenSpaceDAG.pdf", height = 4, width = 8)
plot(dag_simp)
dev.off()



### Will also initialise a function to turn lasso output into a useful summary table
lasso_table <- function(lasso_model) {
  old_covars <- ""
  old_deviance <- 0
  old_varNum <- 0
  old_lambda <- NA
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  #df
  
  for (i in 1:length(lasso_model$df)) {
    #print(i)
    new_covars <- attributes(which(lasso_model$beta[, i] != 0))$names
    new_deviance <- lasso_model$dev.ratio[i]
    new_varNum <- lasso_model$df[i]
    new_lambda <- lasso_model$lambda[i]
    #print(new_covars); print(new_deviance); print(new_varNum)
    
    # See if covars added
    if (new_varNum > old_varNum) {
      change <- setdiff(new_covars, old_covars) # Find the new covariate(s) added
      change <- paste(change, collapse = " ") # Combine variable together, if > 1
      change <- paste0(change, " (+)") # Append a "(+)" sign
      change_drop <- setdiff(old_covars, new_covars) # Make sure no covariates dropped at same time
      if (!identical(change_drop, character(0))) { # If a covar also dropped, then combine with 'change'
        change_drop <- paste(change_drop, collapse = " ") # Combine variable together, if > 1
        change <- paste0(change, " ", change_drop, " (-)")
      }
      dev_diff <- round((new_deviance - old_deviance) * 100, 3) # Diff in deviance between current and previous lambda
      new_dev <- round(new_deviance * 100, 3) # Current deviance value
      vars_noSpaces <- paste(new_covars, collapse = " ") # Combine all variables together
      temp <- cbind(change, new_dev, dev_diff, new_varNum, new_lambda, vars_noSpaces, NA, NA) # Combine values together
      df <- rbind(df, temp) # Merge with template data frame
    }
    
    # See if covars removed
    if (new_varNum < old_varNum) {
      change <- setdiff(old_covars, new_covars) # Find the covariate(s) removed
      change <- paste(change, collapse = " ") # Combine variable together, if > 1
      change <- paste0(change, " (-)") # Append a "(-)" sign
      change_add <- setdiff(new_covars, old_covars) # Make sure no covariates added at same time
      if (!identical(change_add, character(0))) { # If a covar also dropped, then combine with 'change'
        change_add <- paste(change_add, collapse = " ") # Combine variable together, if > 1
        change <- paste0(change, " ", change_add, " (+)")
      }
      dev_diff <- round((new_deviance - old_deviance) * 100, 3) # Diff in deviance between current and previous lambda
      new_dev <- round(new_deviance * 100, 3) # Current deviance value
      vars_noSpaces <- paste(new_covars, collapse = " ") # Combine all variables together
      temp <- cbind(change, new_dev, dev_diff, new_varNum, new_lambda, vars_noSpaces, NA, NA) # Combine values together
      df <- rbind(df, temp) # Merge with template data frame
    }
    
    # See if covars added and removed at the same time (where number of variables stays the same)
    if (new_varNum == old_varNum & setequal(old_covars, new_covars) == FALSE) {
      change_add <- setdiff(new_covars, old_covars) # Find the covariate(s) added
      change_add <- paste(change_add, collapse = " ") # Combine variables together, if > 1
      change_add <- paste0(change_add, " (+)") # Append a "(+)" sign
      change_drop <- setdiff(old_covars, new_covars) # Find the covariate(s) removed
      change_drop <- paste(change_drop, collapse = " ") # Combine variables together, if > 1
      change_drop <- paste0(change_drop, " (-)") # Append a "(-)" sign
      change <- paste0(change_add, " ", change_drop) # Combine the added and dropped variables together
      dev_diff <- round((new_deviance - old_deviance) * 100, 3) # Diff in deviance between current and previous lambda
      new_dev <- round(new_deviance * 100, 3) # Current deviance value
      vars_noSpaces <- paste(new_covars, collapse = " ") # Combine all variables together
      temp <- cbind(change, new_dev, dev_diff, new_varNum, new_lambda, vars_noSpaces, NA, NA) # Combine values together
      df <- rbind(df, temp) # Merge with template data frame
    }
    
    # Rename the old covars, deviance and variable number
    old_covars <- new_covars
    old_deviance <- new_deviance
    old_varNum <- new_varNum
  }
  
  colnames(df) <- c("Variables", "DevRatio", "DevDiff", "VarNum", "Lambda", "model_vars", "aic", "bic")
  #df
  
  # Make a var to show number of steps where variables added, and rename the covars to blank (as are included by default)
  df$steps <- 1:nrow(df)
  df$Variables[df$steps == 1] <- ""
  return(df)
}



####################################################################################################
##### Read in data
data_raw <- read_csv(file = "../Data/B3930_lifeCycleAndALSPAC.csv")

## Keep just relevant variables
data <- data_raw %>%
  dplyr::select(aln, qlet, 
         f7ms026a, f7ms102, f7sa021, f7sa022, 
         greenSpace_preg, greenSpace_4, greenSpace_7,
         greenDist_preg, greenDist_4, greenDist_7,
         NDVI500_preg, NDVI500_4, NDVI500_7,
         a049, h298, m2083,
         c645a, c666a, edu_m_0, eusilc_income, eusilc_income_quintiles, jan1993imd2010q5_M, areaSES_preg,
         kz021, f7003c, c800, ethn3_m)

head(data)
glimpse(data)
summary(data)


#### Clean/process variables as necessary

### Start with cardiovascular outcomes

## BMI of child at F@7
summary(data$f7ms026a)

data <- data %>%
  rename(BMI_f7 = f7ms026a) %>%
  mutate(BMI_f7 = ifelse(BMI_f7 < 0, NA, BMI_f7))

summary(data$BMI_f7)
sd(data$BMI_f7, na.rm = TRUE)
hist(data$BMI_f7)


## Create binary 'overweight/obese' variable - Using standard deviation score from 1990 British Growth reference charts (i.e., adjusted for age and sex), and defining overweight/obese as at or above 85% centile (z-score of 1.036)
summary(data$f7ms102)

data <- data %>%
  rename(BMI_f7_z = f7ms102) %>%
  mutate(BMI_f7_z = ifelse(BMI_f7_z < -100, NA, BMI_f7_z))

summary(data$BMI_f7_z)
hist(data$BMI_f7_z)

data <- data %>%
  mutate(overweight = ifelse(is.na(BMI_f7_z), NA,
                                   ifelse(BMI_f7_z >= qnorm(0.85), 1, 0))) %>%
  mutate(overweight = as.factor(overweight)) %>%
  mutate(overweight = fct_recode(overweight, "Yes" = "1", "No" = "0"))

table(data$overweight, useNA = "ifany")
prop.table(table(data$overweight)) * 100
prop.table(table(data$overweight, useNA = "ifany")) * 100

## Systolic and diastolic blood pressure
summary(data$f7sa021); summary(data$f7sa022)

data <- data %>%
  rename(sysBP = f7sa021, diaBP = f7sa022) %>%
  mutate(sysBP = ifelse(sysBP < 0, NA, sysBP)) %>%
  mutate(diaBP = ifelse(diaBP < 0, NA, diaBP))

summary(data$sysBP); summary(data$diaBP)
sd(data$sysBP, na.rm = TRUE); sd(data$diaBP, na.rm = TRUE)
hist(data$sysBP)
hist(data$diaBP)


### Exposure variables

## Access to green space (with 300m of address) in pregnancy, age 4 and age 7 - First tidy variables, then explore variation over time
table(data$greenSpace_preg, useNA = "ifany")
table(data$greenSpace_4, useNA = "ifany")
table(data$greenSpace_7, useNA = "ifany")

# code as NA if any time-points are NA
data <- data %>%
  mutate(cc = ifelse(complete.cases(greenSpace_preg, greenSpace_4, greenSpace_7), 1, 0)) %>%
  mutate(greenSpace_preg = ifelse(cc == 0, NA, greenSpace_preg)) %>%
  mutate(greenSpace_4 = ifelse(cc == 0, NA, greenSpace_4)) %>%
  mutate(greenSpace_7 = ifelse(cc == 0, NA, greenSpace_7)) %>%
  dplyr::select(-cc)

table(data$greenSpace_preg, useNA = "ifany")
prop.table(table(data$greenSpace_preg)) * 100
prop.table(table(data$greenSpace_preg, useNA = "ifany")) * 100

table(data$greenSpace_4, useNA = "ifany")
prop.table(table(data$greenSpace_4)) * 100
prop.table(table(data$greenSpace_4, useNA = "ifany")) * 100

table(data$greenSpace_7, useNA = "ifany")
prop.table(table(data$greenSpace_7)) * 100
prop.table(table(data$greenSpace_7, useNA = "ifany")) * 100

# Create pattern of green space access over time at all 3 time-points
data <- data %>%
  mutate(greenSpace_combo = ifelse(!complete.cases(greenSpace_preg, greenSpace_4, greenSpace_7), NA, 
                                   ifelse(greenSpace_preg == 0 & greenSpace_4 == 0 & greenSpace_7 == 0, 1,
                                   ifelse(greenSpace_preg == 1 & greenSpace_4 == 0 & greenSpace_7 == 0, 2, 
                                   ifelse(greenSpace_preg == 0 & greenSpace_4 == 1 & greenSpace_7 == 0, 3, 
                                   ifelse(greenSpace_preg == 0 & greenSpace_4 == 0 & greenSpace_7 == 1, 4, 
                                   ifelse(greenSpace_preg == 1 & greenSpace_4 == 1 & greenSpace_7 == 0, 5, 
                                   ifelse(greenSpace_preg == 1 & greenSpace_4 == 0 & greenSpace_7 == 1, 6, 
                                   ifelse(greenSpace_preg == 0 & greenSpace_4 == 1 & greenSpace_7 == 1, 7, 8)))))))))

table(data$greenSpace_combo, useNA = "ifany")
prop.table(table(data$greenSpace_combo)) * 100

table(data$greenSpace_combo[!is.na(data$BMI_f7)])

## Turn this date into a sankey plot - To do this, first convert the data to summary wide format, then convert to long/lode format.
data_temp_green <- data %>%
  dplyr::select(greenSpace_preg, greenSpace_4, greenSpace_7) %>%
  filter(complete.cases(greenSpace_preg, greenSpace_4, greenSpace_7)) %>%
  mutate(greenSpace_preg = as.factor(greenSpace_preg)) %>%
  mutate(greenSpace_preg = recode(greenSpace_preg, "0" = "No", "1" = "Yes")) %>%
  mutate(greenSpace_preg = factor(greenSpace_preg, levels = c("Yes", "No"))) %>%
  mutate(greenSpace_4 = as.factor(greenSpace_4)) %>%
  mutate(greenSpace_4 = recode(greenSpace_4, "0" = "No", "1" = "Yes")) %>%
  mutate(greenSpace_4 = factor(greenSpace_4, levels = c("Yes", "No"))) %>%
  mutate(greenSpace_7 = as.factor(greenSpace_7)) %>%
  mutate(greenSpace_7 = recode(greenSpace_7, "0" = "No", "1" = "Yes")) %>%
  mutate(greenSpace_7 = factor(greenSpace_7, levels = c("Yes", "No"))) %>%
  group_by(greenSpace_preg, greenSpace_4, greenSpace_7) %>%
  summarise(freq = n())

summary(data_temp_green)
head(data_temp_green)

data_temp_lodes_green <- to_lodes_form(data_temp_green, axes = 1:3, id = "traj")
data_temp_lodes_green <- data_temp_lodes_green %>%
  rename(time = x, Response = stratum)
head(data_temp_lodes_green)
summary(data_temp_lodes_green)

sankey_green <- ggplot(data_temp_lodes_green,
                      aes(x = time, stratum = Response, alluvium = traj,
                          y = freq,
                          fill = Response, label = Response)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("Time-point") + ylab("Frequency")

sankey_green

# Save this plot
pdf("sankey_green.pdf", height = 6, width = 10)
plot(sankey_green)
dev.off()


# All possible combinations have some cell counts, but some are very small (e.g., only 22 people have the pattern 0-1-0) - This drops even more (to 10) if restrict to children with F@7 BMI data... This *might* just about be enough to work with, but will also use just 2 time-points (pregnancy and age 7) to ensure sufficient sample sizes and see if get same pattern of results
data <- data %>%
  mutate(greenSpace_combo2 = ifelse(!complete.cases(greenSpace_preg, greenSpace_7), NA,
                                    ifelse(greenSpace_preg == 0 & greenSpace_7 == 0, 1,
                                    ifelse(greenSpace_preg == 1 & greenSpace_7 == 0, 2,
                                    ifelse(greenSpace_preg == 0 & greenSpace_7 == 1, 3, 4)))))

table(data$greenSpace_combo2, useNA = "ifany")
prop.table(table(data$greenSpace_combo2)) * 100

table(data$greenSpace_combo2[!is.na(data$BMI_f7)])

# Sample sizes much better now (min 372 in F@7 sample), but as only 2 time-points can't explore as many life-course trajectories



## Distance to green space in pregnancy, age 4 and age 7 - First tidy variables, then explore variation over time
summary(data$greenDist_preg); summary(data$greenDist_4); summary(data$greenDist_7)

# code as NA if any time-points are NA
data <- data %>%
  mutate(cc = ifelse(complete.cases(greenDist_preg, greenDist_4, greenDist_7), 1, 0)) %>%
  mutate(greenDist_preg = ifelse(cc == 0, NA, greenDist_preg)) %>%
  mutate(greenDist_4 = ifelse(cc == 0, NA, greenDist_4)) %>%
  mutate(greenDist_7 = ifelse(cc == 0, NA, greenDist_7)) %>%
  dplyr::select(-cc)

summary(data$greenDist_preg); summary(data$greenDist_4); summary(data$greenDist_7)
sd(data$greenDist_preg, na.rm = TRUE); sd(data$greenDist_4, na.rm = TRUE); sd(data$greenDist_7, na.rm = TRUE)
hist(data$greenDist_preg)
hist(data$greenDist_4)
hist(data$greenDist_7)

# How correlated are all these exposures - They are correlated, but not perfectly, so hopefully is enough variation to work with
cor(data$greenDist_preg, data$greenDist_4, use = "pairwise.complete.obs")
cor(data$greenDist_preg, data$greenDist_7, use = "pairwise.complete.obs")
cor(data$greenDist_4, data$greenDist_7, use = "pairwise.complete.obs")



## Next exposure is NDVI value within 500m of address - First tidy variables, then explore variation over time
summary(data$NDVI500_preg); summary(data$NDVI500_4); summary(data$NDVI500_7)

# code as NA if any time-points are NA
data <- data %>%
  mutate(cc = ifelse(complete.cases(NDVI500_preg, NDVI500_4, NDVI500_7), 1, 0)) %>%
  mutate(NDVI500_preg = ifelse(cc == 0, NA, NDVI500_preg)) %>%
  mutate(NDVI500_4 = ifelse(cc == 0, NA, NDVI500_4)) %>%
  mutate(NDVI500_7 = ifelse(cc == 0, NA, NDVI500_7)) %>%
  dplyr::select(-cc)

summary(data$NDVI500_preg); summary(data$NDVI500_4); summary(data$NDVI500_7)
sd(data$NDVI500_preg, na.rm = TRUE); sd(data$NDVI500_4, na.rm = TRUE); sd(data$NDVI500_7, na.rm = TRUE)
barplot(table(data$NDVI500_preg))
barplot(table(data$NDVI500_4))
barplot(table(data$NDVI500_7))

# How correlated are all these exposures - They are correlated, but not perfectly, so hopefully is enough variation to work with
cor(data$NDVI500_preg, data$NDVI500_4, use = "pairwise.complete.obs")
cor(data$NDVI500_preg, data$NDVI500_7, use = "pairwise.complete.obs")
cor(data$NDVI500_4, data$NDVI500_7, use = "pairwise.complete.obs")



## Access to garden in pregnancy, age 4 and age 7 (based on ALSPAC, not LifeCycle, data) - First tidy variables, then explore variation over time
table(data$a049, useNA = "ifany")
table(data$h298, useNA = "ifany")
table(data$m2083, useNA = "ifany")

# code as NA if any time-points are NA (and tidy vars more generally)
data <- data %>%
  rename(garden_preg = a049, garden_4 = h298, garden_7 = m2083) %>%
  mutate(garden_preg = ifelse(garden_preg < 0, NA, garden_preg)) %>%
  mutate(garden_4 = ifelse(garden_4 < 0, NA, garden_4)) %>%
  mutate(garden_7 = ifelse(garden_7 < 0, NA, garden_7)) %>%
  mutate(cc = ifelse(complete.cases(garden_preg, garden_4, garden_7), 1, 0)) %>%
  mutate(garden_preg = ifelse(cc == 0, NA, garden_preg)) %>%
  mutate(garden_4 = ifelse(cc == 0, NA, garden_4)) %>%
  mutate(garden_7 = ifelse(cc == 0, NA, garden_7)) %>%
  dplyr::select(-cc)

table(data$garden_preg, useNA = "ifany")
table(data$garden_4, useNA = "ifany")
table(data$garden_7, useNA = "ifany")

# Recode to binary 'access to garden' or not
data <- data %>%
  mutate(garden_preg = ifelse(garden_preg == 2, 1,
                              ifelse(garden_preg == 3, 0, garden_preg))) %>%
  mutate(garden_4 = ifelse(garden_4 == 2, 1,
                              ifelse(garden_4 == 3, 0, garden_4))) %>%
  mutate(garden_7 = ifelse(garden_7 == 2, 1,
                              ifelse(garden_7 == 3, 0, garden_7)))

table(data$garden_preg, useNA = "ifany")
prop.table(table(data$garden_preg)) * 100
prop.table(table(data$garden_preg, useNA = "ifany")) * 100

table(data$garden_4, useNA = "ifany")
prop.table(table(data$garden_4)) * 100
prop.table(table(data$garden_4, useNA = "ifany")) * 100

table(data$garden_7, useNA = "ifany")
prop.table(table(data$garden_7)) * 100
prop.table(table(data$garden_7, useNA = "ifany")) * 100

# Create pattern of green space access over time at all 3 time-points
data <- data %>%
  mutate(garden_combo = ifelse(!complete.cases(garden_preg, garden_4, garden_7), NA, 
                                   ifelse(garden_preg == 0 & garden_4 == 0 & garden_7 == 0, 1,
                                   ifelse(garden_preg == 1 & garden_4 == 0 & garden_7 == 0, 2, 
                                   ifelse(garden_preg == 0 & garden_4 == 1 & garden_7 == 0, 3, 
                                   ifelse(garden_preg == 0 & garden_4 == 0 & garden_7 == 1, 4, 
                                   ifelse(garden_preg == 1 & garden_4 == 1 & garden_7 == 0, 5, 
                                   ifelse(garden_preg == 1 & garden_4 == 0 & garden_7 == 1, 6, 
                                   ifelse(garden_preg == 0 & garden_4 == 1 & garden_7 == 1, 7, 8)))))))))

table(data$garden_combo, useNA = "ifany")
prop.table(table(data$garden_combo)) * 100

table(data$garden_combo[!is.na(data$BMI_f7)])

## Turn this date into a sankey plot - To do this, first convert the data to summary wide format, then convert to long/lode format.
data_temp_garden <- data %>%
  dplyr::select(garden_preg, garden_4, garden_7) %>%
  filter(complete.cases(garden_preg, garden_4, garden_7)) %>%
  mutate(garden_preg = as.factor(garden_preg)) %>%
  mutate(garden_preg = recode(garden_preg, "0" = "No", "1" = "Yes")) %>%
  mutate(garden_preg = factor(garden_preg, levels = c("Yes", "No"))) %>%
  mutate(garden_4 = as.factor(garden_4)) %>%
  mutate(garden_4 = recode(garden_4, "0" = "No", "1" = "Yes")) %>%
  mutate(garden_4 = factor(garden_4, levels = c("Yes", "No"))) %>%
  mutate(garden_7 = as.factor(garden_7)) %>%
  mutate(garden_7 = recode(garden_7, "0" = "No", "1" = "Yes")) %>%
  mutate(garden_7 = factor(garden_7, levels = c("Yes", "No"))) %>%
  group_by(garden_preg, garden_4, garden_7) %>%
  summarise(freq = n())

summary(data_temp_garden)
head(data_temp_garden)

data_temp_lodes_garden <- to_lodes_form(data_temp_garden, axes = 1:3, id = "traj")
data_temp_lodes_garden <- data_temp_lodes_garden %>%
  rename(time = x, Response = stratum)
head(data_temp_lodes_garden)
summary(data_temp_lodes_garden)

sankey_garden <- ggplot(data_temp_lodes_garden,
                       aes(x = time, stratum = Response, alluvium = traj,
                           y = freq,
                           fill = Response, label = Response)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("Time-point") + ylab("Frequency")

sankey_garden

# Save this plot
#pdf("sankey_garden.pdf", height = 6, width = 10)
#plot(sankey_garden)
#dev.off()


# All possible combinations have some cell counts, but some are very small (e.g., only 6 people have the pattern 0-1-0) - This drops even more (to 3) if restrict to children with F@7 BMI data... This *might* just about be enough to work with, but will also use just 2 time-points (pregnancy and age 7) to ensure sufficient sample sizes and see if get same pattern of results
data <- data %>%
  mutate(garden_combo2 = ifelse(!complete.cases(garden_preg, garden_7), NA,
                                    ifelse(garden_preg == 0 & garden_7 == 0, 1,
                                    ifelse(garden_preg == 1 & garden_7 == 0, 2,
                                    ifelse(garden_preg == 0 & garden_7 == 1, 3, 4)))))

table(data$garden_combo2, useNA = "ifany")
prop.table(table(data$garden_combo2)) * 100

table(data$garden_combo2[!is.na(data$BMI_f7)])

# Sample sizes slightly better now, but really not by much (min 8 in F@7 sample), and as only 2 time-points can't explore as many life-course trajectories



### SEP interaction covariates

## Parental education - Will choose higher educational qualification of either parent. First, though, check how similar the maternal education variable is from LifeCycle against the ALSPAC var (are identical, so can just drop the LifeCycle one)
table(data$c645a, useNA = "ifany")
table(data$edu_m_0, useNA = "ifany")
table(data$c645a, data$edu_m_0, useNA = "ifany")

data <- data %>%
  dplyr::select(-edu_m_0)

# Combine highest parental education, then code into binary variable (A level or higher vs O level or lower)
table(data$c645a, useNA = "ifany")
table(data$c666a, useNA = "ifany")

data <- data %>%
  mutate(c645a = ifelse(c645a < 0, NA, c645a)) %>%
  mutate(c666a = ifelse(c666a < 0, NA, c666a)) %>%
  mutate(edu = ifelse(is.na(c645a) & is.na(c666a), NA,
                      ifelse(c645a > 3 | c666a > 3, 1, 0)))

table(data$edu, useNA = "ifany")
prop.table(table(data$edu)) * 100
prop.table(table(data$edu, useNA = "ifany")) * 100


## Deprivation/SES - Compare ALSPAC IMD and LifeCycle SES first - Are highly correlated, but not identical. Will just use the LifeCycle data here.
table(data$jan1993imd2010q5_M, useNA = "ifany")
table(data$areaSES_preg, useNA = "ifany")
table(data$jan1993imd2010q5_M, data$areaSES_preg, useNA = "ifany")

data <- data %>%
  dplyr::select(-jan1993imd2010q5_M)

# Now recode into binary deprivation/SES variable, with quintiles 4 and 5 defined as 'deprived'
data <- data %>%
  mutate(deprived = ifelse(is.na(areaSES_preg), NA,
                           ifelse(areaSES_preg > 3, 1, 0)))

table(data$deprived, useNA = "ifany")
prop.table(table(data$deprived)) * 100
prop.table(table(data$deprived, useNA = "ifany")) * 100


## Log equivalised household income - Will code quintiles 1 and 2 as 'low income'
summary(data$eusilc_income)
table(data$eusilc_income_quintiles, useNA = "ifany")

data <- data %>%
  mutate(lowIncome = ifelse(is.na(eusilc_income_quintiles), NA,
                            ifelse(eusilc_income_quintiles < 3, 1, 0)))

table(data$lowIncome, useNA = "ifany")
prop.table(table(data$lowIncome)) * 100
prop.table(table(data$lowIncome, useNA = "ifany")) * 100


### And finally some of the additional confounders/covariates

## Sex of child
table(data$kz021, useNA = "ifany")

data <- data %>%
  rename(male = kz021) %>%
  mutate(male = ifelse(male == 1, 1, 0))

table(data$male, useNA = "ifany")
prop.table(table(data$male, useNA = "ifany")) * 100


## Age of child at F@7 clinic (months)
summary(data$f7003c)

data <- data %>%
  rename(age_f7 = f7003c)

summary(data$age_f7)
sd(data$age_f7, na.rm = TRUE)
hist(data$age_f7)


## Maternal ethnicity - Code into white vs other than white. First, compare ALSPAC and LifeCycle ethnicity data to make sure are similar (are exactly the same, so will just use the ALSPAC data)
table(data$c800, useNA = "ifany")
table(data$ethn3_m, useNA = "ifany")

data <- data %>%
  dplyr::select(-ethn3_m) %>%
  rename(white = c800) %>%
  mutate(white = ifelse(white < 0, NA, white)) %>%
  mutate(white = ifelse(is.na(white), NA,
                            ifelse(white > 1, 0, white)))

table(data$white, useNA = "ifany")
prop.table(table(data$white)) * 100
prop.table(table(data$white, useNA = "ifany")) * 100



#####################################################################################################
##### Access to green space models

### Derive the different life-course variables and encode as hypotheses
table(data$greenSpace_preg, useNA = "ifany")
table(data$greenSpace_4, useNA = "ifany")
table(data$greenSpace_7, useNA = "ifany")
table(data$edu, useNA = "ifany")

# Make a dataset just for this exposure
data_access <- data


## Start with education as the SEP covariate/interaction term (to try and reduce the correlation between the hypotheses and the interaction terms and improve power, will center all these hypotheses)

# Critical period at first time point only
data_access$crit1 <- data_access$greenSpace_preg
data_access$crit1 <- data_access$crit1 - mean(data_access$crit1, na.rm = TRUE)

# Interaction between SEP and first time point
data_access$int1 <- data_access$edu * data_access$crit1

# Critical period at second time point only
data_access$crit2 <- data_access$greenSpace_4
data_access$crit2 <- data_access$crit2 - mean(data_access$crit2, na.rm = TRUE)

# Interaction between SEP and second time point
data_access$int2 <- data_access$edu * data_access$crit2

# Critical period at third time point only
data_access$crit3 <- data_access$greenSpace_7
data_access$crit3 <- data_access$crit3 - mean(data_access$crit3, na.rm = TRUE)

# Interaction between SEP and third time point
data_access$int3 <- data_access$edu * data_access$crit3

# Linear accumulation of all exposures
data_access$accumulation <- data_access$greenSpace_preg + data_access$greenSpace_4 + data_access$greenSpace_7
data_access$accumulation <- data_access$accumulation - mean(data_access$accumulation, na.rm = TRUE)

# Interaction between SEP and cumulative exposure
data_access$int_accum <- data_access$edu * data_access$accumulation

# Increase in access to green space from time 1 to time 2
data_access$green_inc12 <- (1 - data_access$greenSpace_preg) * data_access$greenSpace_4
data_access$green_inc12 <- data_access$green_inc12 - mean(data_access$green_inc12, na.rm = TRUE)

# Increase in access to green space from time 1 to time 2, with an interaction with SEP
data_access$green_inc12_int <- data_access$green_inc12 * data_access$edu

# Decrease in access to green space from time 1 to time 2
data_access$green_dec12 <- (1 - data_access$greenSpace_4) * data_access$greenSpace_preg
data_access$green_dec12 <- data_access$green_dec12 - mean(data_access$green_dec12, na.rm = TRUE)

# Decrease in access to green space from time 1 to time 2, with an interaction with SEP
data_access$green_dec12_int <- data_access$green_dec12 * data_access$edu

# Increase in access to green space from time 2 to time 3
data_access$green_inc23 <- (1 - data_access$greenSpace_4) * data_access$greenSpace_7
data_access$green_inc23 <- data_access$green_inc23 - mean(data_access$green_inc23, na.rm = TRUE)

# Increase in access to green space from time 2 to time 3, with an interaction with SEP
data_access$green_inc23_int <- data_access$green_inc23 * data_access$edu

# Decrease in access to green space from time 2 to time 3
data_access$green_dec23 <- (1 - data_access$greenSpace_7) * data_access$greenSpace_4
data_access$green_dec23 <- data_access$green_dec23 - mean(data_access$green_dec23, na.rm = TRUE)

# Decrease in access to green space from time 2 to time 3, with an interaction with SEP
data_access$green_dec23_int <- data_access$green_dec23 * data_access$edu


## Make a dataset just with the outcomes, exposure hypotheses and covariates
data_access_edu <- data_access %>%
  dplyr::select(BMI_f7, overweight, sysBP, diaBP, 
         age_f7, male, white, edu, 
         crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_inc12, green_inc12_int, 
         green_dec12, green_dec12_int, green_inc23, green_inc23_int, green_dec23, green_dec23_int)


## Now analyse each outcome in turn, starting with BMI

# Reduce dataset down to just complete cases
data_access_edu_bmi <- data_access_edu %>%
  dplyr::select(-overweight, -sysBP, -diaBP) %>%
  filter(complete.cases(BMI_f7, age_f7, male, white, edu, crit1, int1))

summary(data_access_edu_bmi)
nrow(data_access_edu_bmi)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_edu_bmi %>%
  dplyr::select(-BMI_f7)

x_hypos <- as.matrix(x_hypos)


## Check correlation matrix of all these hypotheses (>0.9 would be a cause for concern)
dim(x_hypos)
cor(x_hypos[,5:20])
cor(x_hypos[,5:20]) > 0.9
cor(x_hypos[,5:20]) > 0.95

# Biggest issues is with the accumulation variable, which is highly correlated with the critical period variables, so will drop this accumulation variable (and its SEP-interaction term) from these analysis due to collinearity and effectively measuring the same thing. After centering, the critical period interaction terms were no longer highly-correlated (r<0.9).
x_hypos <- x_hypos[,!colnames(x_hypos) %in% c("accumulation", "int_accum")]
head(x_hypos)

dim(x_hypos)
cor(x_hypos[,5:18])
cor(x_hypos[,5:18]) > 0.9
cor(x_hypos[,5:18]) > 0.95


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_bmi <- glmnet(x_hypos, data_access_edu_bmi$BMI_f7, 
              alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_edu_bmi


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, use the 'lasso-table' function defined above to pick out the changes in variables and the increment in deviance ratio
df <- lasso_table(mod_access_edu_bmi)
df

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_edu_bmi$log_lambda <- log(mod_access_edu_bmi$lambda)
df$log_lambda <- log(as.numeric(df$Lambda))
df

# Make the plot
plot(mod_access_edu_bmi$log_lambda, mod_access_edu_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_bmi$log_lambda)), ylim = c(0.007, max(mod_access_edu_bmi$dev.ratio)))
text(df$log_lambda, 0.007, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessEduBMI.pdf", height = 6, width = 10)
plot(mod_access_edu_bmi$log_lambda, mod_access_edu_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_bmi$log_lambda)), ylim = c(0.007, max(mod_access_edu_bmi$dev.ratio)))
text(df$log_lambda, 0.007, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will run all combinations of variables in a standard LM/GLM model, and store the AIC and BIC values
for (i in 1:nrow(df)) {
  
  vars_temp <- strsplit(df$model_vars[i], " ")[[1]] # Split the variables at each stage of the lasso
  x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp] # Matrix of the covariates at each step of lasso
  
  # Run the model
  mod_ic <- lm(data_access_edu_bmi$BMI_f7 ~ x_hypos_new)
  print(summary(mod_ic))

  # Store the AIC values
  df$aic[i] <- AIC(mod_ic)
  df$bic[i] <- BIC(mod_ic)
  
}

# Select the models with the lowest AIC and BIC values - Both models selected the baseline covariate-only model
df

df$model_vars[which.min(df$aic)]
df$model_vars[which.min(df$bic)]



## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' (minimum MSE) and '1 SE' models, although the 1 SE model is may be better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model
mod.cv <- cv.glmnet(x_hypos, data_access_edu_bmi$BMI_f7, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters runs along the top of the plot)
plot(mod.cv) 

# The 1SE model contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)

# The optimal' model gives the same result (and as this optimal lasso is mainly for prediction, one may expect an increase in model complexity, suggesting again that green space effects are essentially non-existent)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between access to green space in childhood and BMI at age 7 (when controlling for ethnicity and parental education).



#### As a sensitivity analysis, will check whether get same results using just 2 time-points (pregnancy and age 7).

data_access2 <- data

## Encode the hypotheses (and center variables)

# Critical period at first time point only
data_access2$crit1 <- data_access2$greenSpace_preg
data_access2$crit1 <- data_access2$crit1 - mean(data_access2$crit1, na.rm = TRUE)

# Interaction between SEP and first time point
data_access2$int1 <- data_access2$edu * data_access2$crit1

# Critical period at second time point only
data_access2$crit2 <- data_access2$greenSpace_7
data_access2$crit2 <- data_access2$crit2 - mean(data_access2$crit2, na.rm = TRUE)

# Interaction between SEP and second time point
data_access2$int2 <- data_access2$edu * data_access2$crit2

# Linear accumulation of all exposures
data_access2$accumulation <- data_access2$greenSpace_preg + data_access2$greenSpace_7
data_access2$accumulation <- data_access2$accumulation - mean(data_access2$accumulation, na.rm = TRUE)

# Interaction between SEP and cumulative exposure
data_access2$int_accum <- data_access2$edu * data_access2$accumulation

# Increase in access to green space from time 1 to time 2
data_access2$green_inc12 <- (1 - data_access2$greenSpace_preg) * data_access2$greenSpace_7
data_access2$green_inc12 <- data_access2$green_inc12 - mean(data_access2$green_inc12, na.rm = TRUE)

# Increase in access to green space from time 1 to time 2, with an interaction with SEP
data_access2$green_inc12_int <- data_access2$green_inc12 * data_access2$edu

# Decrease in access to green space from time 1 to time 2
data_access2$green_dec12 <- (1 - data_access2$greenSpace_7) * data_access2$greenSpace_preg
data_access2$green_dec12 <- data_access2$green_dec12 - mean(data_access2$green_dec12, na.rm = TRUE)

# Decrease in access to green space from time 1 to time 2, with an interaction with SEP
data_access2$green_dec12_int <- data_access2$green_dec12 * data_access2$edu


## Make a dataset just with the outcomes, exposure hypotheses and covariates
data_access2_edu <- data_access2 %>%
  dplyr::select(BMI_f7, BMI_f7_z, sysBP, diaBP, 
         age_f7, male, white, edu, 
         crit1, int1, crit2, int2, accumulation, int_accum, green_inc12, green_inc12_int, 
         green_dec12, green_dec12_int)


## Now analyse each outcome in turn, starting with BMI

# Reduce dataset down to just complete cases
data_access2_edu_bmi <- data_access2_edu %>%
  dplyr::select(-BMI_f7_z, -sysBP, -diaBP) %>%
  filter(complete.cases(BMI_f7, age_f7, male, white, edu, crit1, int1))

summary(data_access2_edu_bmi)

# Save the life-course hypotheses and covariates as a matrix
x_hypos2 <- data_access2_edu_bmi %>%
  dplyr::select(-BMI_f7)

x_hypos2 <- as.matrix(x_hypos2)


## Check correlation matrix of all these hypotheses (>0.9 would be a cause for concern)
dim(x_hypos2)
cor(x_hypos2[,5:14])
cor(x_hypos2[,5:14]) > 0.9
cor(x_hypos2[,5:14]) > 0.95

# Biggest issues are with the accumulation variables, which are highly correlated with the critical period variables, so will drop these accumulation variables from these analysis due to collinearity and effectively measuring the same thing.
x_hypos2 <- x_hypos2[,!colnames(x_hypos2) %in% c("accumulation", "int_accum")]
head(x_hypos2)

dim(x_hypos2)
cor(x_hypos2[,5:12])
cor(x_hypos2[,5:12]) > 0.9
cor(x_hypos2[,5:12]) > 0.95


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_bmi2 <- glmnet(x_hypos2, data_access2_edu_bmi$BMI_f7, 
                             alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos2) - 4))))

mod_access_edu_bmi2


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, use the 'lasso-table' function defined above to pick out the changes in variables and the increment in deviance ratio
df <- lasso_table(mod_access_edu_bmi2)
df


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_edu_bmi2$log_lambda <- log(mod_access_edu_bmi2$lambda)
df$log_lambda <- log(as.numeric(df$Lambda))
df

# Make the plot
plot(mod_access_edu_bmi2$log_lambda, mod_access_edu_bmi2$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_bmi2$log_lambda)), ylim = c(0.008, max(mod_access_edu_bmi2$dev.ratio)))
text(df$log_lambda, 0.008, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessEduBMI_2Exposures.pdf", height = 6, width = 10)
plot(mod_access_edu_bmi2$log_lambda, mod_access_edu_bmi2$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_bmi2$log_lambda)), ylim = c(0.008, max(mod_access_edu_bmi2$dev.ratio)))
text(df$log_lambda, 0.008, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## Will run all combinations of variables in a standard LM/GLM model, and store the AIC and BIC values
for (i in 1:nrow(df)) {
  
  vars_temp <- strsplit(df$model_vars[i], " ")[[1]] # Split the variables at each stage of the lasso
  x_hypos_new <- x_hypos2[, colnames(x_hypos2) %in% vars_temp] # Matrix of the covariates at each step of lasso
  
  # Run the model
  mod_ic <- lm(data_access2_edu_bmi$BMI_f7 ~ x_hypos_new)
  print(summary(mod_ic))
  
  # Store the AIC values
  df$aic[i] <- AIC(mod_ic)
  df$bic[i] <- BIC(mod_ic)
  
}

# Select the models with the lowest AIC and BIC values - Both models selected the baseline covariate-only model
df

df$model_vars[which.min(df$aic)]
df$model_vars[which.min(df$bic)]


# Alternative method using cross-validated lasso (again, both give null/covariate only model as best fit)
mod.cv2 <- cv.glmnet(x_hypos2, data_access2_edu_bmi$BMI_f7, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos2) - 4))))
mod.cv2

coef(mod.cv2, s = mod.cv2$lambda.1se)
coef(mod.cv2, s = mod.cv2$lambda.min)



#### Next, will explore binary access to green space exposure, with education as covariate/interaction term, and binary overweight variable as the outcome

# Reduce dataset down to just complete cases
data_access_edu_overweight <- data_access_edu %>%
  dplyr::select(-BMI_f7, -sysBP, -diaBP) %>%
  filter(complete.cases(overweight, age_f7, male, white, edu, crit1, int1))

summary(data_access_edu_overweight)
nrow(data_access_edu_overweight)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_edu_overweight %>%
  dplyr::select(-overweight)

x_hypos <- as.matrix(x_hypos)


## Remove the accumulation variables (as high collinearity with critical period hypotheses)
x_hypos <- x_hypos[,!colnames(x_hypos) %in% c("accumulation", "int_accum")]
head(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_over <- glmnet(x_hypos, data_access_edu_overweight$overweight, family = "binomial",
                             alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_edu_over

# Plot these results
plot(mod_access_edu_over)


### Visual inspection of results 

# First, use the 'lasso-table' function defined above to pick out the changes in variables and the increment in deviance ratio
df <- lasso_table(mod_access_edu_over)
df

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_edu_over$log_lambda <- log(mod_access_edu_over$lambda)
mod_access_edu_over

df$log_lambda <- log(as.numeric(df$Lambda))
df

# Make plot
plot(mod_access_edu_over$log_lambda, mod_access_edu_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_over$log_lambda)), ylim = c(0.001, max(mod_access_edu_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessEduOverweight.pdf", height = 6, width = 10)
plot(mod_access_edu_over$log_lambda, mod_access_edu_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_over$log_lambda)), ylim = c(0.001, max(mod_access_edu_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will run all combinations of variables in a standard LM/GLM model, and store the AIC and BIC values
for (i in 1:nrow(df)) {
  
  vars_temp <- strsplit(df$model_vars[i], " ")[[1]] # Split the variables at each stage of the lasso
  x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp] # Matrix of the covariates at each step of lasso
  
  # Run the model
  mod_ic <- glm(data_access_edu_overweight$overweight ~ x_hypos_new, family = "binomial")
  print(summary(mod_ic))
  
  # Store the AIC values
  df$aic[i] <- AIC(mod_ic)
  df$bic[i] <- BIC(mod_ic)
  
}

# Select the models with the lowest AIC and BIC values - Both models selected the baseline covariate-only model
df

df$model_vars[which.min(df$aic)]
df$model_vars[which.min(df$bic)]


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_edu_overweight$overweight, family = "binomial", 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between access to green space in childhood and being overweight at age 7 (when controlling for ethnicity and parental education).



#### Next, will explore binary access to green space exposure, with education as covariate/interaction term, and continuous systolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_access_edu_sysBP <- data_access_edu %>%
  dplyr::select(-BMI_f7, -overweight, -diaBP) %>%
  filter(complete.cases(sysBP, age_f7, male, white, edu, crit1, int1))

summary(data_access_edu_sysBP)
nrow(data_access_edu_sysBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_edu_sysBP %>%
  dplyr::select(-sysBP)

x_hypos <- as.matrix(x_hypos)

## Remove the accumulation variables (as high collinearity with critical period hypotheses)
x_hypos <- x_hypos[,!colnames(x_hypos) %in% c("accumulation", "int_accum")]
head(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_sysBP <- glmnet(x_hypos, data_access_edu_sysBP$sysBP,
                              alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_edu_sysBP

# Plot these results
plot(mod_access_edu_sysBP)

### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, use the 'lasso-table' function defined above to pick out the changes in variables and the increment in deviance ratio
df <- lasso_table(mod_access_edu_sysBP)
df

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_edu_sysBP$log_lambda <- log(mod_access_edu_sysBP$lambda)
mod_access_edu_sysBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# Make the plot
plot(mod_access_edu_sysBP$log_lambda, mod_access_edu_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_sysBP$log_lambda)), ylim = c(0.005, max(mod_access_edu_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessEduSysBP.pdf", height = 6, width = 10)
plot(mod_access_edu_sysBP$log_lambda, mod_access_edu_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_sysBP$log_lambda)), ylim = c(0.005, max(mod_access_edu_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will run all combinations of variables in a standard LM/GLM model, and store the AIC and BIC values
for (i in 1:nrow(df)) {
  
  vars_temp <- strsplit(df$model_vars[i], " ")[[1]] # Split the variables at each stage of the lasso
  x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp] # Matrix of the covariates at each step of lasso
  
  # Run the model
  mod_ic <- lm(data_access_edu_sysBP$sysBP ~ x_hypos_new)
  print(summary(mod_ic))
  
  # Store the AIC values
  df$aic[i] <- AIC(mod_ic)
  df$bic[i] <- BIC(mod_ic)
  
}

# Select the models with the lowest AIC and BIC values - BIC selected the baseline covariate-only model, while AIC selected a much more complicated model with 5 additional parameters
df

df$model_vars[which.min(df$aic)]
df$model_vars[which.min(df$bic)]


## Run selective inference on the best-fitting AIC lasso model

# select the chosen lambda value
lambda <- as.numeric(df$Lambda[which.min(df$aic)])

# Extract the coefficients (have to exclude the intercept)
beta <- coef(mod_access_edu_sysBP, s = lambda, x = x_hypos, y = data_access_edu_sysBP$sysBP, 
             penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))), exact = TRUE)[-1]

# Perform selective inference
si_res <- fixedLassoInf(x_hypos, data_access_edu_sysBP$sysBP, beta, lambda, alpha = 0.05)
si_res
si_res$vars

# Save results in nicer data frame
si_res <- as.data.frame(cbind(si_res$vars, si_res$coef0, si_res$sd, si_res$pv, si_res$ci))
colnames(si_res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res[,-1] <- lapply(si_res[,-1], function(x) round(as.numeric(as.character(x)), 3))
si_res

## And compare against standard regression results (all very similar, even p-values and confidence intervals)
vars_temp <- strsplit(df$model_vars[which.min(df$aic)], " ")[[1]]
x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp]
head(x_hypos_new)

mod <- lm(data_access_edu_sysBP$sysBP ~ x_hypos_new)
summary(mod)
confint(mod)


## Looking at fitted values to try and understand/interpret results
temp <- as.data.frame(x_hypos_new)
temp$fit <- predict(mod)

# Lowest AIC model suggests an interaction between SEP and decline in green space from time 1 to time 2, such that a decline in green space is associated with lower systolic BP in those from lower SEP backgrounds, while has little association with higher SEP backgrounds. Odd and not particularly realistic effects (why would reducing the amount of green space lower BP, and only in those from lower SEP backgrounds?), so could just be noise.
summary(temp$fit[temp$edu == 0 & temp$green_dec12 == min(temp$green_dec12)]) # Low SEP and no reduction in green space
summary(temp$fit[temp$edu == 0 & temp$green_dec12 == max(temp$green_dec12)]) # Low SEP and reduction in green space
summary(temp$fit[temp$edu == 1 & temp$green_dec12 == min(temp$green_dec12)]) # High SEP and no reduction in green space
summary(temp$fit[temp$edu == 1 & temp$green_dec12 == max(temp$green_dec12)]) # High SEP and reduction in green space


## We can also force the main effects to be included in the reported model, if the interaction term was detected by the lasso - This applies to 'int3' (will include main effect of 'crit3') and 'green_inc23_int' (will include main effect of 'green_inc23'). Note: Have had to force all of the encoded variables from the original model here, in addition to the main effects, as if just force the main effects in then the model changes and some previously-included variables are dropped - Compare 'beta, 'beta_test' and 'beta_mainEffects'.
beta_test <- coef(mod_access_edu_sysBP, s = lambda, x = x_hypos, y = data_access_edu_sysBP$sysBP, 
                         penalty.factor = (c(0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1)), exact = TRUE)[-1]

beta_mainEffects <- coef(mod_access_edu_sysBP, s = lambda, x = x_hypos, y = data_access_edu_sysBP$sysBP, 
             penalty.factor = (c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1)), exact = TRUE)[-1]

beta; beta_test; beta_mainEffects

# Perform selective inference
si_res_mainEffects <- fixedLassoInf(x_hypos, data_access_edu_sysBP$sysBP, beta_mainEffects, lambda, alpha = 0.05)
si_res_mainEffects
si_res_mainEffects$vars

# Save results in nicer data frame
si_res_mainEffects <- as.data.frame(cbind(si_res_mainEffects$vars, si_res_mainEffects$coef0, 
                                          si_res_mainEffects$sd, si_res_mainEffects$pv, si_res_mainEffects$ci))
colnames(si_res_mainEffects) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res_mainEffects[,-1] <- lapply(si_res_mainEffects[,-1], function(x) round(as.numeric(as.character(x)), 3))
si_res_mainEffects



## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_edu_sysBP$sysBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - The 1SE just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure, although the minimum mean squared error value has 15 parameters (most of the hypotheses); however, the mean squared error is practically identical, suggesting little improvement in model fit.
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


## Run selective inference on the best-fitting CV min MSE lasso model

# select the chosen lambda value (as noted in the script for the example simulation of these methods, I've found that the 'penalty.factor' argument doesn't always work if the lambda value is not manually specified - I have no idea why...)
(lambda <- mod.cv$lambda.min)
lambda <- 0.0057156

# Extract the coefficients (have to exclude the intercept)
beta <- coef(mod.cv, s = lambda, x = x_hypos, y = data_access_edu_sysBP$sysBP, 
             penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))), exact = TRUE)[-1]

# Perform selective inference
si_res <- fixedLassoInf(x_hypos, data_access_edu_sysBP$sysBP, beta, lambda, alpha = 0.05)
si_res
si_res$vars

# Save results in nicer data frame
si_res <- as.data.frame(cbind(si_res$vars, si_res$coef0, si_res$sd, si_res$pv, si_res$ci))
colnames(si_res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res[,-1] <- lapply(si_res[,-1], function(x) round(as.numeric(as.character(x)), 3))
si_res


## We can also force the main effects to be included in the reported model, if the interaction term was detected by the lasso - This applies to 'int2' (will include main effect of 'crit2'). Note: Have had to force all of the encoded variables from the original model here, in addition to the main effects, as if just force the main effects in then the model changes and some previously-included variables are dropped - Compare 'beta, 'beta_test' and 'beta_mainEffects'.
beta_test <- coef(mod.cv, s = lambda, x = x_hypos, y = data_access_edu_sysBP$sysBP, 
                  penalty.factor = (c(0, 0, 0, 0, 1, 1, 0, rep(1, 11))), exact = TRUE)[-1]

beta_mainEffects <- coef(mod.cv, s = lambda, x = x_hypos, y = data_access_edu_sysBP$sysBP, 
                         penalty.factor = (c(0, 0, 0, 0, 1, 1, rep(0, 12))), exact = TRUE)[-1]

beta; beta_test; beta_mainEffects

# Perform selective inference (although here these parameters could not be calculated)
si_res_mainEffects <- fixedLassoInf(x_hypos, data_access_edu_sysBP$sysBP, beta_mainEffects, lambda, alpha = 0.05)
si_res_mainEffects
si_res_mainEffects$vars

# Save results in nicer data frame
si_res_mainEffects <- as.data.frame(cbind(si_res_mainEffects$vars, si_res_mainEffects$coef0, 
                                          si_res_mainEffects$sd, si_res_mainEffects$pv, si_res_mainEffects$ci))
colnames(si_res_mainEffects) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res_mainEffects[,-1] <- lapply(si_res_mainEffects[,-1], function(x) round(as.numeric(as.character(x)), 3))
si_res_mainEffects


### Is some variation between the methods, but suggests no strong association between access to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and parental education).



#### Next, will explore binary access to green space exposure, with education as covariate/interaction term, and continuous diastolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_access_edu_diaBP <- data_access_edu %>%
  dplyr::select(-BMI_f7, -overweight, -sysBP) %>%
  filter(complete.cases(diaBP, age_f7, male, white, edu, crit1, int1))

summary(data_access_edu_diaBP)
nrow(data_access_edu_diaBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_edu_diaBP %>%
  dplyr::select(-diaBP)

x_hypos <- as.matrix(x_hypos)

## Remove the accumulation variables (as high collinearity with critical period hypotheses)
x_hypos <- x_hypos[,!colnames(x_hypos) %in% c("accumulation", "int_accum")]
head(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_diaBP <- glmnet(x_hypos, data_access_edu_diaBP$diaBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_edu_diaBP

# Plot these results
plot(mod_access_edu_diaBP)



### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, use the 'lasso-table' function defined above to pick out the changes in variables and the increment in deviance ratio
df <- lasso_table(mod_access_edu_diaBP)
df

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_edu_diaBP$log_lambda <- log(mod_access_edu_diaBP$lambda)
mod_access_edu_diaBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# Make the plot
plot(mod_access_edu_diaBP$log_lambda, mod_access_edu_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_diaBP$log_lambda)), ylim = c(0.004, max(mod_access_edu_diaBP$dev.ratio)))
text(df$log_lambda, 0.004, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessEduDiaBP.pdf", height = 6, width = 10)
plot(mod_access_edu_diaBP$log_lambda, mod_access_edu_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_diaBP$log_lambda)), ylim = c(0.004, max(mod_access_edu_diaBP$dev.ratio)))
text(df$log_lambda, 0.004, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that little is really going on here (although perhaps weak evidence that 'green_dec12' associated with outcome) - Will run all combinations of variables in a standard LM/GLM model, and store the AIC and BIC values
for (i in 1:nrow(df)) {
  
  vars_temp <- strsplit(df$model_vars[i], " ")[[1]] # Split the variables at each stage of the lasso
  x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp] # Matrix of the covariates at each step of lasso
  
  # Run the model
  mod_ic <- lm(data_access_edu_diaBP$diaBP ~ x_hypos_new)
  print(summary(mod_ic))
  
  # Store the AIC values
  df$aic[i] <- AIC(mod_ic)
  df$bic[i] <- BIC(mod_ic)
  
}

# Select the models with the lowest AIC and BIC values - BIC selected the baseline covariate-only model, while AIC selected the additional variable 'green_dec_12'
df

df$model_vars[which.min(df$aic)]
df$model_vars[which.min(df$bic)]

## Run selective inference on the best-fitting AIC lasso model

# select the chosen lambda value
lambda <- as.numeric(df$Lambda[which.min(df$aic)])

# Extract the coefficients (have to exclude the intercept)
beta <- coef(mod_access_edu_diaBP, s = lambda, x = x_hypos, y = data_access_edu_diaBP$diaBP, 
             penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))), exact = TRUE)[-1]

# Perform selective inference
si_res <- fixedLassoInf(x_hypos, data_access_edu_diaBP$diaBP, beta, lambda, alpha = 0.05)
si_res
si_res$vars

# Save results in nicer data frame
si_res <- as.data.frame(cbind(si_res$vars, si_res$coef0, si_res$sd, si_res$pv, si_res$ci))
colnames(si_res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res[,-1] <- lapply(si_res[,-1], function(x) round(as.numeric(as.character(x)), 3))
si_res

## And compare against standard regression results (all very similar, even p-values and confidence intervals)
vars_temp <- strsplit(df$model_vars[which.min(df$aic)], " ")[[1]]
x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp]
head(x_hypos_new)

mod <- lm(data_access_edu_diaBP$diaBP ~ x_hypos_new)
summary(mod)

# Lowest AIC model suggests that decline in green space from time 1 to time 2 is associated with lower diastolic BP. Odd and not particularly realistic effects (why would reducing the amount of green space lower BP?), so could just be noise.


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_edu_diaBP$diaBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - The 1SE model is just the null/covariate-only model, while the 'optimal'/lowest MSE model also contains green_dec12. The difference in MSE is minimal, though, so not especially convincing. Suggests at best a very weak association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


## Run selective inference on the best-fitting CV min MSE lasso model

# select the chosen lambda value
(lambda <- mod.cv$lambda.min)
lambda <- 0.1139

# Extract the coefficients (have to exclude the intercept)
beta <- coef(mod.cv, s = lambda, x = x_hypos, y = data_access_edu_diaBP$diaBP, 
             penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))), exact = TRUE)[-1]

# Perform selective inference
si_res <- fixedLassoInf(x_hypos, data_access_edu_diaBP$diaBP, beta, lambda, alpha = 0.05)
si_res
si_res$vars

# Save results in nicer data frame
si_res <- as.data.frame(cbind(si_res$vars, si_res$coef0, si_res$sd, si_res$pv, si_res$ci))
colnames(si_res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res[,-1] <- lapply(si_res[,-1], function(x) round(as.numeric(as.character(x)), 3))
si_res


### All methods match up relatively well, and indicate little association between access to green space in childhood and diastolic blood pressure at age 7 (when controlling for ethnicity and parental education), although potentially decrease in access to green space between preg and age 4 weakly associated with lower diastolic BP.



