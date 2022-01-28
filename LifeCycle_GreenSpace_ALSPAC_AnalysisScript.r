### LifeCycle Green Space Analysis - ALSPAC analysis script

## Created 27/1/2021 by Dan Major-Smith
## R version 4.0.4


## Make sure workspace is clear
rm(list = ls())

# Install and load packages
#install.packages("tidyverse")
library(tidyverse)

#install.packages("glmnet")
library(glmnet)

#install.packages("dagitty")
library(dagitty)


# Set working directory
setwd("C:\\Users\\ds16565\\OneDrive - University of Bristol\\MyFiles-Migrated\\Documents\\Projects\\Lifecycle\\Results")


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



####################################################################################################
##### Read in data
data_raw <- read_csv(file = "../Data/B3930_lifeCycleAndALSPAC.csv")

## Keep just relevant variables
data <- data_raw %>%
  select(aln, qlet, 
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


## Systolic and diastolic blood pressure
summary(data$f7sa021); summary(data$f7sa022)

data <- data %>%
  rename(sysBP = f7sa021, diaBP = f7sa022) %>%
  mutate(sysBP = ifelse(sysBP < 0, NA, sysBP)) %>%
  mutate(diaBP = ifelse(diaBP < 0, NA, diaBP))

summary(data$sysBP); summary(data$diaBP)
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
  select(-cc)

table(data$greenSpace_preg, useNA = "ifany")
table(data$greenSpace_4, useNA = "ifany")
table(data$greenSpace_7, useNA = "ifany")

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
  select(-cc)

summary(data$greenDist_preg); summary(data$greenDist_4); summary(data$greenDist_7)
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
  select(-cc)

summary(data$NDVI500_preg); summary(data$NDVI500_4); summary(data$NDVI500_7)
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
  select(-cc)

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
table(data$garden_4, useNA = "ifany")
table(data$garden_7, useNA = "ifany")

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
  select(-edu_m_0)

# Combine highest parental education, then code into binary variable (A level or higher vs O level or lower)
table(data$c645a, useNA = "ifany")
table(data$c666a, useNA = "ifany")

data <- data %>%
  mutate(c645a = ifelse(c645a < 0, NA, c645a)) %>%
  mutate(c666a = ifelse(c666a < 0, NA, c666a)) %>%
  mutate(edu = ifelse(is.na(c645a) & is.na(c666a), NA,
                      ifelse(c645a > 3 | c666a > 3, 1, 0)))

table(data$edu, useNA = "ifany")


## Deprivation/SES - Compare ALSPAC IMD and LifeCycle SES first - Are highly correlated, but not identical. Will just use the LifeCycle data here.
table(data$jan1993imd2010q5_M, useNA = "ifany")
table(data$areaSES_preg, useNA = "ifany")
table(data$jan1993imd2010q5_M, data$areaSES_preg, useNA = "ifany")

data <- data %>%
  select(-jan1993imd2010q5_M)

# Now recode into binary deprivation/SES variable, with quintiles 4 and 5 defined as 'deprived'
data <- data %>%
  mutate(deprived = ifelse(is.na(areaSES_preg), NA,
                           ifelse(areaSES_preg > 3, 1, 0)))

table(data$deprived, useNA = "ifany")


## Log equivalised household income - Will code quintiles 1 and 2 as 'low income'
summary(data$eusilc_income)
table(data$eusilc_income_quintiles, useNA = "ifany")

data <- data %>%
  mutate(lowIncome = ifelse(is.na(eusilc_income_quintiles), NA,
                            ifelse(eusilc_income_quintiles < 3, 1, 0)))

table(data$lowIncome, useNA = "ifany")



### And finally some of the additional confounders/covariates

## Sex of child
table(data$kz021, useNA = "ifany")

data <- data %>%
  rename(male = kz021) %>%
  mutate(male = ifelse(male == 1, 1, 0))

table(data$male, useNA = "ifany")


## Age of child at F@7 clinic (months)
summary(data$f7003c)

data <- data %>%
  rename(age_f7 = f7003c)

hist(data$age_f7)


## Maternal ethnicity - Code into white vs other than white. First, compare ALSPAC and LifeCycle ethnicity data to make sure are similar (are exactly the same, so will just use the ALSPAC data)
table(data$c800, useNA = "ifany")
table(data$ethn3_m, useNA = "ifany")

data <- data %>%
  select(-ethn3_m) %>%
  rename(white = c800) %>%
  mutate(white = ifelse(white < 0, NA, white)) %>%
  mutate(white = ifelse(is.na(white), NA,
                            ifelse(white > 1, 0, white)))

table(data$white, useNA = "ifany")



#####################################################################################################
##### Access to green space models

### Derive the different life-course variables and encode as hypotheses
table(data$greenSpace_preg, useNA = "ifany")
table(data$greenSpace_4, useNA = "ifany")
table(data$greenSpace_7, useNA = "ifany")
table(data$edu, useNA = "ifany")

# Make a dataset just for this exposure
data_access <- data


## Start with education as the SEP covariate/interaction term

# Critical period at first time point only
data_access$crit1 <- data_access$greenSpace_preg

# Interaction between SEP and first time point
data_access$int1 <- data_access$edu * data_access$greenSpace_preg

# Critical period at second time point only
data_access$crit2 <- data_access$greenSpace_4

# Interaction between SEP and second time point
data_access$int2 <- data_access$edu * data_access$greenSpace_4

# Critical period at third time point only
data_access$crit3 <- data_access$greenSpace_7

# Interaction between SEP and third time point
data_access$int3 <- data_access$edu * data_access$greenSpace_7

# Linear accumulation of all exposures
data_access$accumulation <- data_access$greenSpace_preg + data_access$greenSpace_4 + data_access$greenSpace_7

# Interaction between SEP and cumulative exposure
data_access$int_accum <- data_access$edu * data_access$accumulation

# Increase in access to green space from time 1 to time 2
data_access$green_inc12 <- (1 - data_access$greenSpace_preg) * data_access$greenSpace_4

# Increase in access to green space from time 1 to time 2, with an interaction with SEP
data_access$green_inc12_int <- (1 - data_access$greenSpace_preg) * data_access$greenSpace_4 * data_access$edu

# Decrease in access to green space from time 1 to time 2
data_access$green_dec12 <- (1 - data_access$greenSpace_4) * data_access$greenSpace_preg

# Decrease in access to green space from time 1 to time 2, with an interaction with SEP
data_access$green_dec12_int <- (1 - data_access$greenSpace_4) * data_access$greenSpace_preg * data_access$edu

# Increase in access to green space from time 2 to time 3
data_access$green_inc23 <- (1 - data_access$greenSpace_4) * data_access$greenSpace_7

# Increase in access to green space from time 2 to time 3, with an interaction with SEP
data_access$green_inc23_int <- (1 - data_access$greenSpace_4) * data_access$greenSpace_7 * data_access$edu

# Decrease in access to green space from time 2 to time 3
data_access$green_dec23 <- (1 - data_access$greenSpace_7) * data_access$greenSpace_4

# Decrease in access to green space from time 2 to time 3, with an interaction with SEP
data_access$green_dec23_int <- (1 - data_access$greenSpace_7) * data_access$greenSpace_4 * data_access$edu


## Make a dataset just with the outcomes, exposure hypotheses and covariates
data_access_edu <- data_access %>%
  select(BMI_f7, overweight, sysBP, diaBP, 
         age_f7, male, white, edu, 
         crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_inc12, green_inc12_int, 
         green_dec12, green_dec12_int, green_inc23, green_inc23_int, green_dec23, green_dec23_int)


## Now analyse each outcome in turn, starting with BMI

# Reduce dataset down to just complete cases
data_access_edu_bmi <- data_access_edu %>%
  select(-overweight, -sysBP, -diaBP) %>%
  filter(complete.cases(BMI_f7, age_f7, male, white, edu, crit1, int1))

summary(data_access_edu_bmi)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_edu_bmi %>%
  select(-BMI_f7)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_bmi <- glmnet(x_hypos, data_access_edu_bmi$BMI_f7, 
              alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_edu_bmi

# Plot these results
plot(mod_access_edu_bmi)

# Look at the variables included at each step
coef(mod_access_edu_bmi, s = max(mod_access_edu_bmi$lambda[mod_access_edu_bmi$df == 4])); min(mod_access_edu_bmi$dev.ratio[mod_access_edu_bmi$df == 4])

coef(mod_access_edu_bmi, s = max(mod_access_edu_bmi$lambda[mod_access_edu_bmi$df == 5])); min(mod_access_edu_bmi$dev.ratio[mod_access_edu_bmi$df == 5])

coef(mod_access_edu_bmi, s = max(mod_access_edu_bmi$lambda[mod_access_edu_bmi$df == 6])); min(mod_access_edu_bmi$dev.ratio[mod_access_edu_bmi$df == 6])

coef(mod_access_edu_bmi, s = max(mod_access_edu_bmi$lambda[mod_access_edu_bmi$df == 7])); min(mod_access_edu_bmi$dev.ratio[mod_access_edu_bmi$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_edu_bmi$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_edu_bmi$beta[, i] != 0))$names
  new_deviance <- mod_access_edu_bmi$dev.ratio[i]
  new_varNum <- mod_access_edu_bmi$df[i]
  new_lambda <- mod_access_edu_bmi$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "int2"] <- "int2/green_dec23"
df$Variables[df$Variables == "int1"] <- "int1/int_accum/inc12_int/dec23_int"
df$Variables[df$Variables == "crit1"] <- "crit1/crit2"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "edu", ]
df <- df[df$Variables != "green_dec23", ]
df <- df[df$Variables != "int_accum", ]
df <- df[df$Variables != "green_inc12_int", ]
df <- df[df$Variables != "green_dec23_int", ]
df <- df[df$Variables != "crit2", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. This makes it clearer when different variables were added (not the neatest plot, though, as variables bunch up and overlap as lambda approaches 0...). Also, need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_edu_bmi$lambda, mod_access_edu_bmi$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_bmi$lambda)), ylim = c(0, max(mod_access_edu_bmi$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_edu_bmi$log_lambda <- log(mod_access_edu_bmi$lambda)
mod_access_edu_bmi

df$log_lambda <- log(as.numeric(df$Lambda))
df


# This looks better, although one misleading aspect is that lasso works by initialising the lambda value just above the threshold where no variables are entered (excluding covariates restrained to be in the model). This means that there will always be a 'first' variable entered early in the model, making it seem like this is the best fit to the data. However, because there always *has* to be one variable entered first, it doesn't mean that this is actually predictive of the outcome. In the plot here, even though 'green_inc23' was entered first, the actual model fit improvement over the null/covariate-only model is minimal (0.007% increase in deviance ratio), which is essentially 0, suggesting there is little/no association between access to green space and BMI.
plot(mod_access_edu_bmi$log_lambda, mod_access_edu_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_bmi$log_lambda)), ylim = c(0.007, max(mod_access_edu_bmi$dev.ratio)))
text(df$log_lambda, 0.007, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessEduBMI.pdf", height = 7, width = 11)
plot(mod_access_edu_bmi$log_lambda, mod_access_edu_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_bmi$log_lambda)), ylim = c(0.007, max(mod_access_edu_bmi$dev.ratio)))
text(df$log_lambda, 0.007, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter added (green_inc23) increases model fit of standard linear regression model
base_mod <- lm(data_access_edu_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "edu"])
summary(base_mod)

param1_mod <- lm(data_access_edu_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "edu"] + x_hypos[, "green_inc23"])
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between access to green space in childhood and BMI at age 7.
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model
mod.cv <- cv.glmnet(x_hypos, data_access_edu_bmi$BMI_f7, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters runs along the top of the plot)
plot(mod.cv) 

# The 1SE model contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)

# The optimal' model gives the same result (and as this optimal lasso is mainly for prediction, one may expect an increase in model complexity, suggesting again that green space effects are essentially non-existant)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between access to green space in childhood and BMI at age 7 (when controlling for ethnicity and parental education).


#### As a sensitivity analysis, will check whether get same results using just 2 time-points (birth and age 7).

data_access2 <- data

## Encode the hypotheses

# Critical period at first time point only
data_access2$crit1 <- data_access2$greenSpace_preg

# Interaction between SEP and first time point
data_access2$int1 <- data_access2$edu * data_access2$greenSpace_preg

# Critical period at second time point only
data_access2$crit2 <- data_access2$greenSpace_7

# Interaction between SEP and second time point
data_access2$int2 <- data_access2$edu * data_access2$greenSpace_7

# Linear accumulation of all exposures
data_access2$accumulation <- data_access2$greenSpace_preg + data_access2$greenSpace_7

# Interaction between SEP and cumulative exposure
data_access2$int_accum <- data_access2$edu * data_access2$accumulation

# Increase in access to green space from time 1 to time 2
data_access2$green_inc12 <- (1 - data_access2$greenSpace_preg) * data_access2$greenSpace_7

# Increase in access to green space from time 1 to time 2, with an interaction with SEP
data_access2$green_inc12_int <- (1 - data_access2$greenSpace_preg) * data_access2$greenSpace_7 * data_access2$edu

# Decrease in access to green space from time 1 to time 2
data_access2$green_dec12 <- (1 - data_access2$greenSpace_7) * data_access2$greenSpace_preg

# Decrease in access to green space from time 1 to time 2, with an interaction with SEP
data_access2$green_dec12_int <- (1 - data_access2$greenSpace_7) * data_access2$greenSpace_preg * data_access2$edu


## Make a dataset just with the outcomes, exposure hypotheses and covariates
data_access2_edu <- data_access2 %>%
  select(BMI_f7, BMI_f7_z, sysBP, diaBP, 
         age_f7, male, white, edu, 
         crit1, int1, crit2, int2, accumulation, int_accum, green_inc12, green_inc12_int, 
         green_dec12, green_dec12_int)


## Now analyse each outcome in turn, starting with BMI

# Reduce dataset down to just complete cases
data_access2_edu_bmi <- data_access2_edu %>%
  select(-BMI_f7_z, -sysBP, -diaBP) %>%
  filter(complete.cases(BMI_f7, age_f7, male, white, edu, crit1, int1))

summary(data_access2_edu_bmi)

# Save the life-course hypotheses and covariates as a matrix
x_hypos2 <- data_access2_edu_bmi %>%
  select(-BMI_f7)

x_hypos2 <- as.matrix(x_hypos2)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_bmi2 <- glmnet(x_hypos2, data_access2_edu_bmi$BMI_f7, 
                             alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos2) - 4))))

mod_access_edu_bmi2

# Look at model coefficients
coef(mod_access_edu_bmi2, s = max(mod_access_edu_bmi2$lambda[mod_access_edu_bmi2$df == 4])); min(mod_access_edu_bmi2$dev.ratio[mod_access_edu_bmi2$df == 4])

coef(mod_access_edu_bmi2, s = max(mod_access_edu_bmi2$lambda[mod_access_edu_bmi2$df == 6])); min(mod_access_edu_bmi2$dev.ratio[mod_access_edu_bmi2$df == 6])

# Increase in model fit when first variables entered is again minimal (0.003% compared to baseline), suggesting little association - Get same conclusion if use LR test or cross-validated lasso
(mod_access_edu_bmi2$dev.ratio[2] - mod_access_edu_bmi2$dev.ratio[1]) * 100

base_mod2 <- lm(data_access2_edu_bmi$BMI_f7 ~ x_hypos2[, "age_f7"] + x_hypos2[, "male"] + x_hypos2[, "white"] + 
                 x_hypos2[, "edu"])
summary(base_mod2)

param1_mod2 <- lm(data_access2_edu_bmi$BMI_f7 ~ x_hypos2[, "age_f7"] + x_hypos2[, "male"] + x_hypos2[, "white"] + 
                   x_hypos2[, "edu"] + x_hypos2[, "crit2"]  + x_hypos2[, "int2"])
summary(param1_mod2)

anova(base_mod2, param1_mod2)


# Alternative method using cross-validated lasso (again, both give null/covariate only model as best fit)
mod.cv2 <- cv.glmnet(x_hypos2, data_access2_edu_bmi$BMI_f7, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos2) - 4))))
mod.cv2

coef(mod.cv2, s = mod.cv2$lambda.1se)
coef(mod.cv2, s = mod.cv2$lambda.min)


## Quick check whether is any association with green space access is not include any of the covariates (using the 3 time-point model)
x_hypos_noCovars <- x_hypos[,!colnames(x_hypos) %in% c("age_f7", "male", "white", "edu")]
head(x_hypos_noCovars)

mod_access_edu_bmi_nocovs <- glmnet(x_hypos_noCovars, data_access_edu_bmi$BMI_f7, alpha = 1)
mod_access_edu_bmi_nocovs
#coef(mod_access_edu_bmi_nocovs)

coef(mod_access_edu_bmi_nocovs, s = max(mod_access_edu_bmi_nocovs$lambda[mod_access_edu_bmi_nocovs$df == 1])); min(mod_access_edu_bmi_nocovs$dev.ratio[mod_access_edu_bmi_nocovs$df == 1])


# Even without other covariates, improvement in model fit is still pretty minimal (0.009%, relative to null model)
(mod_access_edu_bmi_nocovs$dev.ratio[2] - mod_access_edu_bmi_nocovs$dev.ratio[1]) * 100

# Again, same results if use LR test and cross-validated lasso
base_mod_nocovs <- lm(data_access_edu_bmi$BMI_f7 ~ 1)
summary(base_mod_nocovs)

param1_mod_nocovs <- lm(data_access_edu_bmi$BMI_f7 ~ x_hypos_noCovars[, "green_inc23"])
summary(param1_mod_nocovs)

anova(base_mod_nocovs, param1_mod_nocovs)


# Alternative method using cross-validated lasso (again, both give null model as best fit)
mod.cv_nocovs <- cv.glmnet(x_hypos_noCovars, data_access2_edu_bmi$BMI_f7, alpha = 1)
mod.cv_nocovs

coef(mod.cv_nocovs, s = mod.cv_nocovs$lambda.1se)
coef(mod.cv_nocovs, s = mod.cv_nocovs$lambda.min)



#### Next, will explore binary access to green space exposure, with education as covariate/interaction term, and binary overweight variable as the outcome

# Reduce dataset down to just complete cases
data_access_edu_overweight <- data_access_edu %>%
  select(-BMI_f7, -sysBP, -diaBP) %>%
  filter(complete.cases(overweight, age_f7, male, white, edu, crit1, int1))

summary(data_access_edu_overweight)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_edu_overweight %>%
  select(-overweight)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_over <- glmnet(x_hypos, data_access_edu_overweight$overweight, family = "binomial",
                             alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_edu_over

# Plot these results
plot(mod_access_edu_over)

# Look at the variables included at each step
coef(mod_access_edu_over, s = max(mod_access_edu_over$lambda[mod_access_edu_over$df == 4])); min(mod_access_edu_over$dev.ratio[mod_access_edu_over$df == 4])

coef(mod_access_edu_over, s = max(mod_access_edu_over$lambda[mod_access_edu_over$df == 6])); min(mod_access_edu_over$dev.ratio[mod_access_edu_over$df == 6])

coef(mod_access_edu_over, s = max(mod_access_edu_over$lambda[mod_access_edu_over$df == 7])); min(mod_access_edu_over$dev.ratio[mod_access_edu_over$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_edu_over$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_edu_over$beta[, i] != 0))$names
  new_deviance <- mod_access_edu_over$dev.ratio[i]
  new_varNum <- mod_access_edu_over$df[i]
  new_lambda <- mod_access_edu_over$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "crit2"] <- "crit2/green_dec12_int"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "edu", ]
df <- df[df$Variables != "green_dec12_int", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_edu_over$lambda, mod_access_edu_over$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_over$lambda)), ylim = c(0, max(mod_access_edu_over$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_edu_over$log_lambda <- log(mod_access_edu_over$lambda)
mod_access_edu_over

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_access_edu_over$log_lambda, mod_access_edu_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_over$log_lambda)), ylim = c(0.001, max(mod_access_edu_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessEduOverweight.pdf", height = 7, width = 11)
plot(mod_access_edu_over$log_lambda, mod_access_edu_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_over$log_lambda)), ylim = c(0.001, max(mod_access_edu_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (crit2 and green_dec12_int) increases model fit of standard logistic regression model
base_mod <- glm(data_access_edu_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "edu"], family = "binomial")
summary(base_mod)

param1_mod <- glm(data_access_edu_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                  x_hypos[, "white"] + x_hypos[, "edu"] + x_hypos[, "crit2"] + x_hypos[, "green_dec12_int"],
                  family = "binomial")
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between access to green space in childhood and being overweight at age 7.
anova(base_mod, param1_mod, test = "Chisq")


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
  select(-BMI_f7, -overweight, -diaBP) %>%
  filter(complete.cases(sysBP, age_f7, male, white, edu, crit1, int1))

summary(data_access_edu_sysBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_edu_sysBP %>%
  select(-sysBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_sysBP <- glmnet(x_hypos, data_access_edu_sysBP$sysBP,
                              alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_edu_sysBP

# Plot these results
plot(mod_access_edu_sysBP)

# Look at the variables included at each step
coef(mod_access_edu_sysBP, s = max(mod_access_edu_sysBP$lambda[mod_access_edu_sysBP$df == 4])); min(mod_access_edu_sysBP$dev.ratio[mod_access_edu_sysBP$df == 4])

coef(mod_access_edu_sysBP, s = max(mod_access_edu_sysBP$lambda[mod_access_edu_sysBP$df == 5])); min(mod_access_edu_sysBP$dev.ratio[mod_access_edu_sysBP$df == 5])

coef(mod_access_edu_sysBP, s = max(mod_access_edu_sysBP$lambda[mod_access_edu_sysBP$df == 6])); min(mod_access_edu_sysBP$dev.ratio[mod_access_edu_sysBP$df == 6])

coef(mod_access_edu_sysBP, s = max(mod_access_edu_sysBP$lambda[mod_access_edu_sysBP$df == 7])); min(mod_access_edu_sysBP$dev.ratio[mod_access_edu_sysBP$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_edu_sysBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_edu_sysBP$beta[, i] != 0))$names
  new_deviance <- mod_access_edu_sysBP$dev.ratio[i]
  new_varNum <- mod_access_edu_sysBP$df[i]
  new_lambda <- mod_access_edu_sysBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "green_inc12_int"] <- "green_inc12_int/green_dec23"
df$Variables[df$Variables == "crit1"] <- "crit1/green_inc12"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "edu", ]
df <- df[df$Variables != "green_dec23", ]
df <- df[df$Variables != "green_inc12", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_edu_sysBP$lambda, mod_access_edu_sysBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_sysBP$lambda)), ylim = c(0, max(mod_access_edu_sysBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_edu_sysBP$log_lambda <- log(mod_access_edu_sysBP$lambda)
mod_access_edu_sysBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_access_edu_sysBP$log_lambda, mod_access_edu_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_sysBP$log_lambda)), ylim = c(0.005, max(mod_access_edu_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessEduSysBP.pdf", height = 7, width = 11)
plot(mod_access_edu_sysBP$log_lambda, mod_access_edu_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_sysBP$log_lambda)), ylim = c(0.005, max(mod_access_edu_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (green_dec12) increases model fit of standard linear regression model
base_mod <- lm(data_access_edu_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                  x_hypos[, "edu"])
summary(base_mod)

param1_mod <- lm(data_access_edu_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                    x_hypos[, "white"] + x_hypos[, "edu"] + x_hypos[, "green_dec12"])
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between access to green space in childhood and being overweight at age 7.
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_edu_sysBP$sysBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between access to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and parental education).



#### Next, will explore binary access to green space exposure, with education as covariate/interaction term, and continuous diastolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_access_edu_diaBP <- data_access_edu %>%
  select(-BMI_f7, -overweight, -sysBP) %>%
  filter(complete.cases(diaBP, age_f7, male, white, edu, crit1, int1))

summary(data_access_edu_diaBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_edu_diaBP %>%
  select(-diaBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_edu_diaBP <- glmnet(x_hypos, data_access_edu_diaBP$diaBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_edu_diaBP

# Plot these results
plot(mod_access_edu_diaBP)

# Look at the variables included at each step
coef(mod_access_edu_diaBP, s = max(mod_access_edu_diaBP$lambda[mod_access_edu_diaBP$df == 4])); min(mod_access_edu_diaBP$dev.ratio[mod_access_edu_diaBP$df == 4])

coef(mod_access_edu_diaBP, s = max(mod_access_edu_diaBP$lambda[mod_access_edu_diaBP$df == 5])); min(mod_access_edu_diaBP$dev.ratio[mod_access_edu_diaBP$df == 5])

coef(mod_access_edu_diaBP, s = max(mod_access_edu_diaBP$lambda[mod_access_edu_diaBP$df == 6])); min(mod_access_edu_diaBP$dev.ratio[mod_access_edu_diaBP$df == 6])

coef(mod_access_edu_diaBP, s = max(mod_access_edu_diaBP$lambda[mod_access_edu_diaBP$df == 8])); min(mod_access_edu_diaBP$dev.ratio[mod_access_edu_diaBP$df == 8])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_edu_diaBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_edu_diaBP$beta[, i] != 0))$names
  new_deviance <- mod_access_edu_diaBP$dev.ratio[i]
  new_varNum <- mod_access_edu_diaBP$df[i]
  new_lambda <- mod_access_edu_diaBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "green_inc12"] <- "green_inc12/green_dec23"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "edu", ]
df <- df[df$Variables != "green_dec23", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_edu_diaBP$lambda, mod_access_edu_diaBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_diaBP$lambda)), ylim = c(0, max(mod_access_edu_diaBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_edu_diaBP$log_lambda <- log(mod_access_edu_diaBP$lambda)
mod_access_edu_diaBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome). Here, variable 'green_dec12' is entered quite a bit before all other variables, and although the deviance ratio hasn't improved by much (0.02%), it's a larger increase than seen for the models above (normally <0.008%).
plot(mod_access_edu_diaBP$log_lambda, mod_access_edu_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_diaBP$log_lambda)), ylim = c(0.004, max(mod_access_edu_diaBP$dev.ratio)))
text(df$log_lambda, 0.004, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessEduDiaBP.pdf", height = 7, width = 11)
plot(mod_access_edu_diaBP$log_lambda, mod_access_edu_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_edu_diaBP$log_lambda)), ylim = c(0.004, max(mod_access_edu_diaBP$dev.ratio)))
text(df$log_lambda, 0.004, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (green_dec12) increases model fit of standard linear regression model
base_mod <- lm(data_access_edu_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "edu"])
summary(base_mod)

param1_mod <- lm(data_access_edu_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "edu"] + x_hypos[, "green_dec12"])
summary(param1_mod)

# Does seem to be an association here, with a decrease in green space access between pregnancy and age 4 associated with lower diastolic BP. However, given that the effect size is tiny and interpretation not clear (why would reducing access to green space be associated with lower BP?), this could just be the model hooking up to random noise.
anova(base_mod, param1_mod)

# Does adding the next parameter (green_inc23_int) improve model fit? No.
param2_mod <- lm(data_access_edu_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "edu"] + x_hypos[, "green_dec12"] + + x_hypos[, "green_inc23_int"])
summary(param2_mod)

anova(param1_mod, param2_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_edu_diaBP$diaBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - The 1SE model is just the null/covariate-only model, while the 'optimal'/lowest MSE model also contains green_dec12. The difference in MSE is minimal, though, so not especially convincing. Suggests at best a very weak association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between access to green space in childhood and diastolic blood pressure at age 7 (when controlling for ethnicity and parental education), although potentially decrease in access to green space between preg and age 4 linked to lower diastolic BP.



##############################################################################################################
#### Next, explore whether binary access to green space as exposure associated with cardiometabolic outcomes with deprivation as covariate/interaction term (instead of parental education)

### Derive the different life-course variables and encode as hypotheses
table(data$greenSpace_preg, useNA = "ifany")
table(data$greenSpace_4, useNA = "ifany")
table(data$greenSpace_7, useNA = "ifany")
table(data$deprived, useNA = "ifany")

# Make a dataset just for this exposure
data_access <- data


## Here deprivation is the SEP covariate/interaction term

# Critical period at first time point only
data_access$crit1 <- data_access$greenSpace_preg

# Interaction between SEP and first time point
data_access$int1 <- data_access$deprived * data_access$greenSpace_preg

# Critical period at second time point only
data_access$crit2 <- data_access$greenSpace_4

# Interaction between SEP and second time point
data_access$int2 <- data_access$deprived * data_access$greenSpace_4

# Critical period at third time point only
data_access$crit3 <- data_access$greenSpace_7

# Interaction between SEP and third time point
data_access$int3 <- data_access$deprived * data_access$greenSpace_7

# Linear accumulation of all exposures
data_access$accumulation <- data_access$greenSpace_preg + data_access$greenSpace_4 + data_access$greenSpace_7

# Interaction between SEP and cumulative exposure
data_access$int_accum <- data_access$deprived * data_access$accumulation

# Increase in access to green space from time 1 to time 2
data_access$green_inc12 <- (1 - data_access$greenSpace_preg) * data_access$greenSpace_4

# Increase in access to green space from time 1 to time 2, with an interaction with SEP
data_access$green_inc12_int <- (1 - data_access$greenSpace_preg) * data_access$greenSpace_4 * data_access$deprived

# Decrease in access to green space from time 1 to time 2
data_access$green_dec12 <- (1 - data_access$greenSpace_4) * data_access$greenSpace_preg

# Decrease in access to green space from time 1 to time 2, with an interaction with SEP
data_access$green_dec12_int <- (1 - data_access$greenSpace_4) * data_access$greenSpace_preg * data_access$deprived

# Increase in access to green space from time 2 to time 3
data_access$green_inc23 <- (1 - data_access$greenSpace_4) * data_access$greenSpace_7

# Increase in access to green space from time 2 to time 3, with an interaction with SEP
data_access$green_inc23_int <- (1 - data_access$greenSpace_4) * data_access$greenSpace_7 * data_access$deprived

# Decrease in access to green space from time 2 to time 3
data_access$green_dec23 <- (1 - data_access$greenSpace_7) * data_access$greenSpace_4

# Decrease in access to green space from time 2 to time 3, with an interaction with SEP
data_access$green_dec23_int <- (1 - data_access$greenSpace_7) * data_access$greenSpace_4 * data_access$deprived


## Make a dataset just with the outcomes, exposure hypotheses and covariates
data_access_dep <- data_access %>%
  select(BMI_f7, overweight, sysBP, diaBP, 
         age_f7, male, white, deprived, 
         crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_inc12, green_inc12_int, 
         green_dec12, green_dec12_int, green_inc23, green_inc23_int, green_dec23, green_dec23_int)


## Now analyse each outcome in turn, starting with BMI

# Reduce dataset down to just complete cases
data_access_dep_bmi <- data_access_dep %>%
  select(-overweight, -sysBP, -diaBP) %>%
  filter(complete.cases(BMI_f7, age_f7, male, white, deprived, crit1, int1))

summary(data_access_dep_bmi)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_dep_bmi %>%
  select(-BMI_f7)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (dep, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_dep_bmi <- glmnet(x_hypos, data_access_dep_bmi$BMI_f7, 
                             alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_dep_bmi

# Plot these results
plot(mod_access_dep_bmi)

# Look at the variables included at each step
coef(mod_access_dep_bmi, s = max(mod_access_dep_bmi$lambda[mod_access_dep_bmi$df == 4])); min(mod_access_dep_bmi$dev.ratio[mod_access_dep_bmi$df == 4])

coef(mod_access_dep_bmi, s = max(mod_access_dep_bmi$lambda[mod_access_dep_bmi$df == 5])); min(mod_access_dep_bmi$dev.ratio[mod_access_dep_bmi$df == 5])

coef(mod_access_dep_bmi, s = max(mod_access_dep_bmi$lambda[mod_access_dep_bmi$df == 6])); min(mod_access_dep_bmi$dev.ratio[mod_access_dep_bmi$df == 6])

coef(mod_access_dep_bmi, s = max(mod_access_dep_bmi$lambda[mod_access_dep_bmi$df == 7])); min(mod_access_dep_bmi$dev.ratio[mod_access_dep_bmi$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_dep_bmi$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_dep_bmi$beta[, i] != 0))$names
  new_deviance <- mod_access_dep_bmi$dev.ratio[i]
  new_varNum <- mod_access_dep_bmi$df[i]
  new_lambda <- mod_access_dep_bmi$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "crit2"] <- "crit2/green_dec23_int"
df$Variables[df$Variables == "green_inc12"] <- "green_inc12/green_dec12"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "deprived", ]
df <- df[df$Variables != "green_dec23_int", ]
df <- df[df$Variables != "green_dec12", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. This makes it clearer when different variables were added (not the neatest plot, though, as variables bunch up and overlap as lambda approaches 0...). Also, need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_dep_bmi$lambda, mod_access_dep_bmi$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_bmi$lambda)), ylim = c(0, max(mod_access_dep_bmi$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_dep_bmi$log_lambda <- log(mod_access_dep_bmi$lambda)
mod_access_dep_bmi

df$log_lambda <- log(as.numeric(df$Lambda))
df


# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_access_dep_bmi$log_lambda, mod_access_dep_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_bmi$log_lambda)), ylim = c(0.008, max(mod_access_dep_bmi$dev.ratio)))
text(df$log_lambda, 0.008, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessDepBMI.pdf", height = 7, width = 11)
plot(mod_access_dep_bmi$log_lambda, mod_access_dep_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_bmi$log_lambda)), ylim = c(0.008, max(mod_access_dep_bmi$dev.ratio)))
text(df$log_lambda, 0.008, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter added (green_inc23) increases model fit of standard linear regression model
base_mod <- lm(data_access_dep_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "deprived"])
summary(base_mod)

param1_mod <- lm(data_access_dep_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                   x_hypos[, "deprived"] + x_hypos[, "green_inc23"])
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between access to green space in childhood and BMI at age 7.
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model
mod.cv <- cv.glmnet(x_hypos, data_access_dep_bmi$BMI_f7, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at the top of the plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between access to green space in childhood and BMI at age 7 (when controlling for ethnicity and household deprivation).



#### Next, will explore binary access to green space exposure, with deprivation as covariate/interaction term, and binary overweight variable as the outcome

# Reduce dataset down to just complete cases
data_access_dep_overweight <- data_access_dep %>%
  select(-BMI_f7, -sysBP, -diaBP) %>%
  filter(complete.cases(overweight, age_f7, male, white, deprived, crit1, int1))

summary(data_access_dep_overweight)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_dep_overweight %>%
  select(-overweight)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_dep_over <- glmnet(x_hypos, data_access_dep_overweight$overweight, family = "binomial",
                              alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_dep_over

# Plot these results
plot(mod_access_dep_over)

# Look at the variables included at each step
coef(mod_access_dep_over, s = max(mod_access_dep_over$lambda[mod_access_dep_over$df == 4])); min(mod_access_dep_over$dev.ratio[mod_access_dep_over$df == 4])

coef(mod_access_dep_over, s = max(mod_access_dep_over$lambda[mod_access_dep_over$df == 6])); min(mod_access_dep_over$dev.ratio[mod_access_dep_over$df == 6])

coef(mod_access_dep_over, s = max(mod_access_dep_over$lambda[mod_access_dep_over$df == 7])); min(mod_access_dep_over$dev.ratio[mod_access_dep_over$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_dep_over$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_dep_over$beta[, i] != 0))$names
  new_deviance <- mod_access_dep_over$dev.ratio[i]
  new_varNum <- mod_access_dep_over$df[i]
  new_lambda <- mod_access_dep_over$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "crit2"] <- "crit2/green_inc12_int"
df$Variables[df$Variables == "int3"] <- "int3/accumulation/green_inc12"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "deprived", ]
df <- df[df$Variables != "green_inc12_int", ]
df <- df[df$Variables != "accumulation", ]
df <- df[df$Variables != "green_inc12", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_dep_over$lambda, mod_access_dep_over$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_over$lambda)), ylim = c(0, max(mod_access_dep_over$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_dep_over$log_lambda <- log(mod_access_dep_over$lambda)
mod_access_dep_over

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_access_dep_over$log_lambda, mod_access_dep_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_over$log_lambda)), ylim = c(0.001, max(mod_access_dep_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessDepOverweight.pdf", height = 7, width = 11)
plot(mod_access_dep_over$log_lambda, mod_access_dep_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_over$log_lambda)), ylim = c(0.001, max(mod_access_dep_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (crit2 and green_inc12_int) increases model fit of standard logistic regression model
base_mod <- glm(data_access_dep_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                  x_hypos[, "deprived"], family = "binomial")
summary(base_mod)

param1_mod <- glm(data_access_dep_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                    x_hypos[, "white"] + x_hypos[, "deprived"] + x_hypos[, "crit2"] + x_hypos[, "green_inc12_int"],
                  family = "binomial")
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between access to green space in childhood and being overweight at age 7.
anova(base_mod, param1_mod, test = "Chisq")


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_dep_overweight$overweight, family = "binomial", 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between access to green space in childhood and being overweight at age 7 (when controlling for ethnicity and household deprivation).



#### Next, will explore binary access to green space exposure, with deprivation as covariate/interaction term, and continuous systolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_access_dep_sysBP <- data_access_dep %>%
  select(-BMI_f7, -overweight, -diaBP) %>%
  filter(complete.cases(sysBP, age_f7, male, white, deprived, crit1, int1))

summary(data_access_dep_sysBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_dep_sysBP %>%
  select(-sysBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_dep_sysBP <- glmnet(x_hypos, data_access_dep_sysBP$sysBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_dep_sysBP

# Plot these results
plot(mod_access_dep_sysBP)

# Look at the variables included at each step
coef(mod_access_dep_sysBP, s = max(mod_access_dep_sysBP$lambda[mod_access_dep_sysBP$df == 4])); min(mod_access_dep_sysBP$dev.ratio[mod_access_dep_sysBP$df == 4])

coef(mod_access_dep_sysBP, s = max(mod_access_dep_sysBP$lambda[mod_access_dep_sysBP$df == 5])); min(mod_access_dep_sysBP$dev.ratio[mod_access_dep_sysBP$df == 5])

coef(mod_access_dep_sysBP, s = max(mod_access_dep_sysBP$lambda[mod_access_dep_sysBP$df == 6])); min(mod_access_dep_sysBP$dev.ratio[mod_access_dep_sysBP$df == 6])

coef(mod_access_dep_sysBP, s = max(mod_access_dep_sysBP$lambda[mod_access_dep_sysBP$df == 8])); min(mod_access_dep_sysBP$dev.ratio[mod_access_dep_sysBP$df == 8])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_dep_sysBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_dep_sysBP$beta[, i] != 0))$names
  new_deviance <- mod_access_dep_sysBP$dev.ratio[i]
  new_varNum <- mod_access_dep_sysBP$df[i]
  new_lambda <- mod_access_dep_sysBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "accumulation"] <- "accumulation/green_dec23"
df$Variables[df$Variables == "int1"] <- "int1/int2/crit3"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "deprived", ]
df <- df[df$Variables != "green_dec23", ]
df <- df[df$Variables != "int2", ]
df <- df[df$Variables != "crit3", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_dep_sysBP$lambda, mod_access_dep_sysBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_sysBP$lambda)), ylim = c(0, max(mod_access_dep_sysBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_dep_sysBP$log_lambda <- log(mod_access_dep_sysBP$lambda)
mod_access_dep_sysBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_access_dep_sysBP$log_lambda, mod_access_dep_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_sysBP$log_lambda)), ylim = c(0.005, max(mod_access_dep_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessDepSysBP.pdf", height = 7, width = 11)
plot(mod_access_dep_sysBP$log_lambda, mod_access_dep_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_sysBP$log_lambda)), ylim = c(0.005, max(mod_access_dep_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that not much is really going on here, although first variable eneterd ('green_dec12_int') does appear to increase model fit marginally... - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (green_dec12) increases model fit of standard linear regression model
base_mod <- lm(data_access_dep_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "deprived"])
summary(base_mod)

param1_mod <- lm(data_access_dep_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "deprived"] + x_hypos[, "green_dec12_int"])
summary(param1_mod)

# Does seem to be an association here, with a decrease in green space access between pregnancy and age 4 in deprived individuals associated with lower diastolic BP. However, given that the effect size is tiny and interpretation not clear (why would reducing access to green space be associated with lower BP among deprived individuals?), this could just be the model hooking up to random noise.
anova(base_mod, param1_mod)

# Does adding the next parameter (crit2) improve model fit? No.
param2_mod <- lm(data_access_dep_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "deprived"] + x_hypos[, "green_dec12_int"] + x_hypos[, "crit2"])
summary(param2_mod)

anova(param1_mod, param2_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_dep_sysBP$sysBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - The 1SE model is just the null/covariate-only model, while the 'optimal'/lowest MSE model also contains green_dec12_int. The difference in MSE is minimal, though, so not especially convincing. Suggests at best a very weak association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between access to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and household deprivation), although potentially decrease in access to green space between preg and age 4 linked to lower systolic BP, but only in those from more deprived areas.



#### Next, will explore binary access to green space exposure, with deprivation as covariate/interaction term, and continuous diastolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_access_dep_diaBP <- data_access_dep %>%
  select(-BMI_f7, -overweight, -sysBP) %>%
  filter(complete.cases(diaBP, age_f7, male, white, deprived, crit1, int1))

summary(data_access_dep_diaBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_dep_diaBP %>%
  select(-diaBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_dep_diaBP <- glmnet(x_hypos, data_access_dep_diaBP$diaBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_dep_diaBP

# Plot these results
plot(mod_access_dep_diaBP)

# Look at the variables included at each step
coef(mod_access_dep_diaBP, s = max(mod_access_dep_diaBP$lambda[mod_access_dep_diaBP$df == 4])); min(mod_access_dep_diaBP$dev.ratio[mod_access_dep_diaBP$df == 4])

coef(mod_access_dep_diaBP, s = max(mod_access_dep_diaBP$lambda[mod_access_dep_diaBP$df == 5])); min(mod_access_dep_diaBP$dev.ratio[mod_access_dep_diaBP$df == 5])

coef(mod_access_dep_diaBP, s = max(mod_access_dep_diaBP$lambda[mod_access_dep_diaBP$df == 7])); min(mod_access_dep_diaBP$dev.ratio[mod_access_dep_diaBP$df == 7])

coef(mod_access_dep_diaBP, s = max(mod_access_dep_diaBP$lambda[mod_access_dep_diaBP$df == 8])); min(mod_access_dep_diaBP$dev.ratio[mod_access_dep_diaBP$df == 8])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_dep_diaBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_dep_diaBP$beta[, i] != 0))$names
  new_deviance <- mod_access_dep_diaBP$dev.ratio[i]
  new_varNum <- mod_access_dep_diaBP$df[i]
  new_lambda <- mod_access_dep_diaBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "green_inc12_int"] <- "green_inc12_int/green_dec12_int"
df$Variables[df$Variables == "int2"] <- "int2/int3"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "deprived", ]
df <- df[df$Variables != "green_dec12_int", ]
df <- df[df$Variables != "int3", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_dep_diaBP$lambda, mod_access_dep_diaBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_diaBP$lambda)), ylim = c(0, max(mod_access_dep_diaBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_dep_diaBP$log_lambda <- log(mod_access_dep_diaBP$lambda)
mod_access_dep_diaBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome). Here, variable 'green_dec12' is entered quite a bit before all other variables, and although the deviance ratio hasn't improved by much (0.02%), it's a larger increase than seen for the models above (normally <0.008%).
plot(mod_access_dep_diaBP$log_lambda, mod_access_dep_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_diaBP$log_lambda)), ylim = c(0.002, max(mod_access_dep_diaBP$dev.ratio)))
text(df$log_lambda, 0.002, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessDepDiaBP.pdf", height = 7, width = 11)
plot(mod_access_dep_diaBP$log_lambda, mod_access_dep_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_diaBP$log_lambda)), ylim = c(0.002, max(mod_access_dep_diaBP$dev.ratio)))
text(df$log_lambda, 0.002, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (green_dec12) increases model fit of standard linear regression model
base_mod <- lm(data_access_dep_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "deprived"])
summary(base_mod)

param1_mod <- lm(data_access_dep_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "deprived"] + x_hypos[, "green_dec12"])
summary(param1_mod)

# Does seem to be an association here, with a decrease in green space access between pregnancy and age 4 associated with lower diastolic BP. However, given that the effect size is tiny and interpretation not clear (why would reducing access to green space be associated with lower BP?), this could just be the model hooking up to random noise.
anova(base_mod, param1_mod)

# Does adding the next parameter(s) (green_inc12_int and green_dec12_int) improve model fit? No.
param2_mod <- lm(data_access_dep_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "deprived"] + x_hypos[, "green_dec12"] + 
                   x_hypos[, "green_inc23_int"] + x_hypos[, "green_dec12_int"])
summary(param2_mod)

anova(param1_mod, param2_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_dep_diaBP$diaBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - The 1SE model is just the null/covariate-only model, while the 'optimal'/lowest MSE model also contains green_dec12, green_inc12_int and green_dec12_int. The difference in MSE is minimal, though, so not especially convincing. Suggests at best a very weak association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between access to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and household deprivation), although potentially decrease in access to green space between preg and age 4 linked to lower diastolic BP.




##############################################################################################################
#### Next, explore whether binary access to green space as exposure associated with cardiometabolic outcomes with disposable household income as covariate/interaction term (instead of parental education)

### Derive the different life-course variables and encode as hypotheses
table(data$greenSpace_preg, useNA = "ifany")
table(data$greenSpace_4, useNA = "ifany")
table(data$greenSpace_7, useNA = "ifany")
table(data$lowIncome, useNA = "ifany")

# Make a dataset just for this exposure
data_access <- data


## Here low household income is the SEP covariate/interaction term

# Critical period at first time point only
data_access$crit1 <- data_access$greenSpace_preg

# Interaction between SEP and first time point
data_access$int1 <- data_access$lowIncome * data_access$greenSpace_preg

# Critical period at second time point only
data_access$crit2 <- data_access$greenSpace_4

# Interaction between SEP and second time point
data_access$int2 <- data_access$lowIncome * data_access$greenSpace_4

# Critical period at third time point only
data_access$crit3 <- data_access$greenSpace_7

# Interaction between SEP and third time point
data_access$int3 <- data_access$lowIncome * data_access$greenSpace_7

# Linear accumulation of all exposures
data_access$accumulation <- data_access$greenSpace_preg + data_access$greenSpace_4 + data_access$greenSpace_7

# Interaction between SEP and cumulative exposure
data_access$int_accum <- data_access$lowIncome * data_access$accumulation

# Increase in access to green space from time 1 to time 2
data_access$green_inc12 <- (1 - data_access$greenSpace_preg) * data_access$greenSpace_4

# Increase in access to green space from time 1 to time 2, with an interaction with SEP
data_access$green_inc12_int <- (1 - data_access$greenSpace_preg) * data_access$greenSpace_4 * data_access$lowIncome

# Decrease in access to green space from time 1 to time 2
data_access$green_dec12 <- (1 - data_access$greenSpace_4) * data_access$greenSpace_preg

# Decrease in access to green space from time 1 to time 2, with an interaction with SEP
data_access$green_dec12_int <- (1 - data_access$greenSpace_4) * data_access$greenSpace_preg * data_access$lowIncome

# Increase in access to green space from time 2 to time 3
data_access$green_inc23 <- (1 - data_access$greenSpace_4) * data_access$greenSpace_7

# Increase in access to green space from time 2 to time 3, with an interaction with SEP
data_access$green_inc23_int <- (1 - data_access$greenSpace_4) * data_access$greenSpace_7 * data_access$lowIncome

# Decrease in access to green space from time 2 to time 3
data_access$green_dec23 <- (1 - data_access$greenSpace_7) * data_access$greenSpace_4

# Decrease in access to green space from time 2 to time 3, with an interaction with SEP
data_access$green_dec23_int <- (1 - data_access$greenSpace_7) * data_access$greenSpace_4 * data_access$lowIncome


## Make a dataset just with the outcomes, exposure hypotheses and covariates
data_access_inc <- data_access %>%
  select(BMI_f7, overweight, sysBP, diaBP, 
         age_f7, male, white, lowIncome, 
         crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_inc12, green_inc12_int, 
         green_dec12, green_dec12_int, green_inc23, green_inc23_int, green_dec23, green_dec23_int)


## Now analyse each outcome in turn, starting with BMI

# Reduce dataset down to just complete cases
data_access_inc_bmi <- data_access_inc %>%
  select(-overweight, -sysBP, -diaBP) %>%
  filter(complete.cases(BMI_f7, age_f7, male, white, lowIncome, crit1, int1))

summary(data_access_inc_bmi)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_inc_bmi %>%
  select(-BMI_f7)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (income, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_inc_bmi <- glmnet(x_hypos, data_access_inc_bmi$BMI_f7, 
                             alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_inc_bmi

# Plot these results
plot(mod_access_inc_bmi)

# Look at the variables included at each step
coef(mod_access_inc_bmi, s = max(mod_access_inc_bmi$lambda[mod_access_inc_bmi$df == 4])); min(mod_access_inc_bmi$dev.ratio[mod_access_inc_bmi$df == 4])

coef(mod_access_inc_bmi, s = max(mod_access_inc_bmi$lambda[mod_access_inc_bmi$df == 5])); min(mod_access_inc_bmi$dev.ratio[mod_access_inc_bmi$df == 5])

coef(mod_access_inc_bmi, s = max(mod_access_inc_bmi$lambda[mod_access_inc_bmi$df == 7])); min(mod_access_inc_bmi$dev.ratio[mod_access_inc_bmi$df == 7])

coef(mod_access_inc_bmi, s = max(mod_access_inc_bmi$lambda[mod_access_inc_bmi$df == 8])); min(mod_access_inc_bmi$dev.ratio[mod_access_inc_bmi$df == 8])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_inc_bmi$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_inc_bmi$beta[, i] != 0))$names
  new_deviance <- mod_access_inc_bmi$dev.ratio[i]
  new_varNum <- mod_access_inc_bmi$df[i]
  new_lambda <- mod_access_inc_bmi$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "crit1"] <- "crit1/green_dec23_int"
df$Variables[df$Variables == "int2"] <- "int2/int3"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "lowIncome", ]
df <- df[df$Variables != "green_dec23_int", ]
df <- df[df$Variables != "int3", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. This makes it clearer when different variables were added (not the neatest plot, though, as variables bunch up and overlap as lambda approaches 0...). Also, need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_inc_bmi$lambda, mod_access_inc_bmi$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_inc_bmi$lambda)), ylim = c(0, max(mod_access_inc_bmi$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_inc_bmi$log_lambda <- log(mod_access_inc_bmi$lambda)
mod_access_inc_bmi

df$log_lambda <- log(as.numeric(df$Lambda))
df


# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_access_inc_bmi$log_lambda, mod_access_inc_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_inc_bmi$log_lambda)), ylim = c(0.006, max(mod_access_inc_bmi$dev.ratio)))
text(df$log_lambda, 0.006, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessIncomeBMI.pdf", height = 7, width = 11)
plot(mod_access_inc_bmi$log_lambda, mod_access_inc_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_inc_bmi$log_lambda)), ylim = c(0.006, max(mod_access_inc_bmi$dev.ratio)))
text(df$log_lambda, 0.006, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter added (green_inc23) increases model fit of standard linear regression model
base_mod <- lm(data_access_inc_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "lowIncome"])
summary(base_mod)

param1_mod <- lm(data_access_inc_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                   x_hypos[, "lowIncome"] + x_hypos[, "green_inc23"])
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between access to green space in childhood and BMI at age 7.
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model
mod.cv <- cv.glmnet(x_hypos, data_access_inc_bmi$BMI_f7, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at the top of the plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between access to green space in childhood and BMI at age 7 (when controlling for ethnicity and household deprivation).



#### Next, will explore binary access to green space exposure, with income as covariate/interaction term, and binary overweight variable as the outcome

# Reduce dataset down to just complete cases
data_access_inc_overweight <- data_access_inc %>%
  select(-BMI_f7, -sysBP, -diaBP) %>%
  filter(complete.cases(overweight, age_f7, male, white, lowIncome, crit1, int1))

summary(data_access_inc_overweight)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_inc_overweight %>%
  select(-overweight)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_inc_over <- glmnet(x_hypos, data_access_inc_overweight$overweight, family = "binomial",
                              alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_inc_over

# Plot these results
plot(mod_access_inc_over)

# Look at the variables included at each step
coef(mod_access_inc_over, s = max(mod_access_inc_over$lambda[mod_access_inc_over$df == 4])); min(mod_access_inc_over$dev.ratio[mod_access_inc_over$df == 4])

coef(mod_access_inc_over, s = max(mod_access_inc_over$lambda[mod_access_inc_over$df == 6])); min(mod_access_inc_over$dev.ratio[mod_access_inc_over$df == 6])

coef(mod_access_inc_over, s = max(mod_access_inc_over$lambda[mod_access_inc_over$df == 7])); min(mod_access_inc_over$dev.ratio[mod_access_inc_over$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_inc_over$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_inc_over$beta[, i] != 0))$names
  new_deviance <- mod_access_inc_over$dev.ratio[i]
  new_varNum <- mod_access_inc_over$df[i]
  new_lambda <- mod_access_inc_over$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "crit1"] <- "crit1/green_inc12"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "lowIncome", ]
df <- df[df$Variables != "green_inc12", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_inc_over$lambda, mod_access_inc_over$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_inc_over$lambda)), ylim = c(0, max(mod_access_inc_over$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_inc_over$log_lambda <- log(mod_access_inc_over$lambda)
mod_access_inc_over

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome). (Have to add a bit of faff here to get y-axis to not display in scientific notation)
plot(mod_access_inc_over$log_lambda, mod_access_inc_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", yaxt = "n",
     xlim = rev(range(mod_access_inc_over$log_lambda)), ylim = c(0.000, max(mod_access_inc_over$dev.ratio)))
axis(2, at = c(0, 0.0002, 0.0004, 0.0006, 0.0008), labels = format(c(0, 0.0002, 0.0004, 0.0006, 0.0008), scientific = FALSE))
text(df$log_lambda, 0.000, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessIncomeOverweight.pdf", height = 7, width = 11)
plot(mod_access_inc_over$log_lambda, mod_access_inc_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", yaxt = "n",
     xlim = rev(range(mod_access_inc_over$log_lambda)), ylim = c(0.000, max(mod_access_inc_over$dev.ratio)))
axis(2, at = c(0, 0.0002, 0.0004, 0.0006, 0.0008), labels = format(c(0, 0.0002, 0.0004, 0.0006, 0.0008), scientific = FALSE))
text(df$log_lambda, 0.000, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (crit2 and green_inc12_int) increases model fit of standard logistic regression model
base_mod <- glm(data_access_inc_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                  x_hypos[, "lowIncome"], family = "binomial")
summary(base_mod)

param1_mod <- glm(data_access_inc_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                    x_hypos[, "white"] + x_hypos[, "lowIncome"] + x_hypos[, "crit1"] + x_hypos[, "green_inc12_int"],
                  family = "binomial")
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between access to green space in childhood and being overweight at age 7.
anova(base_mod, param1_mod, test = "Chisq")


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_inc_overweight$overweight, family = "binomial", 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between access to green space in childhood and being overweight at age 7 (when controlling for ethnicity and household income).



#### Next, will explore binary access to green space exposure, with income as covariate/interaction term, and continuous systolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_access_inc_sysBP <- data_access_inc %>%
  select(-BMI_f7, -overweight, -diaBP) %>%
  filter(complete.cases(sysBP, age_f7, male, white, lowIncome, crit1, int1))

summary(data_access_inc_sysBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_inc_sysBP %>%
  select(-sysBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (income, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_inc_sysBP <- glmnet(x_hypos, data_access_inc_sysBP$sysBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_inc_sysBP

# Plot these results
plot(mod_access_inc_sysBP)

# Look at the variables included at each step
coef(mod_access_inc_sysBP, s = max(mod_access_inc_sysBP$lambda[mod_access_inc_sysBP$df == 4])); min(mod_access_inc_sysBP$dev.ratio[mod_access_inc_sysBP$df == 4])

coef(mod_access_inc_sysBP, s = max(mod_access_inc_sysBP$lambda[mod_access_inc_sysBP$df == 6])); min(mod_access_inc_sysBP$dev.ratio[mod_access_inc_sysBP$df == 6])

coef(mod_access_inc_sysBP, s = max(mod_access_inc_sysBP$lambda[mod_access_inc_sysBP$df == 7])); min(mod_access_inc_sysBP$dev.ratio[mod_access_inc_sysBP$df == 7])

coef(mod_access_inc_sysBP, s = max(mod_access_inc_sysBP$lambda[mod_access_inc_sysBP$df == 8])); min(mod_access_inc_sysBP$dev.ratio[mod_access_inc_sysBP$df == 8])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_inc_sysBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_inc_sysBP$beta[, i] != 0))$names
  new_deviance <- mod_access_inc_sysBP$dev.ratio[i]
  new_varNum <- mod_access_inc_sysBP$df[i]
  new_lambda <- mod_access_inc_sysBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "green_dec12"] <- "green_dec12/green_dec12_int"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "lowIncome", ]
df <- df[df$Variables != "green_dec12_int", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_inc_sysBP$lambda, mod_access_inc_sysBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_inc_sysBP$lambda)), ylim = c(0, max(mod_access_inc_sysBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_inc_sysBP$log_lambda <- log(mod_access_inc_sysBP$lambda)
mod_access_inc_sysBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_access_inc_sysBP$log_lambda, mod_access_inc_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_inc_sysBP$log_lambda)), ylim = c(0.000, max(mod_access_inc_sysBP$dev.ratio)))
text(df$log_lambda, 0.000, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessIncSysBP.pdf", height = 7, width = 11)
plot(mod_access_inc_sysBP$log_lambda, mod_access_inc_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_inc_sysBP$log_lambda)), ylim = c(0.000, max(mod_access_inc_sysBP$dev.ratio)))
text(df$log_lambda, 0.000, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that not much is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (green_dec12 and green_dec12_int) increases model fit of standard linear regression model
base_mod <- lm(data_access_inc_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "lowIncome"])
summary(base_mod)

param1_mod <- lm(data_access_inc_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "lowIncome"] + x_hypos[, "green_dec12"] + 
                   x_hypos[, "green_dec12_int"])
summary(param1_mod)

# No real association here
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_inc_sysBP$sysBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv)

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between access to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and household income).



#### Next, will explore binary access to green space exposure, with income as covariate/interaction term, and continuous diastolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_access_inc_diaBP <- data_access_inc %>%
  select(-BMI_f7, -overweight, -sysBP) %>%
  filter(complete.cases(diaBP, age_f7, male, white, lowIncome, crit1, int1))

summary(data_access_inc_diaBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_access_inc_diaBP %>%
  select(-diaBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (income, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_access_inc_diaBP <- glmnet(x_hypos, data_access_inc_diaBP$diaBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_access_inc_diaBP

# Plot these results
plot(mod_access_inc_diaBP)

# Look at the variables included at each step
coef(mod_access_inc_diaBP, s = max(mod_access_inc_diaBP$lambda[mod_access_inc_diaBP$df == 4])); min(mod_access_inc_diaBP$dev.ratio[mod_access_inc_diaBP$df == 4])

coef(mod_access_inc_diaBP, s = max(mod_access_inc_diaBP$lambda[mod_access_inc_diaBP$df == 5])); min(mod_access_inc_diaBP$dev.ratio[mod_access_inc_diaBP$df == 5])

coef(mod_access_inc_diaBP, s = max(mod_access_inc_diaBP$lambda[mod_access_inc_diaBP$df == 6])); min(mod_access_inc_diaBP$dev.ratio[mod_access_inc_diaBP$df == 6])

coef(mod_access_inc_diaBP, s = max(mod_access_inc_diaBP$lambda[mod_access_inc_diaBP$df == 7])); min(mod_access_inc_diaBP$dev.ratio[mod_access_inc_diaBP$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_access_inc_diaBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_access_inc_diaBP$beta[, i] != 0))$names
  new_deviance <- mod_access_inc_diaBP$dev.ratio[i]
  new_varNum <- mod_access_inc_diaBP$df[i]
  new_lambda <- mod_access_inc_diaBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "crit1"] <- "crit1/green_inc23"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "lowIncome", ]
df <- df[df$Variables != "green_inc23", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_access_inc_diaBP$lambda, mod_access_inc_diaBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_inc_diaBP$lambda)), ylim = c(0, max(mod_access_inc_diaBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_access_inc_diaBP$log_lambda <- log(mod_access_inc_diaBP$lambda)
mod_access_inc_diaBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_access_inc_diaBP$log_lambda, mod_access_inc_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_inc_diaBP$log_lambda)), ylim = c(0.002, max(mod_access_inc_diaBP$dev.ratio)))
text(df$log_lambda, 0.002, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_accessIncomeDiaBP.pdf", height = 7, width = 11)
plot(mod_access_dep_diaBP$log_lambda, mod_access_dep_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_access_dep_diaBP$log_lambda)), ylim = c(0.002, max(mod_access_dep_diaBP$dev.ratio)))
text(df$log_lambda, 0.002, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (green_inc12_int) increases model fit of standard linear regression model
base_mod <- lm(data_access_inc_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "lowIncome"])
summary(base_mod)

param1_mod <- lm(data_access_inc_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "lowIncome"] + x_hypos[, "green_inc12_int"])
summary(param1_mod)

# Does seem to be an association here, with an increase in green space access between pregnancy and age 4 associated with lower diastolic BP among those with low incomes. While direction may be as predicted (increase in green space associated with lower BP), given that the effect size is tiny and interpretation not clear (why would reducing access to green space be associated with lower BP in previous models, but the opposite effect here?), this could just be the model hooking up to random noise.
anova(base_mod, param1_mod)

# Does adding the next parameter(s) (green_dec12) improve model fit? Yes, apparently. Now, an overall decrease in green space from pregnancy to age 4 associated with lower BP (this effect was found for the other diastolic BP results above, but doens't make a huge amount of biological sense.)
param2_mod <- lm(data_access_inc_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "lowIncome"] + x_hypos[, "green_inc12_int"] +
                   x_hypos[, "green_dec12"])
summary(param2_mod)

anova(param1_mod, param2_mod)

# Does adding the next parameter(s) (green_dec23) improve model fit? No.
param3_mod <- lm(data_access_inc_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "lowIncome"] + x_hypos[, "green_inc12_int"] +
                   x_hypos[, "green_dec12"] + x_hypos[, "green_dec23"])
summary(param3_mod)

anova(param2_mod, param3_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_access_inc_diaBP$diaBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - The 1SE model is just the null/covariate-only model, while the 'optimal'/lowest MSE model also contains green_dec12, green_inc12_int and green_dec12. The difference in MSE is minimal, though, so not especially convincing. Suggests at best a very weak association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between access to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and household income), although potentially decrease in access to green space between preg and age 4 linked to lower diastolic BP (but also increase in green space over this time if lower income...).




##############################################################################################################
##### Now, repeat all the above, but this time using 'distance to nearest green space' as the exposure

### Derive the different life-course variables and encode as hypotheses
summary(data$greenDist_preg)
summary(data$greenDist_4)
summary(data$greenDist_7)
table(data$edu, useNA = "ifany")

# Make a dataset just for this exposure
data_dist <- data


## Start with education as the SEP covariate/interaction term

# Critical period at first time point only
data_dist$crit1 <- data_dist$greenDist_preg

# Interaction between SEP and first time point
data_dist$int1 <- data_dist$edu * data_dist$greenDist_preg

# Critical period at second time point only
data_dist$crit2 <- data_dist$greenDist_4

# Interaction between SEP and second time point
data_dist$int2 <- data_dist$edu * data_dist$greenDist_4

# Critical period at third time point only
data_dist$crit3 <- data_dist$greenDist_7

# Interaction between SEP and third time point
data_dist$int3 <- data_dist$edu * data_dist$greenDist_7

# Linear accumulation of all exposures
data_dist$accumulation <- (data_dist$greenDist_preg + data_dist$greenDist_4 + data_dist$greenDist_7) / 3

# Interaction between SEP and cumulative exposure
data_dist$int_accum <- data_dist$edu * data_dist$accumulation

# Change in distance to green space between time 1 to time 2
data_dist$green_ch12 <- data_dist$greenDist_4 - data_dist$greenDist_preg

# Change in distance to green space between time 1 to time 2, with an interaction with SEP
data_dist$green_ch12_int <- (data_dist$greenDist_4 - data_dist$greenDist_preg) * data_dist$edu

# Change in distance to green space between time 2 to time 3
data_dist$green_ch23 <- data_dist$greenDist_7 - data_dist$greenDist_4

# Change in distance to green space between time 2 to time 3, with an interaction with SEP
data_dist$green_ch23_int <- (data_dist$greenDist_7 - data_dist$greenDist_4) * data_dist$edu


## Make a dataset just with the outcomes, exposure hypotheses and covariates
data_dist_edu <- data_dist %>%
  select(BMI_f7, overweight, sysBP, diaBP, 
         age_f7, male, white, edu, 
         crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_ch12, green_ch12_int, 
         green_ch23, green_ch23_int)


## Now analyse each outcome in turn, starting with BMI

# Reduce dataset down to just complete cases
data_dist_edu_bmi <- data_dist_edu %>%
  select(-overweight, -sysBP, -diaBP) %>%
  filter(complete.cases(BMI_f7, age_f7, male, white, edu, crit1, int1))

summary(data_dist_edu_bmi)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_edu_bmi %>%
  select(-BMI_f7)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_edu_bmi <- glmnet(x_hypos, data_dist_edu_bmi$BMI_f7, 
                             alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_edu_bmi

# Plot these results
plot(mod_dist_edu_bmi)

# Look at the variables included at each step
coef(mod_dist_edu_bmi, s = max(mod_dist_edu_bmi$lambda[mod_dist_edu_bmi$df == 4])); min(mod_dist_edu_bmi$dev.ratio[mod_dist_edu_bmi$df == 4])

coef(mod_dist_edu_bmi, s = max(mod_dist_edu_bmi$lambda[mod_dist_edu_bmi$df == 7])); min(mod_dist_edu_bmi$dev.ratio[mod_dist_edu_bmi$df == 7])

coef(mod_dist_edu_bmi, s = max(mod_dist_edu_bmi$lambda[mod_dist_edu_bmi$df == 8])); min(mod_dist_edu_bmi$dev.ratio[mod_dist_edu_bmi$df == 8])

coef(mod_dist_edu_bmi, s = max(mod_dist_edu_bmi$lambda[mod_dist_edu_bmi$df == 9])); min(mod_dist_edu_bmi$dev.ratio[mod_dist_edu_bmi$df == 9])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_edu_bmi$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_edu_bmi$beta[, i] != 0))$names
  new_deviance <- mod_dist_edu_bmi$dev.ratio[i]
  new_varNum <- mod_dist_edu_bmi$df[i]
  new_lambda <- mod_dist_edu_bmi$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "crit2"] <- "crit2/crit3/green_ch12"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "edu", ]
df <- df[df$Variables != "crit3", ]
df <- df[df$Variables != "green_ch12", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. This makes it clearer when different variables were added (not the neatest plot, though, as variables bunch up and overlap as lambda approaches 0...). Also, need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_edu_bmi$lambda, mod_dist_edu_bmi$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_bmi$lambda)), ylim = c(0, max(mod_dist_edu_bmi$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_edu_bmi$log_lambda <- log(mod_dist_edu_bmi$lambda)
mod_dist_edu_bmi

df$log_lambda <- log(as.numeric(df$Lambda))
df


# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_edu_bmi$log_lambda, mod_dist_edu_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_bmi$log_lambda)), ylim = c(0.007, max(mod_dist_edu_bmi$dev.ratio)))
text(df$log_lambda, 0.007, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceEduBMI.pdf", height = 7, width = 11)
plot(mod_dist_edu_bmi$log_lambda, mod_dist_edu_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_bmi$log_lambda)), ylim = c(0.007, max(mod_dist_edu_bmi$dev.ratio)))
text(df$log_lambda, 0.007, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(2) added (crit2, crit3 and green_ch12) increases model fit of standard linear regression model
base_mod <- lm(data_dist_edu_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "edu"])
summary(base_mod)

param1_mod <- lm(data_dist_edu_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                   x_hypos[, "edu"] + x_hypos[, "crit2"] + x_hypos[, "crit3"] + x_hypos[, "green_ch12"])
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between distance to green space in childhood and BMI at age 7.
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model
mod.cv <- cv.glmnet(x_hypos, data_dist_edu_bmi$BMI_f7, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at the top of the plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between distance to green space in childhood and BMI at age 7 (when controlling for ethnicity and parental education).



#### Next, will explore continuous distance to green space exposure, with education as covariate/interaction term, and binary overweight variable as the outcome

# Reduce dataset down to just complete cases
data_dist_edu_overweight <- data_dist_edu %>%
  select(-BMI_f7, -sysBP, -diaBP) %>%
  filter(complete.cases(overweight, age_f7, male, white, edu, crit1, int1))

summary(data_dist_edu_overweight)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_edu_overweight %>%
  select(-overweight)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_edu_over <- glmnet(x_hypos, data_dist_edu_overweight$overweight, family = "binomial",
                              alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_edu_over

# Plot these results
plot(mod_dist_edu_over)

# Look at the variables included at each step
coef(mod_dist_edu_over, s = max(mod_dist_edu_over$lambda[mod_dist_edu_over$df == 4])); min(mod_dist_edu_over$dev.ratio[mod_dist_edu_over$df == 4])

coef(mod_dist_edu_over, s = max(mod_dist_edu_over$lambda[mod_dist_edu_over$df == 5])); min(mod_dist_edu_over$dev.ratio[mod_dist_edu_over$df == 5])

coef(mod_dist_edu_over, s = max(mod_dist_edu_over$lambda[mod_dist_edu_over$df == 6])); min(mod_dist_edu_over$dev.ratio[mod_dist_edu_over$df == 6])

coef(mod_dist_edu_over, s = max(mod_dist_edu_over$lambda[mod_dist_edu_over$df == 7])); min(mod_dist_edu_over$dev.ratio[mod_dist_edu_over$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_edu_over$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_edu_over$beta[, i] != 0))$names
  new_deviance <- mod_dist_edu_over$dev.ratio[i]
  new_varNum <- mod_dist_edu_over$df[i]
  new_lambda <- mod_dist_edu_over$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "edu", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_edu_over$lambda, mod_dist_edu_over$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_over$lambda)), ylim = c(0, max(mod_dist_edu_over$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_edu_over$log_lambda <- log(mod_dist_edu_over$lambda)
mod_dist_edu_over

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_edu_over$log_lambda, mod_dist_edu_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_over$log_lambda)), ylim = c(0.001, max(mod_dist_edu_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceEduOverweight.pdf", height = 7, width = 11)
plot(mod_dist_edu_over$log_lambda, mod_dist_edu_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_over$log_lambda)), ylim = c(0.001, max(mod_dist_edu_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (crit2) increases model fit of standard logistic regression model
base_mod <- glm(data_dist_edu_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                  x_hypos[, "edu"], family = "binomial")
summary(base_mod)

param1_mod <- glm(data_dist_edu_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                    x_hypos[, "white"] + x_hypos[, "edu"] + x_hypos[, "crit2"], family = "binomial")
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between distance to green space in childhood and being overweight at age 7.
anova(base_mod, param1_mod, test = "Chisq")


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_dist_edu_overweight$overweight, family = "binomial", 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between distance to green space in childhood and being overweight at age 7 (when controlling for ethnicity and parental education).



#### Next, will explore continuous distance to green space exposure, with education as covariate/interaction term, and continuous systolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_dist_edu_sysBP <- data_dist_edu %>%
  select(-BMI_f7, -overweight, -diaBP) %>%
  filter(complete.cases(sysBP, age_f7, male, white, edu, crit1, int1))

summary(data_dist_edu_sysBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_edu_sysBP %>%
  select(-sysBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_edu_sysBP <- glmnet(x_hypos, data_dist_edu_sysBP$sysBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_edu_sysBP

# Plot these results
plot(mod_dist_edu_sysBP)

# Look at the variables included at each step
coef(mod_dist_edu_sysBP, s = max(mod_dist_edu_sysBP$lambda[mod_dist_edu_sysBP$df == 4])); min(mod_dist_edu_sysBP$dev.ratio[mod_dist_edu_sysBP$df == 4])

coef(mod_dist_edu_sysBP, s = max(mod_dist_edu_sysBP$lambda[mod_dist_edu_sysBP$df == 5])); min(mod_dist_edu_sysBP$dev.ratio[mod_dist_edu_sysBP$df == 5])

coef(mod_dist_edu_sysBP, s = max(mod_dist_edu_sysBP$lambda[mod_dist_edu_sysBP$df == 6])); min(mod_dist_edu_sysBP$dev.ratio[mod_dist_edu_sysBP$df == 6])

coef(mod_dist_edu_sysBP, s = max(mod_dist_edu_sysBP$lambda[mod_dist_edu_sysBP$df == 7])); min(mod_dist_edu_sysBP$dev.ratio[mod_dist_edu_sysBP$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_edu_sysBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_edu_sysBP$beta[, i] != 0))$names
  new_deviance <- mod_dist_edu_sysBP$dev.ratio[i]
  new_varNum <- mod_dist_edu_sysBP$df[i]
  new_lambda <- mod_dist_edu_sysBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "edu", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_edu_sysBP$lambda, mod_dist_edu_sysBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_sysBP$lambda)), ylim = c(0, max(mod_dist_edu_sysBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_edu_sysBP$log_lambda <- log(mod_dist_edu_sysBP$lambda)
mod_dist_edu_sysBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_edu_sysBP$log_lambda, mod_dist_edu_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_sysBP$log_lambda)), ylim = c(0.005, max(mod_dist_edu_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distEduSysBP.pdf", height = 7, width = 11)
plot(mod_dist_edu_sysBP$log_lambda, mod_dist_edu_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_sysBP$log_lambda)), ylim = c(0.005, max(mod_dist_edu_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (accumulation) increases model fit of standard linear regression model
base_mod <- lm(data_dist_edu_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "edu"])
summary(base_mod)

param1_mod <- lm(data_dist_edu_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "edu"] + x_hypos[, "accumulation"])
summary(param1_mod)

# Is potentially a weak effect of accumulation, with higher average distance associated with ever-so-marginally lower systolic blood pressure (but effect is so marginal, and in opposite direction to what would be expected, may just be latching on to random noise in the data).
anova(base_mod, param1_mod)

# What about adding next variable (int3)? No.
param2_mod <- lm(data_dist_edu_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "edu"] + x_hypos[, "accumulation"] + x_hypos[, "int3"])
summary(param2_mod)

anova(param1_mod, param2_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_dist_edu_sysBP$sysBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with distance to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between distance to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and parental education), although perhaps very weak association with accumulation/greater average distance and lower systolic BP.



#### Next, will explore continuous distance to green space exposure, with education as covariate/interaction term, and continuous diastolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_dist_edu_diaBP <- data_dist_edu %>%
  select(-BMI_f7, -overweight, -sysBP) %>%
  filter(complete.cases(diaBP, age_f7, male, white, edu, crit1, int1))

summary(data_dist_edu_diaBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_edu_diaBP %>%
  select(-diaBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_edu_diaBP <- glmnet(x_hypos, data_dist_edu_diaBP$diaBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_edu_diaBP

# Plot these results
plot(mod_dist_edu_diaBP)

# Look at the variables included at each step
coef(mod_dist_edu_diaBP, s = max(mod_dist_edu_diaBP$lambda[mod_dist_edu_diaBP$df == 4])); min(mod_dist_edu_diaBP$dev.ratio[mod_dist_edu_diaBP$df == 4])

coef(mod_dist_edu_diaBP, s = max(mod_dist_edu_diaBP$lambda[mod_dist_edu_diaBP$df == 5])); min(mod_dist_edu_diaBP$dev.ratio[mod_dist_edu_diaBP$df == 5])

coef(mod_dist_edu_diaBP, s = max(mod_dist_edu_diaBP$lambda[mod_dist_edu_diaBP$df == 6])); min(mod_dist_edu_diaBP$dev.ratio[mod_dist_edu_diaBP$df == 6])

coef(mod_dist_edu_diaBP, s = max(mod_dist_edu_diaBP$lambda[mod_dist_edu_diaBP$df == 7])); min(mod_dist_edu_diaBP$dev.ratio[mod_dist_edu_diaBP$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_edu_diaBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_edu_diaBP$beta[, i] != 0))$names
  new_deviance <- mod_dist_edu_diaBP$dev.ratio[i]
  new_varNum <- mod_dist_edu_diaBP$df[i]
  new_lambda <- mod_dist_edu_diaBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "edu", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_edu_diaBP$lambda, mod_dist_edu_diaBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_diaBP$lambda)), ylim = c(0, max(mod_dist_edu_diaBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_edu_diaBP$log_lambda <- log(mod_dist_edu_diaBP$lambda)
mod_dist_edu_diaBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome). No obvious strong effects
plot(mod_dist_edu_diaBP$log_lambda, mod_dist_edu_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_diaBP$log_lambda)), ylim = c(0.005, max(mod_dist_edu_diaBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distEduDiaBP.pdf", height = 7, width = 11)
plot(mod_dist_edu_diaBP$log_lambda, mod_dist_edu_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_edu_diaBP$log_lambda)), ylim = c(0.005, max(mod_dist_edu_diaBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (green_ch12) increases model fit of standard linear regression model
base_mod <- lm(data_dist_edu_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "edu"])
summary(base_mod)

param1_mod <- lm(data_dist_edu_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "edu"] + x_hypos[, "green_ch12"])
summary(param1_mod)

# No association here
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_dist_edu_diaBP$diaBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between distance to green space in childhood and diastolic blood pressure at age 7 (when controlling for ethnicity and parental education).



##############################################################################################################
#### Next, explore whether continuous distance to green space as exposure associated with cardiometabolic outcomes with deprivation as covariate/interaction term (instead of parental education)

### Derive the different life-course variables and encode as hypotheses
summary(data$greenDist_preg)
summary(data$greenDist_4)
summary(data$greenDist_7)
table(data$deprived, useNA = "ifany")

# Make a dataset just for this exposure
data_dist <- data


## Here deprivation is the SEP covariate/interaction term

# Critical period at first time point only
data_dist$crit1 <- data_dist$greenDist_preg

# Interaction between SEP and first time point
data_dist$int1 <- data_dist$deprived * data_dist$greenDist_preg

# Critical period at second time point only
data_dist$crit2 <- data_dist$greenDist_4

# Interaction between SEP and second time point
data_dist$int2 <- data_dist$deprived * data_dist$greenDist_4

# Critical period at third time point only
data_dist$crit3 <- data_dist$greenDist_7

# Interaction between SEP and third time point
data_dist$int3 <- data_dist$deprived * data_dist$greenDist_7

# Linear accumulation of all exposures
data_dist$accumulation <- (data_dist$greenDist_preg + data_dist$greenDist_4 + data_dist$greenDist_7) / 3

# Interaction between SEP and cumulative exposure
data_dist$int_accum <- data_dist$deprived * data_dist$accumulation

# Change in distance to green space between time 1 to time 2
data_dist$green_ch12 <- data_dist$greenDist_4 - data_dist$greenDist_preg

# Change in distance to green space between time 1 to time 2, with an interaction with SEP
data_dist$green_ch12_int <- (data_dist$greenDist_4 - data_dist$greenDist_preg) * data_dist$deprived

# Change in distance to green space between time 2 to time 3
data_dist$green_ch23 <- data_dist$greenDist_7 - data_dist$greenDist_4

# Change in distance to green space between time 2 to time 3, with an interaction with SEP
data_dist$green_ch23_int <- (data_dist$greenDist_7 - data_dist$greenDist_4) * data_dist$deprived


## Make a dataset just with the outcomes, exposure hypotheses and covariates
data_dist_dep <- data_dist %>%
  select(BMI_f7, overweight, sysBP, diaBP, 
         age_f7, male, white, deprived, 
         crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_ch12, green_ch12_int, 
         green_ch23, green_ch23_int)


## Now analyse each outcome in turn, starting with BMI

# Reduce dataset down to just complete cases
data_dist_dep_bmi <- data_dist_dep %>%
  select(-overweight, -sysBP, -diaBP) %>%
  filter(complete.cases(BMI_f7, age_f7, male, white, deprived, crit1, int1))

summary(data_dist_dep_bmi)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_dep_bmi %>%
  select(-BMI_f7)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (dep, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_dep_bmi <- glmnet(x_hypos, data_dist_dep_bmi$BMI_f7, 
                             alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_dep_bmi

# Plot these results
plot(mod_dist_dep_bmi)

# Look at the variables included at each step
coef(mod_dist_dep_bmi, s = max(mod_dist_dep_bmi$lambda[mod_dist_dep_bmi$df == 4])); min(mod_dist_dep_bmi$dev.ratio[mod_dist_dep_bmi$df == 4])

coef(mod_dist_dep_bmi, s = max(mod_dist_dep_bmi$lambda[mod_dist_dep_bmi$df == 5])); min(mod_dist_dep_bmi$dev.ratio[mod_dist_dep_bmi$df == 5])

coef(mod_dist_dep_bmi, s = max(mod_dist_dep_bmi$lambda[mod_dist_dep_bmi$df == 6])); min(mod_dist_dep_bmi$dev.ratio[mod_dist_dep_bmi$df == 6])

coef(mod_dist_dep_bmi, s = max(mod_dist_dep_bmi$lambda[mod_dist_dep_bmi$df == 7])); min(mod_dist_dep_bmi$dev.ratio[mod_dist_dep_bmi$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_dep_bmi$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_dep_bmi$beta[, i] != 0))$names
  new_deviance <- mod_dist_dep_bmi$dev.ratio[i]
  new_varNum <- mod_dist_dep_bmi$df[i]
  new_lambda <- mod_dist_dep_bmi$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "int2"] <- "int2/green_ch23_int"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "deprived", ]
df <- df[df$Variables != "green_ch23_int", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. This makes it clearer when different variables were added (not the neatest plot, though, as variables bunch up and overlap as lambda approaches 0...). Also, need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_dep_bmi$lambda, mod_dist_dep_bmi$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_bmi$lambda)), ylim = c(0, max(mod_dist_dep_bmi$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_dep_bmi$log_lambda <- log(mod_dist_dep_bmi$lambda)
mod_dist_dep_bmi

df$log_lambda <- log(as.numeric(df$Lambda))
df


# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_dep_bmi$log_lambda, mod_dist_dep_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_bmi$log_lambda)), ylim = c(0.008, max(mod_dist_dep_bmi$dev.ratio)))
text(df$log_lambda, 0.008, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceDepBMI.pdf", height = 7, width = 11)
plot(mod_dist_dep_bmi$log_lambda, mod_dist_dep_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_bmi$log_lambda)), ylim = c(0.008, max(mod_dist_dep_bmi$dev.ratio)))
text(df$log_lambda, 0.008, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter added (crit3) increases model fit of standard linear regression model
base_mod <- lm(data_dist_dep_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "deprived"])
summary(base_mod)

param1_mod <- lm(data_dist_dep_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                   x_hypos[, "deprived"] + x_hypos[, "crit3"])
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between distance to green space in childhood and BMI at age 7.
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model
mod.cv <- cv.glmnet(x_hypos, data_dist_dep_bmi$BMI_f7, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at the top of the plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with distance to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between distance to green space in childhood and BMI at age 7 (when controlling for ethnicity and household deprivation).



#### Next, will explore continuous distance to green space exposure, with deprivation as covariate/interaction term, and binary overweight variable as the outcome

# Reduce dataset down to just complete cases
data_dist_dep_overweight <- data_dist_dep %>%
  select(-BMI_f7, -sysBP, -diaBP) %>%
  filter(complete.cases(overweight, age_f7, male, white, deprived, crit1, int1))

summary(data_dist_dep_overweight)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_dep_overweight %>%
  select(-overweight)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_dep_over <- glmnet(x_hypos, data_dist_dep_overweight$overweight, family = "binomial",
                              alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_dep_over

# Plot these results
plot(mod_dist_dep_over)

# Look at the variables included at each step
coef(mod_dist_dep_over, s = max(mod_dist_dep_over$lambda[mod_dist_dep_over$df == 4])); min(mod_dist_dep_over$dev.ratio[mod_dist_dep_over$df == 4])

coef(mod_dist_dep_over, s = max(mod_dist_dep_over$lambda[mod_dist_dep_over$df == 5])); min(mod_dist_dep_over$dev.ratio[mod_dist_dep_over$df == 5])

coef(mod_dist_dep_over, s = max(mod_dist_dep_over$lambda[mod_dist_dep_over$df == 6])); min(mod_dist_dep_over$dev.ratio[mod_dist_dep_over$df == 6])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_dep_over$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_dep_over$beta[, i] != 0))$names
  new_deviance <- mod_dist_dep_over$dev.ratio[i]
  new_varNum <- mod_dist_dep_over$df[i]
  new_lambda <- mod_dist_dep_over$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "int1"] <- "int1/crit2"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "deprived", ]
df <- df[df$Variables != "crit2", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_dep_over$lambda, mod_dist_dep_over$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_over$lambda)), ylim = c(0, max(mod_dist_dep_over$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_dep_over$log_lambda <- log(mod_dist_dep_over$lambda)
mod_dist_dep_over

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_dep_over$log_lambda, mod_dist_dep_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_over$log_lambda)), ylim = c(0.001, max(mod_dist_dep_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceDepOverweight.pdf", height = 7, width = 11)
plot(mod_dist_dep_over$log_lambda, mod_dist_dep_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_over$log_lambda)), ylim = c(0.001, max(mod_dist_dep_over$dev.ratio)))
text(df$log_lambda, 0.001, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (crit3) increases model fit of standard logistic regression model
base_mod <- glm(data_dist_dep_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                  x_hypos[, "deprived"], family = "binomial")
summary(base_mod)

param1_mod <- glm(data_dist_dep_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                    x_hypos[, "white"] + x_hypos[, "deprived"] + x_hypos[, "crit3"],
                  family = "binomial")
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between access to green space in childhood and being overweight at age 7.
anova(base_mod, param1_mod, test = "Chisq")


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_dist_dep_overweight$overweight, family = "binomial", 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between distance to green space in childhood and being overweight at age 7 (when controlling for ethnicity and household deprivation).



#### Next, will explore continuous distance to green space exposure, with deprivation as covariate/interaction term, and continuous systolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_dist_dep_sysBP <- data_dist_dep %>%
  select(-BMI_f7, -overweight, -diaBP) %>%
  filter(complete.cases(sysBP, age_f7, male, white, deprived, crit1, int1))

summary(data_dist_dep_sysBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_dep_sysBP %>%
  select(-sysBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_dep_sysBP <- glmnet(x_hypos, data_dist_dep_sysBP$sysBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_dep_sysBP

# Plot these results
plot(mod_dist_dep_sysBP)

# Look at the variables included at each step
coef(mod_dist_dep_sysBP, s = max(mod_dist_dep_sysBP$lambda[mod_dist_dep_sysBP$df == 4])); min(mod_dist_dep_sysBP$dev.ratio[mod_dist_dep_sysBP$df == 4])

coef(mod_dist_dep_sysBP, s = max(mod_dist_dep_sysBP$lambda[mod_dist_dep_sysBP$df == 5])); min(mod_dist_dep_sysBP$dev.ratio[mod_dist_dep_sysBP$df == 5])

coef(mod_dist_dep_sysBP, s = max(mod_dist_dep_sysBP$lambda[mod_dist_dep_sysBP$df == 6])); min(mod_dist_dep_sysBP$dev.ratio[mod_dist_dep_sysBP$df == 6])

coef(mod_dist_dep_sysBP, s = max(mod_dist_dep_sysBP$lambda[mod_dist_dep_sysBP$df == 7])); min(mod_dist_dep_sysBP$dev.ratio[mod_dist_dep_sysBP$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_dep_sysBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_dep_sysBP$beta[, i] != 0))$names
  new_deviance <- mod_dist_dep_sysBP$dev.ratio[i]
  new_varNum <- mod_dist_dep_sysBP$df[i]
  new_lambda <- mod_dist_dep_sysBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "deprived", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_dep_sysBP$lambda, mod_dist_dep_sysBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_sysBP$lambda)), ylim = c(0, max(mod_dist_dep_sysBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_dep_sysBP$log_lambda <- log(mod_dist_dep_sysBP$lambda)
mod_dist_dep_sysBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_dep_sysBP$log_lambda, mod_dist_dep_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_sysBP$log_lambda)), ylim = c(0.005, max(mod_dist_dep_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceDepSysBP.pdf", height = 7, width = 11)
plot(mod_dist_dep_sysBP$log_lambda, mod_dist_dep_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_sysBP$log_lambda)), ylim = c(0.005, max(mod_dist_dep_sysBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that not much is really going on here, although first variable eneterd ('accumulation') does appear to increase model fit marginally... - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (green_dec12) increases model fit of standard linear regression model
base_mod <- lm(data_dist_dep_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "deprived"])
summary(base_mod)

param1_mod <- lm(data_dist_dep_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "deprived"] + x_hypos[, "accumulation"])
summary(param1_mod)

# Does seem to be a marginal association here, with an increase in average distance from space over all time-points associated with lower systolic BP. However, given that the effect size is tiny and interpretation not clear (why would increasing distance to green space be associated with lower BP?), this could just be the model hooking up to random noise.
anova(base_mod, param1_mod)

# Does adding the next parameter (crit2) improve model fit? No.
param2_mod <- lm(data_dist_dep_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "deprived"] + x_hypos[, "accumulation"] + x_hypos[, "crit2"])
summary(param2_mod)

anova(param1_mod, param2_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_dist_dep_sysBP$sysBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - The 1SE model is just the null/covariate-only model, while the 'optimal'/lowest MSE model also contains green_dec12_int. The difference in MSE is minimal, though, so not especially convincing. Suggests at best a very weak association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between distance to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and household deprivation), although potentially increase in average distance to green space (accumulation) and increased distanceat age 4 linked to lower systolic BP.



#### Next, will explore continuous distance to green space exposure, with deprivation as covariate/interaction term, and continuous diastolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_dist_dep_diaBP <- data_dist_dep %>%
  select(-BMI_f7, -overweight, -sysBP) %>%
  filter(complete.cases(diaBP, age_f7, male, white, deprived, crit1, int1))

summary(data_dist_dep_diaBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_dep_diaBP %>%
  select(-diaBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_dep_diaBP <- glmnet(x_hypos, data_dist_dep_diaBP$diaBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_dep_diaBP

# Plot these results
plot(mod_dist_dep_diaBP)

# Look at the variables included at each step
coef(mod_dist_dep_diaBP, s = max(mod_dist_dep_diaBP$lambda[mod_dist_dep_diaBP$df == 4])); min(mod_dist_dep_diaBP$dev.ratio[mod_dist_dep_diaBP$df == 4])

coef(mod_dist_dep_diaBP, s = max(mod_dist_dep_diaBP$lambda[mod_dist_dep_diaBP$df == 5])); min(mod_dist_dep_diaBP$dev.ratio[mod_dist_dep_diaBP$df == 5])

coef(mod_dist_dep_diaBP, s = max(mod_dist_dep_diaBP$lambda[mod_dist_dep_diaBP$df == 6])); min(mod_dist_dep_diaBP$dev.ratio[mod_dist_dep_diaBP$df == 6])

coef(mod_dist_dep_diaBP, s = max(mod_dist_dep_diaBP$lambda[mod_dist_dep_diaBP$df == 7])); min(mod_dist_dep_diaBP$dev.ratio[mod_dist_dep_diaBP$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_dep_diaBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_dep_diaBP$beta[, i] != 0))$names
  new_deviance <- mod_dist_dep_diaBP$dev.ratio[i]
  new_varNum <- mod_dist_dep_diaBP$df[i]
  new_lambda <- mod_dist_dep_diaBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "deprived", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_dep_diaBP$lambda, mod_dist_dep_diaBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_diaBP$lambda)), ylim = c(0, max(mod_dist_dep_diaBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_dep_diaBP$log_lambda <- log(mod_dist_dep_diaBP$lambda)
mod_dist_dep_diaBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_dep_diaBP$log_lambda, mod_dist_dep_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_diaBP$log_lambda)), ylim = c(0.005, max(mod_dist_dep_diaBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceDepDiaBP.pdf", height = 7, width = 11)
plot(mod_dist_dep_diaBP$log_lambda, mod_dist_dep_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_diaBP$log_lambda)), ylim = c(0.005, max(mod_dist_dep_diaBP$dev.ratio)))
text(df$log_lambda, 0.005, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (green_ch12) increases model fit of standard linear regression model
base_mod <- lm(data_dist_dep_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "deprived"])
summary(base_mod)

param1_mod <- lm(data_dist_dep_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "deprived"] + x_hypos[, "green_ch12"])
summary(param1_mod)

# No association here
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_dist_dep_diaBP$diaBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with distance to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between distance to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and household deprivation), although potentially increase in average distance to green space (accumulation) linked to lower systolic BP.




##############################################################################################################
#### Next, explore whether continuous distance to green space as exposure associated with cardiometabolic outcomes with disposable household income as covariate/interaction term (instead of parental education)

### Derive the different life-course variables and encode as hypotheses
summary(data$greenDist_preg)
summary(data$greenDist_4)
summary(data$greenDist_7)
table(data$lowIncome, useNA = "ifany")

# Make a dataset just for this exposure
data_dist <- data


## Here low household income is the SEP covariate/interaction term

# Critical period at first time point only
data_dist$crit1 <- data_dist$greenDist_preg

# Interaction between SEP and first time point
data_dist$int1 <- data_dist$lowIncome * data_dist$greenDist_preg

# Critical period at second time point only
data_dist$crit2 <- data_dist$greenDist_4

# Interaction between SEP and second time point
data_dist$int2 <- data_dist$lowIncome * data_dist$greenDist_4

# Critical period at third time point only
data_dist$crit3 <- data_dist$greenDist_7

# Interaction between SEP and third time point
data_dist$int3 <- data_dist$lowIncome * data_dist$greenDist_7

# Linear accumulation of all exposures
data_dist$accumulation <- (data_dist$greenDist_preg + data_dist$greenDist_4 + data_dist$greenDist_7) / 3

# Interaction between SEP and cumulative exposure
data_dist$int_accum <- data_dist$lowIncome * data_dist$accumulation

# Change in distance to green space between time 1 to time 2
data_dist$green_ch12 <- data_dist$greenDist_4 - data_dist$greenDist_preg

# Change in distance to green space between time 1 to time 2, with an interaction with SEP
data_dist$green_ch12_int <- (data_dist$greenDist_4 - data_dist$greenDist_preg) * data_dist$lowIncome

# Change in distance to green space between time 2 to time 3
data_dist$green_ch23 <- data_dist$greenDist_7 - data_dist$greenDist_4

# Change in distance to green space between time 2 to time 3, with an interaction with SEP
data_dist$green_ch23_int <- (data_dist$greenDist_7 - data_dist$greenDist_4) * data_dist$lowIncome


## Make a dataset just with the outcomes, exposure hypotheses and covariates
data_dist_inc <- data_dist %>%
  select(BMI_f7, overweight, sysBP, diaBP, 
         age_f7, male, white, lowIncome, 
         crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_ch12, green_ch12_int, 
         green_ch23, green_ch23_int)


## Now analyse each outcome in turn, starting with BMI

# Reduce dataset down to just complete cases
data_dist_inc_bmi <- data_dist_inc %>%
  select(-overweight, -sysBP, -diaBP) %>%
  filter(complete.cases(BMI_f7, age_f7, male, white, lowIncome, crit1, int1))

summary(data_dist_inc_bmi)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_inc_bmi %>%
  select(-BMI_f7)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (income, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_inc_bmi <- glmnet(x_hypos, data_dist_inc_bmi$BMI_f7, 
                             alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_inc_bmi

# Plot these results
plot(mod_dist_inc_bmi)

# Look at the variables included at each step
coef(mod_dist_inc_bmi, s = max(mod_dist_inc_bmi$lambda[mod_dist_inc_bmi$df == 4])); min(mod_dist_inc_bmi$dev.ratio[mod_dist_inc_bmi$df == 4])

coef(mod_dist_inc_bmi, s = max(mod_dist_inc_bmi$lambda[mod_dist_inc_bmi$df == 5])); min(mod_dist_inc_bmi$dev.ratio[mod_dist_inc_bmi$df == 5])

coef(mod_dist_inc_bmi, s = max(mod_dist_inc_bmi$lambda[mod_dist_inc_bmi$df == 6])); min(mod_dist_inc_bmi$dev.ratio[mod_dist_inc_bmi$df == 6])

coef(mod_dist_inc_bmi, s = max(mod_dist_inc_bmi$lambda[mod_dist_inc_bmi$df == 7])); min(mod_dist_inc_bmi$dev.ratio[mod_dist_inc_bmi$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_inc_bmi$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_inc_bmi$beta[, i] != 0))$names
  new_deviance <- mod_dist_inc_bmi$dev.ratio[i]
  new_varNum <- mod_dist_inc_bmi$df[i]
  new_lambda <- mod_dist_inc_bmi$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "lowIncome", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. This makes it clearer when different variables were added (not the neatest plot, though, as variables bunch up and overlap as lambda approaches 0...). Also, need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_inc_bmi$lambda, mod_dist_inc_bmi$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_inc_bmi$lambda)), ylim = c(0, max(mod_dist_inc_bmi$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)

## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_inc_bmi$log_lambda <- log(mod_dist_inc_bmi$lambda)
mod_dist_inc_bmi

df$log_lambda <- log(as.numeric(df$Lambda))
df


# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_inc_bmi$log_lambda, mod_dist_inc_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_inc_bmi$log_lambda)), ylim = c(0.006, max(mod_dist_inc_bmi$dev.ratio)))
text(df$log_lambda, 0.006, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceIncomeBMI.pdf", height = 7, width = 11)
plot(mod_dist_inc_bmi$log_lambda, mod_dist_inc_bmi$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_inc_bmi$log_lambda)), ylim = c(0.006, max(mod_dist_inc_bmi$dev.ratio)))
text(df$log_lambda, 0.006, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter added (crit2) increases model fit of standard linear regression model
base_mod <- lm(data_dist_inc_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "lowIncome"])
summary(base_mod)

param1_mod <- lm(data_dist_inc_bmi$BMI_f7 ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                   x_hypos[, "lowIncome"] + x_hypos[, "crit2"])
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between distance to green space in childhood and BMI at age 7.
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model
mod.cv <- cv.glmnet(x_hypos, data_dist_inc_bmi$BMI_f7, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at the top of the plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with distance to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between distance to green space in childhood and BMI at age 7 (when controlling for ethnicity and household income).



#### Next, will explore continuous distance to green space exposure, with income as covariate/interaction term, and binary overweight variable as the outcome

# Reduce dataset down to just complete cases
data_dist_inc_overweight <- data_dist_inc %>%
  select(-BMI_f7, -sysBP, -diaBP) %>%
  filter(complete.cases(overweight, age_f7, male, white, lowIncome, crit1, int1))

summary(data_dist_inc_overweight)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_inc_overweight %>%
  select(-overweight)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (edu, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_inc_over <- glmnet(x_hypos, data_dist_inc_overweight$overweight, family = "binomial",
                              alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_inc_over

# Plot these results
plot(mod_dist_inc_over)

# Look at the variables included at each step
coef(mod_dist_inc_over, s = max(mod_dist_inc_over$lambda[mod_dist_inc_over$df == 4])); min(mod_dist_inc_over$dev.ratio[mod_dist_inc_over$df == 4])

coef(mod_dist_inc_over, s = max(mod_dist_inc_over$lambda[mod_dist_inc_over$df == 5])); min(mod_dist_inc_over$dev.ratio[mod_dist_inc_over$df == 5])

coef(mod_dist_inc_over, s = max(mod_dist_inc_over$lambda[mod_dist_inc_over$df == 6])); min(mod_dist_inc_over$dev.ratio[mod_dist_inc_over$df == 6])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_inc_over$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_inc_over$beta[, i] != 0))$names
  new_deviance <- mod_dist_inc_over$dev.ratio[i]
  new_varNum <- mod_dist_inc_over$df[i]
  new_lambda <- mod_dist_inc_over$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "lowIncome", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_inc_over$lambda, mod_dist_inc_over$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_inc_over$lambda)), ylim = c(0, max(mod_dist_inc_over$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_inc_over$log_lambda <- log(mod_dist_inc_over$lambda)
mod_dist_inc_over

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome). (Have to add a bit of faff here to get y-axis to not display in scientific notation)
plot(mod_dist_inc_over$log_lambda, mod_dist_inc_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", yaxt = "n",
     xlim = rev(range(mod_dist_inc_over$log_lambda)), ylim = c(0.000, max(mod_dist_inc_over$dev.ratio)))
axis(2, at = c(0, 0.0001, 0.0002, 0.0003, 0.0004), labels = format(c(0, 0.0001, 0.0002, 0.0003, 0.0004), scientific = FALSE))
text(df$log_lambda, 0.000, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceIncomeOverweight.pdf", height = 7, width = 11)
plot(mod_dist_inc_over$log_lambda, mod_dist_inc_over$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", yaxt = "n",
     xlim = rev(range(mod_dist_inc_over$log_lambda)), ylim = c(0.000, max(mod_dist_inc_over$dev.ratio)))
axis(2, at = c(0, 0.0001, 0.0002, 0.0003, 0.0004), labels = format(c(0, 0.0001, 0.0002, 0.0003, 0.0004), scientific = FALSE))
text(df$log_lambda, 0.000, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (crit2) increases model fit of standard logistic regression model
base_mod <- glm(data_dist_inc_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                  x_hypos[, "lowIncome"], family = "binomial")
summary(base_mod)

param1_mod <- glm(data_dist_inc_overweight$overweight ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                    x_hypos[, "white"] + x_hypos[, "lowIncome"] + x_hypos[, "crit1"],
                  family = "binomial")
summary(param1_mod)

# Nope, is pretty much no association here, suggesting no association between access to green space in childhood and being overweight at age 7.
anova(base_mod, param1_mod, test = "Chisq")


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_dist_inc_overweight$overweight, family = "binomial", 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up pretty well, and indicate very little association between distance to green space in childhood and being overweight at age 7 (when controlling for ethnicity and household income).



#### Next, will explore continuous distance to green space exposure, with income as covariate/interaction term, and continuous systolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_dist_inc_sysBP <- data_dist_inc %>%
  select(-BMI_f7, -overweight, -diaBP) %>%
  filter(complete.cases(sysBP, age_f7, male, white, lowIncome, crit1, int1))

summary(data_dist_inc_sysBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_inc_sysBP %>%
  select(-sysBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (income, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_inc_sysBP <- glmnet(x_hypos, data_dist_inc_sysBP$sysBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_inc_sysBP

# Plot these results
plot(mod_dist_inc_sysBP)

# Look at the variables included at each step
coef(mod_dist_inc_sysBP, s = max(mod_dist_inc_sysBP$lambda[mod_dist_inc_sysBP$df == 4])); min(mod_dist_inc_sysBP$dev.ratio[mod_dist_inc_sysBP$df == 4])

coef(mod_dist_inc_sysBP, s = max(mod_dist_inc_sysBP$lambda[mod_dist_inc_sysBP$df == 6])); min(mod_dist_inc_sysBP$dev.ratio[mod_dist_inc_sysBP$df == 6])

coef(mod_dist_inc_sysBP, s = max(mod_dist_inc_sysBP$lambda[mod_dist_inc_sysBP$df == 7])); min(mod_dist_inc_sysBP$dev.ratio[mod_dist_inc_sysBP$df == 7])

coef(mod_dist_inc_sysBP, s = max(mod_dist_inc_sysBP$lambda[mod_dist_inc_sysBP$df == 8])); min(mod_dist_inc_sysBP$dev.ratio[mod_dist_inc_sysBP$df == 8])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_inc_sysBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_inc_sysBP$beta[, i] != 0))$names
  new_deviance <- mod_dist_inc_sysBP$dev.ratio[i]
  new_varNum <- mod_dist_inc_sysBP$df[i]
  new_lambda <- mod_dist_inc_sysBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df$Variables[df$Variables == "crit2"] <- "crit2/accumulation"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "lowIncome", ]
df <- df[df$Variables != "accumulation", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_inc_sysBP$lambda, mod_dist_inc_sysBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_inc_sysBP$lambda)), ylim = c(0, max(mod_dist_inc_sysBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_inc_sysBP$log_lambda <- log(mod_dist_inc_sysBP$lambda)
mod_dist_inc_sysBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_inc_sysBP$log_lambda, mod_dist_inc_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_inc_sysBP$log_lambda)), ylim = c(0.004, max(mod_dist_inc_sysBP$dev.ratio)))
text(df$log_lambda, 0.004, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceIncSysBP.pdf", height = 7, width = 11)
plot(mod_dist_inc_sysBP$log_lambda, mod_dist_inc_sysBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_inc_sysBP$log_lambda)), ylim = c(0.004, max(mod_dist_inc_sysBP$dev.ratio)))
text(df$log_lambda, 0.004, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that not much is really going on here (although potentially v. weak effect of crit2/accumulation) - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (crit2 and accumulation) increases model fit of standard linear regression model
base_mod <- lm(data_dist_inc_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "lowIncome"])
summary(base_mod)

param1_mod <- lm(data_dist_inc_sysBP$sysBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "lowIncome"] + x_hypos[, "crit2"] + 
                   x_hypos[, "accumulation"])
summary(param1_mod)

# No real association here
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_dist_inc_sysBP$sysBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv)

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with distance to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between distance to green space in childhood and systolic blood pressure at age 7 (when controlling for ethnicity and household income).



#### Next, will explore continuous distance to green space exposure, with income as covariate/interaction term, and continuous diastolic BP variable as the outcome

# Reduce dataset down to just complete cases
data_dist_inc_diaBP <- data_dist_inc %>%
  select(-BMI_f7, -overweight, -sysBP) %>%
  filter(complete.cases(diaBP, age_f7, male, white, lowIncome, crit1, int1))

summary(data_dist_inc_diaBP)

# Save the life-course hypotheses and covariates as a matrix
x_hypos <- data_dist_inc_diaBP %>%
  select(-diaBP)

x_hypos <- as.matrix(x_hypos)


## Run the Lasso model using GLMNET. alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives covariates (income, age, sex and white) weightings of '0', so are always included in the model (default is 1)
mod_dist_inc_diaBP <- glmnet(x_hypos, data_dist_inc_diaBP$diaBP,
                               alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod_dist_inc_diaBP

# Plot these results
plot(mod_dist_inc_diaBP)

# Look at the variables included at each step
coef(mod_dist_inc_diaBP, s = max(mod_dist_inc_diaBP$lambda[mod_dist_inc_diaBP$df == 4])); min(mod_dist_inc_diaBP$dev.ratio[mod_dist_inc_diaBP$df == 4])

coef(mod_dist_inc_diaBP, s = max(mod_dist_inc_diaBP$lambda[mod_dist_inc_diaBP$df == 5])); min(mod_dist_inc_diaBP$dev.ratio[mod_dist_inc_diaBP$df == 5])

coef(mod_dist_inc_diaBP, s = max(mod_dist_inc_diaBP$lambda[mod_dist_inc_diaBP$df == 6])); min(mod_dist_inc_diaBP$dev.ratio[mod_dist_inc_diaBP$df == 6])

coef(mod_dist_inc_diaBP, s = max(mod_dist_inc_diaBP$lambda[mod_dist_inc_diaBP$df == 7])); min(mod_dist_inc_diaBP$dev.ratio[mod_dist_inc_diaBP$df == 7])


### Visual inspection of results (although just looking at the deviance ratios there doesn't seem to be much of an effect of access to green space at all)

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 4, nrow = 0))
df

for (i in 1:length(mod_dist_inc_diaBP$df)) {
  #print(i)
  new_covars <- attributes(which(mod_dist_inc_diaBP$beta[, i] != 0))$names
  new_deviance <- mod_dist_inc_diaBP$dev.ratio[i]
  new_varNum <- mod_dist_inc_diaBP$df[i]
  new_lambda <- mod_dist_inc_diaBP$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if new covars different to old covars
  if (new_varNum != old_varNum) {
    new_addition <- setdiff(new_covars, old_covars) # Find the new covariate that's been added
    dev_diff <- (new_deviance - old_deviance) * 100 # Diff in deviance value between current and previous lambda value
    dev_diff <- round(dev_diff, 3)
    temp <- cbind(new_addition, dev_diff, new_varNum, new_lambda) # Merge with template data frame
    df <- rbind(df, temp)
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevDiff", "VarNum", "Lambda")
df

# If two variables added at once, combine names together (else impossible to read on plot)
df$Variables[df$Variables == "age_f7"] <- "covars"
df

df <- df[df$Variables != "male", ]
df <- df[df$Variables != "white", ]
df <- df[df$Variables != "lowIncome", ]
df

# Make a var to show number of steps where variables added
df$steps <- 1:nrow(df)
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. Need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod_dist_inc_diaBP$lambda, mod_dist_inc_diaBP$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_inc_diaBP$lambda)), ylim = c(0, max(mod_dist_inc_diaBP$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod_dist_inc_diaBP$log_lambda <- log(mod_dist_inc_diaBP$lambda)
mod_dist_inc_diaBP

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, although same issues as above re. interpreting the first variable added (doesn't mean is a strong predictor of outcome).
plot(mod_dist_inc_diaBP$log_lambda, mod_dist_inc_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_inc_diaBP$log_lambda)), ylim = c(0.004, max(mod_dist_inc_diaBP$dev.ratio)))
text(df$log_lambda, 0.004, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_distanceIncomeDiaBP.pdf", height = 7, width = 11)
plot(mod_dist_dep_diaBP$log_lambda, mod_dist_dep_diaBP$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod_dist_dep_diaBP$log_lambda)), ylim = c(0.004, max(mod_dist_dep_diaBP$dev.ratio)))
text(df$log_lambda, 0.004, labels = df$Variables, srt = 90, adj = 0)
dev.off()


## From these results, would seem to be that nothing is really going on here - Will do a likelihood ratio test to see whether inclusion of first parameter(s) added (int1) increases model fit of standard linear regression model
base_mod <- lm(data_dist_inc_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + x_hypos[, "white"] + 
                 x_hypos[, "lowIncome"])
summary(base_mod)

param1_mod <- lm(data_dist_inc_diaBP$diaBP ~ x_hypos[, "age_f7"] + x_hypos[, "male"] + 
                   x_hypos[, "white"] + x_hypos[, "lowIncome"] + x_hypos[, "int1"])
summary(param1_mod)

# No association here.
anova(base_mod, param1_mod)


## Alternative method using cross-validated lasso to find find 'optimal' model for out-of-sample prediction.
mod.cv <- cv.glmnet(x_hypos, data_dist_inc_diaBP$diaBP, 
                    alpha = 1, penalty.factor = (c(0, 0, 0, 0, rep(1, ncol(x_hypos) - 4))))
mod.cv

# Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters is at top of plot)
plot(mod.cv) 

# Parameters in the 1SE and 'optimal' models - Both just contains just the 4 parameters included by default (age, sex, ethnicity and education), suggesting no association with access to green space exposure
coef(mod.cv, s = mod.cv$lambda.1se)
coef(mod.cv, s = mod.cv$lambda.min)


### All methods match up relatively well, and indicate little association between distance to green space in childhood and diastolic blood pressure at age 7 (when controlling for ethnicity and household income).



############## Next, need to repeat for the two remaining exposures: NVDI value and access to garden





