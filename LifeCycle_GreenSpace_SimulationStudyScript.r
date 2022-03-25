### Simulation study script for structured lifecourse models with confounders and interactions

### Created 17/3/2022 by Dan Major-Smith
### R version 4.0.4

### Aim: To test the power of this structured life course approach with interactions to select the 'true' model under varying conditions. Will explore all combinations of the various parameters:
#   - Sample sizes: 1,000 vs 10,000
#   - Exposures: Binary (access to green space; yes/no) vs continuous (distance to green space), and centered vs uncentered. Also want to vary the correlation between exposures to see how collinearity impacts the power of the lasso to detect the true model
#   - Outcome: Binary (overweight/obese) vs continuous (BMI)

## For these simulations, will use the same set-up as in the example simulation script, with 'access to green space' as the exposure measured at three time points, cardiometabolic health as the outcome (BMI/obesity), and SEP as a confounder/interaction term. SEP causes access to green space and the outcome (lower BMI/obesity if higher SEP), while the interaction between SEP and the first green space time-point also causes the outcome (access to green space causes lower BMI/obesity, but only in those from lower SEP backgrounds).


#######################################################################################################
#### Clear workspace and install/load packages

rm(list = ls())

#install.packages("glmnet")
#install.packages("dagitty")

library(glmnet)
library(dagitty)


## Working directory
setwd("C:\\Users\\ds16565\\OneDrive - University of Bristol\\MyFiles-Migrated\\Documents\\Projects\\Lifecycle\\SimResults")


## Hypothetical DAG - The model here is that SEP is a confounder, but also interacts with first green space exposure in pregnancy to impact BMI at age 7
dag <- dagitty('dag {
                SEP [pos = "0,0"]
                Green1 [pos = "1,1"]
                Green2 [pos = "1,1.5"]
                Green3 [pos = "1,2"]
                BMI [pos = "2,0"]
                SEP_Green1_int [pos = "1.5,0.5"]
                
                SEP -> Green1
                SEP -> Green2
                SEP -> Green3
                SEP -> BMI
                SEP -> SEP_Green1_int
                Green1 -> SEP_Green1_int
                SEP_Green1_int -> BMI
                Green1 -> Green2
                Green2 -> Green3
                }')
plot(dag)




#######################################################################################################
#### First, need to set the parameters for the simulation study. 

## This will require a bit of trial-and-error to get the correlations between the exposures right - Will aim for r = 0.3 for the 'less correlated' simulation and r = 0.9 for the more correlated simulation (with this set of values, the correlations between the critical period interaction terms should be even higher, as anything >0.9 could be considered a high degree of collinearity). To fix these values, will use a single large simulated dataset with 1 million observations, to remove random variability.


#### Start with working out the binary exposures
set.seed(4321)
n <- 1000000

# SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
high_sep <- rbinom(n, 1, 0.5)
table(high_sep, useNA = "ifany")

# First green space exposure in pregnancy - Caused by SEP (higher SEP = more exposure to green space; three times the odds) - Approx. 50% near green space
green1_p <- plogis(log(0.6) + (log(3) * high_sep))
green1 <- rbinom(n, 1, green1_p)
table(green1, useNA = "ifany")

table(high_sep, green1)

rm(green1_p)

# Second green space exposure at age 4 - Caused by Green1 and SEP
green2_p <- plogis(log(0.3) + (log(3) * high_sep) + (log(3) * green1))
green2 <- rbinom(n, 1, green2_p)
table(green2, useNA = "ifany")

table(high_sep, green2)
table(green1, green2)

rm(green2_p)

# Third green space exposure at age 7 - Caused by Green2 and SEP
green3_p <- plogis(log(0.3) + (log(3) * high_sep) + (log(3) * green2))
green3 <- rbinom(n, 1, green3_p)
table(green3, useNA = "ifany")

table(high_sep, green3)
table(green1, green3)
table(green2, green3)

rm(green3_p)

## How correlated are the green space exposures - These are ~0.3, so an OR of 3 would be sensible for this 'lower' rate of exposure correlation (and means that none of the correlations between the interactions is above 0.9; see below)
cor(green1, green2)
cor(green2, green3)

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green1 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green1 here.
bmi <- 25 + (-4 * high_sep) + (-2 * green1) + (2 * high_sep * green1) + rnorm(n, 0, 3)
summary(bmi)

## descriptive stats split by possible combinations of SEP and green1 to check simulation worked
summary(bmi[high_sep == 1 & green1 == 1]) # High SEP and access to green space in preg - Lowest BMI
summary(bmi[high_sep == 1 & green1 == 0]) # High SEP and no access to green space in preg - Lowest BMI
summary(bmi[high_sep == 0 & green1 == 1]) # Low SEP and access to green space in preg - Middle BMI
summary(bmi[high_sep == 0 & green1 == 0]) # Low SEP and no access to green space in preg - Highest BMI

# Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be pretty much null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))


### Encode the life course hypotheses

# Critical period at first time point only
crit1 <- green1

# Interaction between SEP and first time point
int1 <- crit1 * high_sep

# Critical period at second time point only
crit2 <- green2

# Interaction between SEP and second time point
int2 <- crit2 * high_sep 

# Critical period at third time point only
crit3 <- green3

# Interaction between SEP and third time point
int3 <- crit3 * high_sep

# Linear accumulation of all exposures
accumulation <- green1 + green2 + green3

# Interaction between SEP and cumulative exposure
int_accum <- high_sep * accumulation

# Increase in access to green space from time 1 to time 2
green_inc12 <- (1 - green1) * green2

# Increase in access to green space from time 1 to time 2, with an interaction with SEP
green_inc12_int <- green_inc12 * high_sep

# Decrease in access to green space from time 1 to time 2
green_dec12 <- (1 - green2) * green1

# Decrease in access to green space from time 1 to time 2, with an interaction with SEP
green_dec12_int <- green_dec12 * high_sep

# Increase in access to green space from time 2 to time 3
green_inc23 <- (1 - green2) * green3

# Increase in access to green space from time 2 to time 3, with an interaction with SEP
green_inc23_int <- green_inc23 * high_sep

# Decrease in access to green space from time 2 to time 3
green_dec23 <- (1 - green3) * green2

# Decrease in access to green space from time 2 to time 3, with an interaction with SEP
green_dec23_int <- green_dec23 * high_sep


## Combine all these into one matrix, including SEP as a confounder
x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_inc12, 
                 green_inc12_int, green_dec12, green_dec12_int, green_inc23, green_inc23_int, green_dec23,
                 green_dec23_int)
head(x_hypos)


## Check correlation matrix of all these hypotheses - Using the parameters above, all correlations are below 0.9, so hopefully not a cause for concern
dim(x_hypos)
cor(x_hypos[,2:17])
cor(x_hypos[,2:17]) > 0.9


### Next, repeat the above but aiming for an r = 0.9, so that exposures are highly correlated, and that the interaction terms will be highly correlated too.
set.seed(4322)

# SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
high_sep <- rbinom(n, 1, 0.5)
table(high_sep, useNA = "ifany")

# First green space exposure in pregnancy - Caused by SEP (higher SEP = more exposure to green space; three times the odds) - Approx. 50% near green space
green1_p <- plogis(log(0.6) + (log(3) * high_sep))
green1 <- rbinom(n, 1, green1_p)
table(green1, useNA = "ifany")

table(high_sep, green1)

rm(green1_p)

# Second green space exposure at age 4 - Caused by Green1 and SEP
green2_p <- plogis(log(0.01) + (log(3) * high_sep) + (log(500) * green1))
green2 <- rbinom(n, 1, green2_p)
table(green2, useNA = "ifany")

table(high_sep, green2)
table(green1, green2)

rm(green2_p)

# Third green space exposure at age 7 - Caused by Green2 and SEP
green3_p <- plogis(log(0.01) + (log(3) * high_sep) + (log(500) * green2))
green3 <- rbinom(n, 1, green3_p)
table(green3, useNA = "ifany")

table(high_sep, green3)
table(green1, green3)
table(green2, green3)

rm(green3_p)

## How correlated are the green space exposures - These are ~0.9, so an OR of 500 would appear sensible for this 'higher' rate of exposure correlation (and means that quite a few of the correlations between the interactions are above 0.9; see below). This is a *very* extreme scenario, and assumes that there is relatively little change in access to green space between the time-points.
cor(green1, green2)
cor(green2, green3)

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green1 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green1 here.
bmi <- 25 + (-4 * high_sep) + (-2 * green1) + (2 * high_sep * green1) + rnorm(n, 0, 3)
summary(bmi)

## descriptive stats split by possible combinations of SEP and green1 to check simulation worked
summary(bmi[high_sep == 1 & green1 == 1]) # High SEP and access to green space in preg - Lowest BMI
summary(bmi[high_sep == 1 & green1 == 0]) # High SEP and no access to green space in preg - Lowest BMI
summary(bmi[high_sep == 0 & green1 == 1]) # Low SEP and access to green space in preg - Middle BMI
summary(bmi[high_sep == 0 & green1 == 0]) # Low SEP and no access to green space in preg - Highest BMI

# Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be pretty much null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))


### Encode the life course hypotheses

# Critical period at first time point only
crit1 <- green1

# Interaction between SEP and first time point
int1 <- crit1 * high_sep

# Critical period at second time point only
crit2 <- green2

# Interaction between SEP and second time point
int2 <- crit2 * high_sep 

# Critical period at third time point only
crit3 <- green3

# Interaction between SEP and third time point
int3 <- crit3 * high_sep

# Linear accumulation of all exposures
accumulation <- green1 + green2 + green3

# Interaction between SEP and cumulative exposure
int_accum <- high_sep * accumulation

# Increase in access to green space from time 1 to time 2
green_inc12 <- (1 - green1) * green2

# Increase in access to green space from time 1 to time 2, with an interaction with SEP
green_inc12_int <- green_inc12 * high_sep

# Decrease in access to green space from time 1 to time 2
green_dec12 <- (1 - green2) * green1

# Decrease in access to green space from time 1 to time 2, with an interaction with SEP
green_dec12_int <- green_dec12 * high_sep

# Increase in access to green space from time 2 to time 3
green_inc23 <- (1 - green2) * green3

# Increase in access to green space from time 2 to time 3, with an interaction with SEP
green_inc23_int <- green_inc23 * high_sep

# Decrease in access to green space from time 2 to time 3
green_dec23 <- (1 - green3) * green2

# Decrease in access to green space from time 2 to time 3, with an interaction with SEP
green_dec23_int <- green_dec23 * high_sep


## Combine all these into one matrix, including SEP as a confounder
x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_inc12, 
                 green_inc12_int, green_dec12, green_dec12_int, green_inc23, green_inc23_int, green_dec23,
                 green_dec23_int)
head(x_hypos)


## Check correlation matrix of all these hypotheses - Quite a few of the correlations are above 0.9, especially for accumulation, accumulation interaction, and critical period interactions (similar to the real-world ALSPAC scenario).
dim(x_hypos)
cor(x_hypos[,2:17])
cor(x_hypos[,2:17]) > 0.9



#### Next, want to repeat the above for continuous 'distance to green space' exposures (establishing 'low' and 'high' correlation scenarios)
set.seed(4323)

# SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
high_sep <- rbinom(n, 1, 0.5)
table(high_sep, useNA = "ifany")

# First distance to green space exposure in pregnancy - Caused by SEP (higher SEP = more exposure to green space; 100m closer if higher SEP) - Mean distance at ~300m
green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40)
summary(green1)

# Second green space exposure at age 4 - Caused by Green1 and SEP
green2 <- 240 + (-50 * high_sep) + (0.3 * green1) + rnorm(n, 0, 40)
summary(green2)

# Third green space exposure at age 7 - Caused by Green2 and SEP
green3 <- 240 + (-50 * high_sep) + (0.3 * green2) + rnorm(n, 0, 40)
summary(green3)

## How correlated are the green space exposures - These are ~0.5, should be fine for this 'lower' rate of exposure correlation (however, given that this is a continuous exposure many of the interaction terms still have an correlation above 0.9; see below)
cor(green1, green2)
cor(green2, green3)

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green1 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green1 here. Have chosen 'green1' parameter to be 0.02 increase in BMI per unit increase in green space, as the green1 standard deviation is ~50, and two times this should cover most of the variation in green space distance, making parameters broadly comparable to to binary green space effect of 2 BMI units (see: Gelman, A. (2008). Scaling regression inputs by dividing by two standard deviations. Statistics in medicine, 27(15), 2865-2873.)
bmi <- 25 + (-4 * high_sep) + (-0.02 * green1) + (0.02 * high_sep * green1) + rnorm(n, 0, 3)
summary(bmi)

# Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be pretty much null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))


### Encode the life course hypotheses

# Critical period at first time point only
crit1 <- green1

# Interaction between SEP and first time point
int1 <- crit1 * high_sep

# Critical period at second time point only
crit2 <- green2

# Interaction between SEP and second time point
int2 <- crit2 * high_sep 

# Critical period at third time point only
crit3 <- green3

# Interaction between SEP and third time point
int3 <- crit3 * high_sep

# Linear accumulation of all exposures
accumulation <- green1 + green2 + green3 / 3

# Interaction between SEP and cumulative exposure
int_accum <- high_sep * accumulation

# Change in distance to green space between time 1 to time 2
green_ch12 <- green2 - green1

# Change in distance to green space between time 1 to time 2, with an interaction with SEP
green_ch12_int <- green_ch12 * high_sep

# Change in distance to green space between time 2 to time 3
green_ch23 <- green3 - green2

# Change in distance to green space between time 2 to time 3, with an interaction with SEP
green_ch23_int <- green_ch23 * high_sep



## Combine all these into one matrix, including SEP as a confounder
x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_ch12, 
                 green_ch12_int, green_ch23, green_ch23_int)
head(x_hypos)


## Check correlation matrix of all these hypotheses - Using the parameters above, even in this 'low' correlation scenario, as the exposures are continuous many of the interactions have correlations above 0.9
dim(x_hypos)
cor(x_hypos[,2:13])
cor(x_hypos[,2:13]) > 0.9


## If we centre the hypotheses, then this should drop considerably

# Critical period at first time point only
crit1 <- green1
crit1 <- crit1 - mean(crit1)

# Interaction between SEP and first time point
int1 <- crit1 * high_sep

# Critical period at second time point only
crit2 <- green2
crit2 <- crit2 - mean(crit2)

# Interaction between SEP and second time point
int2 <- crit2 * high_sep 

# Critical period at third time point only
crit3 <- green3
crit3 <- crit3 - mean(crit3)

# Interaction between SEP and third time point
int3 <- crit3 * high_sep

# Linear accumulation of all exposures
accumulation <- green1 + green2 + green3 / 3
accumulation <- accumulation - mean(accumulation)

# Interaction between SEP and cumulative exposure
int_accum <- high_sep * accumulation

# Change in distance to green space between time 1 to time 2
green_ch12 <- green2 - green1
green_ch12 <- green_ch12 - mean(green_ch12)

# Change in distance to green space between time 1 to time 2, with an interaction with SEP
green_ch12_int <- green_ch12 * high_sep

# Change in distance to green space between time 2 to time 3
green_ch23 <- green3 - green2
green_ch23 <- green_ch23 - mean(green_ch23)

# Change in distance to green space between time 2 to time 3, with an interaction with SEP
green_ch23_int <- green_ch23 * high_sep

## Combine all these into one matrix, including SEP as a confounder
x_hypos_c <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_ch12, 
                 green_ch12_int, green_ch23, green_ch23_int)
head(x_hypos_c)


## Check correlation matrix of all these centered hypotheses - Now all the correlations are < 0.9, so the centering is effective.
dim(x_hypos_c)
cor(x_hypos_c[,2:13])
cor(x_hypos_c[,2:13]) > 0.9



### Finally, we want to repeat the above continuous 'distance to green space' exposures, but this time for the the 'high' correlation scenarios
set.seed(4324)

# SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
high_sep <- rbinom(n, 1, 0.5)
table(high_sep, useNA = "ifany")

# First distance to green space exposure in pregnancy - Caused by SEP (higher SEP = more exposure to green space; 100m closer if higher SEP) - Mean distance at ~300m
green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40)
summary(green1)

# Second green space exposure at age 4 - Caused by Green1 and SEP
green2 <- 50 + (-50 * high_sep) + (0.9 * green1) + rnorm(n, 0, 20)
summary(green2)

# Third green space exposure at age 7 - Caused by Green2 and SEP
green3 <- 50 + (-50 * high_sep) + (0.9 * green2) + rnorm(n, 0, 20)
summary(green3)

## How correlated are the green space exposures - These are ~0.9, should be fine for this 'higher' rate of exposure correlation
cor(green1, green2)
cor(green2, green3)

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green1 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green1 here.
bmi <- 25 + (-4 * high_sep) + (-0.02 * green1) + (0.02 * high_sep * green1) + rnorm(n, 0, 3)
summary(bmi)

# Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be pretty much null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))


### Encode the life course hypotheses

# Critical period at first time point only
crit1 <- green1

# Interaction between SEP and first time point
int1 <- crit1 * high_sep

# Critical period at second time point only
crit2 <- green2

# Interaction between SEP and second time point
int2 <- crit2 * high_sep 

# Critical period at third time point only
crit3 <- green3

# Interaction between SEP and third time point
int3 <- crit3 * high_sep

# Linear accumulation of all exposures
accumulation <- green1 + green2 + green3 / 3

# Interaction between SEP and cumulative exposure
int_accum <- high_sep * accumulation

# Change in distance to green space between time 1 to time 2
green_ch12 <- green2 - green1

# Change in distance to green space between time 1 to time 2, with an interaction with SEP
green_ch12_int <- green_ch12 * high_sep

# Change in distance to green space between time 2 to time 3
green_ch23 <- green3 - green2

# Change in distance to green space between time 2 to time 3, with an interaction with SEP
green_ch23_int <- green_ch23 * high_sep


## Combine all these into one matrix, including SEP as a confounder
x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_ch12, 
                 green_ch12_int, green_ch23, green_ch23_int)
head(x_hypos)

## Check correlation matrix of all these hypotheses - All the accumulation and critical period interaction terms are really highly correlated, just as we'd expect
dim(x_hypos)
cor(x_hypos[,2:13])
cor(x_hypos[,2:13]) > 0.9


## Centering should lower some of the correlations, but should still be relatively high
# Critical period at first time point only
crit1 <- green1
crit1 <- crit1 - mean(crit1)

# Interaction between SEP and first time point
int1 <- crit1 * high_sep

# Critical period at second time point only
crit2 <- green2
crit2 <- crit2 - mean(crit2)

# Interaction between SEP and second time point
int2 <- crit2 * high_sep 

# Critical period at third time point only
crit3 <- green3
crit3 <- crit3 - mean(crit3)

# Interaction between SEP and third time point
int3 <- crit3 * high_sep

# Linear accumulation of all exposures
accumulation <- green1 + green2 + green3 / 3
accumulation <- accumulation - mean(accumulation)

# Interaction between SEP and cumulative exposure
int_accum <- high_sep * accumulation

# Change in distance to green space between time 1 to time 2
green_ch12 <- green2 - green1
green_ch12 <- green_ch12 - mean(green_ch12)

# Change in distance to green space between time 1 to time 2, with an interaction with SEP
green_ch12_int <- green_ch12 * high_sep

# Change in distance to green space between time 2 to time 3
green_ch23 <- green3 - green2
green_ch23 <- green_ch23 - mean(green_ch23)

# Change in distance to green space between time 2 to time 3, with an interaction with SEP
green_ch23_int <- green_ch23 * high_sep

## Combine all these into one matrix, including SEP as a confounder
x_hypos_c <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_ch12, 
                   green_ch12_int, green_ch23, green_ch23_int)
head(x_hypos_c)


## Check correlation matrix of all these centered hypotheses - They are lower than in the uncentered data, but still quite high, and a few ar >0.9 (esp. for accumulation and int_accum with the critical period terms and their interactions).
dim(x_hypos_c)
cor(x_hypos_c[,2:13])
cor(x_hypos_c[,2:13]) > 0.9


### Now that this has established the parameters to use for the simulations, we can start running the actual simulations!



###########################################################################################################
#### Actual simulations

## First, set-up a function to perform the simulations

# n_sims = Number of simulations (any integer; default = 1000)
# sampleSize = Sample size for each simulation (any integer; default = 1000)
# Exposure = Binary or continuous exposures (either "Binary" or "Cont"; default = "Binary")
# Centered = Whether to center the exposures or not (either "No" or "Yes"; default = "No")
# Collinear = How collinear the exposures are (either "Low" or "High"; default = "Low")
# Outcome = Binary or continuous outcome (either "Binary" or "Cont"; default = "Cont")
# Output = Whether to print the model output or not (TRUE or FALSE; default = FALSE)
lasso_sim <- function(n_sims = 1000, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                      Outcome = "Cont", Output = FALSE) {
  
  # Initiate vectors to save results from this simulation to (i.e., whether method identified correct model or not)
  LR_res_temp <- rep(NA, n_sims)
  LR_res_temp_crit1 <- rep(NA, n_sims) # Store if 'crit1' main effect in final model
  LR_res_temp_int1 <- rep(NA, n_sims) # Store if 'int1' interaction effect in final model
  LR_res_temp_crit1int1 <- rep(NA, n_sims) # Store if 'crit1' and 'int1' in final model
  LR_res_temp_crit1int1extra <- rep(NA, n_sims) # Store if 'crit1' and 'int1', plus extra vars, in final model
  CV_res_temp <- rep(NA, n_sims)
  CV_res_temp_crit1 <- rep(NA, n_sims) # Store if 'crit1' main effect in final model
  CV_res_temp_int1 <- rep(NA, n_sims) # Store if 'int1' interaction effect in final model
  CV_res_temp_crit1int1 <- rep(NA, n_sims) # Store if 'crit1' and 'int1' in final model
  CV_res_temp_crit1int1extra <- rep(NA, n_sims) # Store if 'crit1' and 'int1', plus extra vars, in final model
  
  for (i in 1:n_sims) {
    
    print(paste0("Simulation number: ", i))
    
    ## Start by simulating the data
    n <- sampleSize
    
    ## If exposure is binary and collinearity low
    if (Exposure == "Binary" & Collinear == "Low") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1_p <- plogis(log(0.6) + (log(3) * high_sep)) # First green space exposure in pregnancy
      green1 <- rbinom(n, 1, green1_p)
      green2_p <- plogis(log(0.3) + (log(3) * high_sep) + (log(3) * green1)) # Second green space exposure
      green2 <- rbinom(n, 1, green2_p)
      green3_p <- plogis(log(0.3) + (log(3) * high_sep) + (log(3) * green2)) # Third green space exposure
      green3 <- rbinom(n, 1, green3_p)
      bmi <- 25 + (-4 * high_sep) + (-2 * green1) + (2 * high_sep * green1) + rnorm(n, 0, 3) # Cont. BMI outcome
    }
    
    ## If exposure is binary and collinearity high
    if (Exposure == "Binary" & Collinear == "High") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1_p <- plogis(log(0.6) + (log(3) * high_sep)) # First green space exposure in pregnancy
      green1 <- rbinom(n, 1, green1_p)
      green2_p <- plogis(log(0.01) + (log(3) * high_sep) + (log(500) * green1)) # Second green space exposure
      green2 <- rbinom(n, 1, green2_p)
      green3_p <- plogis(log(0.01) + (log(3) * high_sep) + (log(500) * green2)) # Third green space exposure
      green3 <- rbinom(n, 1, green3_p)
      bmi <- 25 + (-4 * high_sep) + (-2 * green1) + (2 * high_sep * green1) + rnorm(n, 0, 3) # Cont. BMI outcome
    }
    
    ## If exposure is continuous and collinearity low
    if (Exposure == "Cont" & Collinear == "Low") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40) # First green space distance exposure in pregnancy
      green2 <- 240 + (-50 * high_sep) + (0.3 * green1) + rnorm(n, 0, 40) # Second green space distance exposure
      green3 <- 240 + (-50 * high_sep) + (0.3 * green2) + rnorm(n, 0, 40) # Third green space distance exposure
      bmi <- 25 + (-4 * high_sep) + (-0.02 * green1) + (0.02 * high_sep * green1) + rnorm(n, 0, 3) # Cont. BMI outcome
    }
    
    ## If exposure is continuous and collinearity high
    if (Exposure == "Cont" & Collinear == "High") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40) # First green space distance exposure in pregnancy
      green2 <- 50 + (-50 * high_sep) + (0.9 * green1) + rnorm(n, 0, 20) # Second green space distance exposure
      green3 <- 50 + (-50 * high_sep) + (0.9 * green2) + rnorm(n, 0, 20) # Third green space distance exposure
      bmi <- 25 + (-4 * high_sep) + (-0.02 * green1) + (0.02 * high_sep * green1) + rnorm(n, 0, 3) # Cont. BMI outcome
    }
    
    # If outcome is binary
    if (Outcome == "Binary") {
      overweight <- ifelse(bmi > 25, 1, 0) # Code BMI to binary 'overweight' variable if > 25
    }
    
    
    ## Encode the life course hypotheses
    
    # Critical periods (same for binary and continuous exposures)
    crit1 <- green1 # Critical period at first time point only
    if (Centered == "Yes") {
      crit1 <- crit1 - mean(crit1) # Center this variable
    }
    int1 <- crit1 * high_sep # Interaction between SEP and first time point
    
    crit2 <- green2 # Critical period at second time point only
    if (Centered == "Yes") {
      crit2 <- crit2 - mean(crit2) # Center this variable
    }
    int2 <- crit2 * high_sep # Interaction between SEP and second time point
    
    crit3 <- green3 # Critical period at third time point only
    if (Centered == "Yes") {
      crit3 <- crit3 - mean(crit3) # Center this variable
    }
    int3 <- crit3 * high_sep # Interaction between SEP and third time point
    
    # Accumulation (different for binary and continuous exposures)
    if (Exposure == "Binary") {
      accumulation <- green1 + green2 + green3 # Linear accumulation of all exposures
      if (Centered == "Yes") {
        accumulation <- accumulation - mean(accumulation) # Center this variable
      }
      int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
    }
    
    if (Exposure == "Cont") {
      accumulation <- green1 + green2 + green3 / 3 # Linear accumulation of all exposures
      if (Centered == "Yes") {
        accumulation <- accumulation - mean(accumulation) # Center this variable
      }
      int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
    }
    
    # Change (different for binary and continuous exposures)
    if (Exposure == "Binary") {
      green_inc12 <- (1 - green1) * green2 # Increase from time 1 to time 2
      if (Centered == "Yes") {
        green_inc12 <- green_inc12 - mean(green_inc12) # Center this variable
      }
      green_inc12_int <- green_inc12 * high_sep # Increase from time 1 to time 2, with an interaction with SEP
      
      green_dec12 <- (1 - green2) * green1 # Decrease from time 1 to time 2
      if (Centered == "Yes") {
        green_dec12 <- green_dec12 - mean(green_dec12) # Center this variable
      }
      green_dec12_int <- green_dec12 * high_sep # Decrease from time 1 to time 2, with an interaction with SEP
      
      green_inc23 <- (1 - green2) * green3 # Increase from time 2 to time 3
      if (Centered == "Yes") {
        green_inc23 <- green_inc23 - mean(green_inc23) # Center this variable
      }
      green_inc23_int <- green_inc23 * high_sep # Increase from time 2 to time 3, with an interaction with SEP
      
      green_dec23 <- (1 - green3) * green2 # Decrease from time 2 to time 3
      if (Centered == "Yes") {
        green_dec23 <- green_dec23 - mean(green_dec23) # Center this variable
      }
      green_dec23_int <- green_dec23 * high_sep # Decrease from time 2 to time 3, with an interaction with SEP
    }
    
    if (Exposure == "Cont") {
      green_ch12 <- green2 - green1 # Change from time 1 to time 2
      if (Centered == "Yes") {
        green_ch12 <- green_ch12 - mean(green_ch12) # Center this variable
      }
      green_ch12_int <- green_ch12 * high_sep # Change from time 1 to time 2, with an interaction with SEP
      
      green_ch23 <- green3 - green2 # Change from time 2 to time 3
      if (Centered == "Yes") {
        green_ch23 <- green_ch23 - mean(green_ch23) # Center this variable
      }
      green_ch23_int <- green_ch23 * high_sep # Change from time 2 to time 3, with an interaction with SEP
    }
    
    
    ## Combine all these into one matrix, including SEP as a confounder
    
    # If binary exposures
    if (Exposure == "Binary") {
      x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_inc12, 
                       green_inc12_int, green_dec12, green_dec12_int, green_inc23, green_inc23_int, green_dec23,
                       green_dec23_int)
    }
    
    # If continuous exposures
    if (Exposure == "Cont") {
      x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_ch12, 
                       green_ch12_int, green_ch23, green_ch23_int)
    }
    
    ## Run the model using glmnet
    
    # If continuous outcome
    if (Outcome == "Cont") {
      mod <- glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    # If binary outcome
    if (Outcome == "Binary") {
      mod <- glmnet(x_hypos, overweight, alpha = 1, family = "binomial", 
                    penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    ## Loop over each step of the lasso, running a likelihood ratio test against the previous model when there is a change (skipping the first step, as this is the SEP-only model). Stop this loop when p > 0.05
    old_covars <- "high_sep" # Set up the SEP-only covariate list
    x_hypos_old <- as.matrix(x_hypos[, colnames(x_hypos) %in% old_covars]) # Matrix of just baseline SEP covariate
    p <- 0 # Initialise the p-value to be 0
    
    for (j in 2:length(mod$df)) {
      new_covars <- attributes(which(mod$beta[, j] != 0))$names # Extract the covariates at each time-point
      
      # See whether covariates have changed at this step
      if (setequal(old_covars, new_covars) == FALSE) {
        x_hypos_new <- x_hypos[, colnames(x_hypos) %in% new_covars] # Matrix of the covariates in updated lasso
        
        # See whether these new covariates improve model fit
        
        # If outcome is continuous
        if (Outcome == "Cont") {
          mod_old <- lm(bmi ~ x_hypos_old)
          mod_new <- lm(bmi ~ x_hypos_new)
          p <- as.data.frame(anova(mod_old, mod_new))[2, 6] # Extracting the p-value from the LR test
        }
        
        # If outcome is binary
        if (Outcome == "Binary") {
          mod_old <- glm(overweight ~ x_hypos_old, family = "binomial")
          mod_new <- glm(overweight ~ x_hypos_new, family = "binomial")
          p <- as.data.frame(anova(mod_old, mod_new, test = "Chisq"))[2, 5] # Extracting the p-value from the LR test
        }
        
        # If p-value is > 0.05, exit the loop (or if p is NA, which can happen if the old and new model have the same number of parameters)
        if (p > 0.05 | is.na(p)) {
          break
        }
        
        # Update the covariate list and the 'old' covariate matrix
        x_hypos_old <- x_hypos_new
        old_covars <- new_covars
      }
      
      # Again, if p-value is > 0.05, exit the loop
      if (p > 0.05 | is.na(p)) {
        break
      }
      
    }
    
    # Print the covariates in the best-fitting LR model if want to print this
    if (Output == TRUE) {
      print(paste0("LR covariates:"))
      print(old_covars)
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so (or NA if unable to calculate p-value)
    LR_res_temp[i] <- ifelse(setequal(old_covars, target_covars) == TRUE, 1, 0)
    LR_res_temp[i] <- ifelse(is.na(p), NA, LR_res_temp[i])
    
    # Store if 'crit1', 'int1' and both in final model
    LR_res_temp_crit1[i] <- ifelse("crit1" %in% old_covars == TRUE, 1, 0)
    LR_res_temp_crit1[i] <- ifelse(is.na(p), NA, LR_res_temp_crit1[i])
    
    LR_res_temp_int1[i] <- ifelse("int1" %in% old_covars == TRUE, 1, 0)
    LR_res_temp_int1[i] <- ifelse(is.na(p), NA, LR_res_temp_int1[i])
    
    LR_res_temp_crit1int1[i] <- ifelse("crit1" %in% old_covars == TRUE & "int1" %in% old_covars == TRUE, 1, 0)
    LR_res_temp_crit1int1[i] <- ifelse(is.na(p), NA, LR_res_temp_crit1int1[i])
    
    LR_res_temp_crit1int1extra[i] <- ifelse(LR_res_temp_crit1int1[i] == 1 & LR_res_temp[i] == 0, 1, 0)
    LR_res_temp_crit1int1extra[i] <- ifelse(is.na(p), NA, LR_res_temp_crit1int1extra[i])
    
    
    ### Next, want to summarise the 1SE cross-validated lasso model, and see whether that corresponds to the correct model or not
    
    # If outcome is continuous
    if (Outcome == "Cont") {
      mod.cv <- cv.glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    # If outcome is binary
    if (Outcome == "Binary") {
      mod.cv <- cv.glmnet(x_hypos, overweight, alpha = 1, family = "binomial", 
                          penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    # Extract the non-zero covariates into a vector
    cv_covars_temp <- as.matrix(coef(mod.cv, s = mod.cv$lambda.1se))
    cv_covars_temp <- as.matrix(cv_covars_temp[!rownames(cv_covars_temp) %in% "(Intercept)", ]) # Drop the intercept
    cv_covars_temp <- cv_covars_temp[cv_covars_temp != 0, ] # Drop all zero coefficients
    
    cv_covars <- attributes(cv_covars_temp)$names # Store all the non-zero hypotheses
    
    # Print the covariates in the best-fitting cross-validated model
    if (Output == TRUE) {
      print(paste0("CV covariates:"))
      print(cv_covars)
      print("")
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
    CV_res_temp[i] <- ifelse(setequal(cv_covars, target_covars) == TRUE, 1, 0)
    
    # Store if 'crit1', 'int1' and both in final model
    CV_res_temp_crit1[i] <- ifelse("crit1" %in% cv_covars == TRUE, 1, 0)
    CV_res_temp_int1[i] <- ifelse("int1" %in% cv_covars == TRUE, 1, 0)
    CV_res_temp_crit1int1[i] <- ifelse("crit1" %in% cv_covars == TRUE & "int1" %in% cv_covars == TRUE, 1, 0)
    CV_res_temp_crit1int1extra[i] <- ifelse(CV_res_temp_crit1int1[i] == 1 & CV_res_temp[i] == 0, 1, 0)
    
  }
  
  # Store the summaries of these results to transfer to the main results table
  res <- data.frame(LR_Nworked = sum(!is.na(LR_res_temp)),
                    LR_Ncorrect = sum(LR_res_temp, na.rm = TRUE),
                    LR_propcorrect = round(sum(LR_res_temp, na.rm = TRUE) / sum(!is.na(LR_res_temp)) * 100, 2),
                    LR_crit1correct = round(sum(LR_res_temp_crit1, na.rm = TRUE) / sum(!is.na(LR_res_temp)) * 100, 2),
                    LR_int1correct = round(sum(LR_res_temp_int1, na.rm = TRUE) / sum(!is.na(LR_res_temp)) * 100, 2),
                    LR_crit1int1correct = round(sum(LR_res_temp_crit1int1, na.rm = TRUE) / 
                                                  sum(!is.na(LR_res_temp)) * 100, 2),
                    LR_crit1int1extra = round(sum(LR_res_temp_crit1int1extra, na.rm = TRUE) / 
                                                sum(!is.na(LR_res_temp)) * 100, 2),
                    CV_propcorrect = round(sum(CV_res_temp) / n_sims * 100, 2),
                    CV_crit1correct = round(sum(CV_res_temp_crit1) / n_sims * 100, 2),
                    CV_int1correct = round(sum(CV_res_temp_int1) / n_sims * 100, 2),
                    CV_crit1int1correct = round(sum(CV_res_temp_crit1int1) / n_sims * 100, 2),
                    CV_crit1int1extra = round(sum(CV_res_temp_crit1int1extra) / n_sims * 100, 2))
  return(res)
  
}


## Next, the number of simulations per combination of parameters (1,000), the target hypotheses we simulated are the true model, and set up a data frame to store the results in
n_sims <- 1000
set.seed(9876)

target_covars <- c("high_sep", "crit1", "int1")

results <- data.frame(sampleSize = rep(c(1000, 10000), 16),
                      Exposure = rep(c(rep("Binary", 8), rep("Cont", 8)), 2),
                      Centered = rep(c("No", "No", "Yes", "Yes"), 8),
                      Collinear = rep(c(rep("Low", 4), rep("High", 4)), 4),
                      Outcome = c(rep("Cont", 16), rep("Binary", 16)),
                      LR_Nworked = rep(NA, nrow(results)),
                      LR_Ncorrect = rep(NA, nrow(results)),
                      LR_propcorrect = rep(NA, nrow(results)),
                      LR_crit1correct = rep(NA, nrow(results)),
                      LR_int1correct = rep(NA, nrow(results)),
                      LR_crit1int1correct = rep(NA, nrow(results)),
                      LR_crit1int1extra = rep(NA, nrow(results)),
                      CV_propcorrect = rep(NA, nrow(results)),
                      CV_crit1correct = rep(NA, nrow(results)),
                      CV_int1correct = rep(NA, nrow(results)),
                      CV_crit1int1correct = rep(NA, nrow(results)),
                      CV_crit1int1extra = rep(NA, nrow(results)))

results


# Also calculate time taken to run script
start_time <- Sys.time()


### First simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 1
results[k, 6:17] <- res


### Second simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 2
results[k, 6:17] <- res


### Third simulation: Sample size = 1000; binary exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 3
results[k, 6:17] <- res


### Fourth simulation: Sample size = 10000; binary exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 4
results[k, 6:17] <- res


### Fifth simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 5
results[k, 6:17] <- res


### Sixth simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 6
results[k, 6:17] <- res


### Seventh simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 7
results[k, 6:17] <- res


### Eighth simulation: Sample size = 10000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 8
results[k, 6:17] <- res


### Ninth simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 9
results[k, 6:17] <- res


### Tenth simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 10
results[k, 6:17] <- res


### Eleventh simulation: Sample size = 1000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 11
results[k, 6:17] <- res


### Twelfth simulation: Sample size = 10000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 12
results[k, 6:17] <- res


### 13th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 13
results[k, 6:17] <- res


### 14th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 14
results[k, 6:17] <- res


### 15th simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 15
results[k, 6:17] <- res


### 16th simulation: Sample size = 10000; continuous exposure; centered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 16
results[k, 6:17] <- res


### 17th simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 17
results[k, 6:17] <- res


### 18th simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 18
results[k, 6:17] <- res


### 19th simulation: Sample size = 1000; binary exposure; centered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 19
results[k, 6:17] <- res


### 20th simulation: Sample size = 10000; binary exposure; centered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 20
results[k, 6:17] <- res


### 21st simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 21
results[k, 6:17] <- res


### 22nd simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 22
results[k, 6:17] <- res


### 23rd simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 23
results[k, 6:17] <- res


### 24th simulation: Sample size = 10000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 24
results[k, 6:17] <- res


### 25th simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 25
results[k, 6:17] <- res


### 26th simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 26
results[k, 6:17] <- res


### 27th simulation: Sample size = 1000; continuous exposure; centered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 27
results[k, 6:17] <- res


### 28th simulation: Sample size = 10000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 28
results[k, 6:17] <- res


### 29th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 29
results[k, 6:17] <- res


### 30th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 30
results[k, 6:17] <- res


### 31st simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 31
results[k, 6:17] <- res


### 32nd simulation: Sample size = 10000; continuous exposure; centered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 32
results[k, 6:17] <- res


# Time taken to run script (took about 40 mins for 100 simulations, so should be about 6 hours for the full 1,000 simulations)
end_time <- Sys.time()

end_time - start_time


## Save the results table
results

write.csv(results, file = "simulationResults.csv", row.names = FALSE)


### Get some results from this
#results <- read.csv("simulationResults.csv")
#head(results)

## Plotting the percentage correct for each model
plot(results$LR_propcorrect, results$CV_propcorrect, pch = 16, xlim = c(0, 100), ylim = c(0, 100),  
     xlab = "% of correct models (LR test)", ylab = "% of correct models (1SE CV lasso)")

pdf(file = "LRbyCV_simulationResults_WithChangeVars.pdf", height = 5, width = 7)
plot(results$LR_propcorrect, results$CV_propcorrect, pch = 16, xlim = c(0, 100), ylim = c(0, 100), 
     xlab = "% of correct models (LR test)", ylab = "% of correct models (1SE CV lasso)")
dev.off()

## Top performing models
# LR method
results[order(-results$LR_propcorrect), c(1:10)]

# 1SE method
results[order(-results$CV_propcorrect), c(1:5, 13:15)]

## Worst performing methods
# LR method
results[order(results$LR_propcorrect), c(1:10)]

# 1SE method
results[order(results$CV_propcorrect), c(1:5, 13:15)]


### Summary stats, split by factors we've varied

## Sample size
# LR test
summary(results[results$sampleSize == 1000, 8:12])
summary(results[results$sampleSize == 10000, 8:12])

# 1SE method
summary(results[results$sampleSize == 1000, 13:17])
summary(results[results$sampleSize == 10000, 13:17])

## Binary vs continuous exposure
# LR test
summary(results[results$Exposure == "Binary", 8:12])
summary(results[results$Exposure == "Cont", 8:12])

# 1SE method
summary(results[results$Exposure == "Binary", 13:17])
summary(results[results$Exposure == "Cont", 13:17])

## Exposures centered or not
# LR test
summary(results[results$Centered == "No", 8:12])
summary(results[results$Centered == "Yes", 8:12])

# 1SE method
summary(results[results$Centered == "No", 13:17])
summary(results[results$Centered == "Yes", 13:17])

## Differences in centering by whether exposure is continuous or binary? Some improvement if center binary variables, but larger improvement if center continuous variables (unsurprisingly!)
# LR test
summary(results[results$Exposure == "Binary" & results$Centered == "No", 8:12])
summary(results[results$Exposure == "Binary" & results$Centered == "Yes", 8:12])
summary(results[results$Exposure == "Cont" & results$Centered == "No", 8:12])
summary(results[results$Exposure == "Cont" & results$Centered == "Yes", 8:12])

# 1SE method
summary(results[results$Exposure == "Binary" & results$Centered == "No", 13:17])
summary(results[results$Exposure == "Binary" & results$Centered == "Yes", 13:17])
summary(results[results$Exposure == "Cont" & results$Centered == "No", 13:17])
summary(results[results$Exposure == "Cont" & results$Centered == "Yes", 13:17])

## Exposures collinear or not
# LR test
summary(results[results$Collinear == "Low", 8:12])
summary(results[results$Collinear == "High", 8:12])

# 1SE method
summary(results[results$Collinear == "Low", 13:17])
summary(results[results$Collinear == "High", 13:17])

## Differences in collinearity by whether exposure is continuous or binary? Differences in collinearity approximately the same, regardless of whether continuous or binary exposure
# LR test
summary(results[results$Exposure == "Binary" & results$Collinear == "Low", 8:12])
summary(results[results$Exposure == "Binary" & results$Collinear == "High", 8:12])
summary(results[results$Exposure == "Cont" & results$Collinear == "Low", 8:12])
summary(results[results$Exposure == "Cont" & results$Collinear == "High", 8:12])

# 1SE method
summary(results[results$Exposure == "Binary" & results$Collinear == "Low", 13:17])
summary(results[results$Exposure == "Binary" & results$Collinear == "High", 13:17])
summary(results[results$Exposure == "Cont" & results$Collinear == "Low", 13:17])
summary(results[results$Exposure == "Cont" & results$Collinear == "High", 13:17])

## Binary or continuous outcome
# LR test
summary(results[results$Outcome == "Binary", 8:12])
summary(results[results$Outcome == "Cont", 8:12])

# 1SE method
summary(results[results$Outcome == "Binary", 13:17])
summary(results[results$Outcome == "Cont", 13:17])

## Combination of binary or continuous outcome and binary or continuous exposure? Continuous exposure and binary outcome performs much worse than all other methods.
# LR test
summary(results[results$Exposure == "Binary" & results$Outcome == "Binary", 8:12])
summary(results[results$Exposure == "Binary" & results$Outcome == "Cont", 8:12])
summary(results[results$Exposure == "Cont" & results$Outcome == "Binary", 8:12])
summary(results[results$Exposure == "Cont" & results$Outcome == "Cont", 8:12])

# 1SE method
summary(results[results$Exposure == "Binary" & results$Outcome == "Binary", 13:17])
summary(results[results$Exposure == "Binary" & results$Outcome == "Cont", 13:17])
summary(results[results$Exposure == "Cont" & results$Outcome == "Binary", 13:17])
summary(results[results$Exposure == "Cont" & results$Outcome == "Cont", 13:17])



########################################################################################################
########################################################################################################
### Will also run some further simulations to explore whether removing some hypotheses could make the results less biased. Here, I will drop all of the 'change' variables - This does assume that 'change' is not associated with the outcome (which we know it isn't in this simulated example), but it is a big assumption if we do not know the true causal model

## First, set-up a function to perform the simulations

# n_sims = Number of simulations (any integer; default = 1000)
# sampleSize = Sample size for each simulation (any integer; default = 1000)
# Exposure = Binary or continuous exposures (either "Binary" or "Cont"; default = "Binary")
# Centered = Whether to center the exposures or not (either "No" or "Yes"; default = "No")
# Collinear = How collinear the exposures are (either "Low" or "High"; default = "Low")
# Outcome = Binary or continuous outcome (either "Binary" or "Cont"; default = "Cont")
# Output = Whether to print the model output or not (TRUE or FALSE; default = FALSE)
lasso_sim_reduced <- function(n_sims = 1000, sampleSize = 1000, Exposure = "Binary", Centered = "No", 
                              Collinear = "Low", Outcome = "Cont", Output = FALSE) {
  
  # Initiate vectors to save results from this simulation to (i.e., whether method identified correct model or not)
  LR_res_temp <- rep(NA, n_sims)
  LR_res_temp_crit1 <- rep(NA, n_sims) # Store if 'crit1' main effect in final model
  LR_res_temp_int1 <- rep(NA, n_sims) # Store if 'int1' interaction effect in final model
  LR_res_temp_crit1int1 <- rep(NA, n_sims) # Store if 'crit1' and 'int1' in final model
  LR_res_temp_crit1int1extra <- rep(NA, n_sims) # Store if 'crit1' and 'int1', plus extra vars, in final model
  CV_res_temp <- rep(NA, n_sims)
  CV_res_temp_crit1 <- rep(NA, n_sims) # Store if 'crit1' main effect in final model
  CV_res_temp_int1 <- rep(NA, n_sims) # Store if 'int1' interaction effect in final model
  CV_res_temp_crit1int1 <- rep(NA, n_sims) # Store if 'crit1' and 'int1' in final model
  CV_res_temp_crit1int1extra <- rep(NA, n_sims) # Store if 'crit1' and 'int1', plus extra vars, in final model
  
  for (i in 1:n_sims) {
    
    print(paste0("Simulation number: ", i))
    
    ## Start by simulating the data
    n <- sampleSize
    
    ## If exposure is binary and collinearity low
    if (Exposure == "Binary" & Collinear == "Low") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1_p <- plogis(log(0.6) + (log(3) * high_sep)) # First green space exposure in pregnancy
      green1 <- rbinom(n, 1, green1_p)
      green2_p <- plogis(log(0.3) + (log(3) * high_sep) + (log(3) * green1)) # Second green space exposure
      green2 <- rbinom(n, 1, green2_p)
      green3_p <- plogis(log(0.3) + (log(3) * high_sep) + (log(3) * green2)) # Third green space exposure
      green3 <- rbinom(n, 1, green3_p)
      bmi <- 25 + (-4 * high_sep) + (-2 * green1) + (2 * high_sep * green1) + rnorm(n, 0, 3) # Cont. BMI outcome
    }
    
    ## If exposure is binary and collinearity high
    if (Exposure == "Binary" & Collinear == "High") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1_p <- plogis(log(0.6) + (log(3) * high_sep)) # First green space exposure in pregnancy
      green1 <- rbinom(n, 1, green1_p)
      green2_p <- plogis(log(0.01) + (log(3) * high_sep) + (log(500) * green1)) # Second green space exposure
      green2 <- rbinom(n, 1, green2_p)
      green3_p <- plogis(log(0.01) + (log(3) * high_sep) + (log(500) * green2)) # Third green space exposure
      green3 <- rbinom(n, 1, green3_p)
      bmi <- 25 + (-4 * high_sep) + (-2 * green1) + (2 * high_sep * green1) + rnorm(n, 0, 3) # Cont. BMI outcome
    }
    
    ## If exposure is continuous and collinearity low
    if (Exposure == "Cont" & Collinear == "Low") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40) # First green space distance exposure in pregnancy
      green2 <- 240 + (-50 * high_sep) + (0.3 * green1) + rnorm(n, 0, 40) # Second green space distance exposure
      green3 <- 240 + (-50 * high_sep) + (0.3 * green2) + rnorm(n, 0, 40) # Third green space distance exposure
      bmi <- 25 + (-4 * high_sep) + (-0.02 * green1) + (0.02 * high_sep * green1) + rnorm(n, 0, 3) # Cont. BMI outcome
    }
    
    ## If exposure is continuous and collinearity high
    if (Exposure == "Cont" & Collinear == "High") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40) # First green space distance exposure in pregnancy
      green2 <- 50 + (-50 * high_sep) + (0.9 * green1) + rnorm(n, 0, 20) # Second green space distance exposure
      green3 <- 50 + (-50 * high_sep) + (0.9 * green2) + rnorm(n, 0, 20) # Third green space distance exposure
      bmi <- 25 + (-4 * high_sep) + (-0.02 * green1) + (0.02 * high_sep * green1) + rnorm(n, 0, 3) # Cont. BMI outcome
    }
    
    # If outcome is binary
    if (Outcome == "Binary") {
      overweight <- ifelse(bmi > 25, 1, 0) # Code BMI to binary 'overweight' variable if > 25
    }
    
    
    ## Encode the life course hypotheses
    
    # Critical periods (same for binary and continuous exposures)
    crit1 <- green1 # Critical period at first time point only
    if (Centered == "Yes") {
      crit1 <- crit1 - mean(crit1) # Center this variable
    }
    int1 <- crit1 * high_sep # Interaction between SEP and first time point
    
    crit2 <- green2 # Critical period at second time point only
    if (Centered == "Yes") {
      crit2 <- crit2 - mean(crit2) # Center this variable
    }
    int2 <- crit2 * high_sep # Interaction between SEP and second time point
    
    crit3 <- green3 # Critical period at third time point only
    if (Centered == "Yes") {
      crit3 <- crit3 - mean(crit3) # Center this variable
    }
    int3 <- crit3 * high_sep # Interaction between SEP and third time point
    
    # Accumulation (different for binary and continuous exposures)
    if (Exposure == "Binary") {
      accumulation <- green1 + green2 + green3 # Linear accumulation of all exposures
      if (Centered == "Yes") {
        accumulation <- accumulation - mean(accumulation) # Center this variable
      }
      int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
    }
    
    if (Exposure == "Cont") {
      accumulation <- green1 + green2 + green3 / 3 # Linear accumulation of all exposures
      if (Centered == "Yes") {
        accumulation <- accumulation - mean(accumulation) # Center this variable
      }
      int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
    }
    
    ## Combine all these into one matrix, including SEP as a confounder
    
    # If binary exposures
    if (Exposure == "Binary") {
      x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum)
    }
    
    # If continuous exposures
    if (Exposure == "Cont") {
      x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum)
    }
    
    ## Run the model using glmnet
    
    # If continuous outcome
    if (Outcome == "Cont") {
      mod <- glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    # If binary outcome
    if (Outcome == "Binary") {
      mod <- glmnet(x_hypos, overweight, alpha = 1, family = "binomial", 
                    penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    ## Loop over each step of the lasso, running a likelihood ratio test against the previous model when there is a change (skipping the first step, as this is the SEP-only model). Stop this loop when p > 0.05
    old_covars <- "high_sep" # Set up the SEP-only covariate list
    x_hypos_old <- as.matrix(x_hypos[, colnames(x_hypos) %in% old_covars]) # Matrix of just baseline SEP covariate
    p <- 0 # Initialise the p-value to be 0
    
    for (j in 2:length(mod$df)) {
      new_covars <- attributes(which(mod$beta[, j] != 0))$names # Extract the covariates at each time-point
      
      # See whether covariates have changed at this step
      if (setequal(old_covars, new_covars) == FALSE) {
        x_hypos_new <- x_hypos[, colnames(x_hypos) %in% new_covars] # Matrix of the covariates in updated lasso
        
        # See whether these new covariates improve model fit
        
        # If outcome is continuous
        if (Outcome == "Cont") {
          mod_old <- lm(bmi ~ x_hypos_old)
          mod_new <- lm(bmi ~ x_hypos_new)
          p <- as.data.frame(anova(mod_old, mod_new))[2, 6] # Extracting the p-value from the LR test
        }
        
        # If outcome is binary
        if (Outcome == "Binary") {
          mod_old <- glm(overweight ~ x_hypos_old, family = "binomial")
          mod_new <- glm(overweight ~ x_hypos_new, family = "binomial")
          p <- as.data.frame(anova(mod_old, mod_new, test = "Chisq"))[2, 5] # Extracting the p-value from the LR test
        }
        
        # If p-value is > 0.05, exit the loop (or if p is NA, which can happen if the old and new model have the same number of parameters)
        if (p > 0.05 | is.na(p)) {
          break
        }
        
        # Update the covariate list and the 'old' covariate matrix
        x_hypos_old <- x_hypos_new
        old_covars <- new_covars
      }
      
      # Again, if p-value is > 0.05, exit the loop
      if (p > 0.05 | is.na(p)) {
        break
      }
      
    }
    
    # Print the covariates in the best-fitting LR model if want to print this
    if (Output == TRUE) {
      print(paste0("LR covariates:"))
      print(old_covars)
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so (or NA if unable to calculate p-value)
    LR_res_temp[i] <- ifelse(setequal(old_covars, target_covars) == TRUE, 1, 0)
    LR_res_temp[i] <- ifelse(is.na(p), NA, LR_res_temp[i])
    
    # Store if 'crit1', 'int1' and both in final model
    LR_res_temp_crit1[i] <- ifelse("crit1" %in% old_covars == TRUE, 1, 0)
    LR_res_temp_crit1[i] <- ifelse(is.na(p), NA, LR_res_temp_crit1[i])
    
    LR_res_temp_int1[i] <- ifelse("int1" %in% old_covars == TRUE, 1, 0)
    LR_res_temp_int1[i] <- ifelse(is.na(p), NA, LR_res_temp_int1[i])
    
    LR_res_temp_crit1int1[i] <- ifelse("crit1" %in% old_covars == TRUE & "int1" %in% old_covars == TRUE, 1, 0)
    LR_res_temp_crit1int1[i] <- ifelse(is.na(p), NA, LR_res_temp_crit1int1[i])
    
    LR_res_temp_crit1int1extra[i] <- ifelse(LR_res_temp_crit1int1[i] == 1 & LR_res_temp[i] == 0, 1, 0)
    LR_res_temp_crit1int1extra[i] <- ifelse(is.na(p), NA, LR_res_temp_crit1int1extra[i])
    
    
    ### Next, want to summarise the 1SE cross-validated lasso model, and see whether that corresponds to the correct model or not
    
    # If outcome is continuous
    if (Outcome == "Cont") {
      mod.cv <- cv.glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    # If outcome is binary
    if (Outcome == "Binary") {
      mod.cv <- cv.glmnet(x_hypos, overweight, alpha = 1, family = "binomial", 
                          penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    # Extract the non-zero covariates into a vector
    cv_covars_temp <- as.matrix(coef(mod.cv, s = mod.cv$lambda.1se))
    cv_covars_temp <- as.matrix(cv_covars_temp[!rownames(cv_covars_temp) %in% "(Intercept)", ]) # Drop the intercept
    cv_covars_temp <- cv_covars_temp[cv_covars_temp != 0, ] # Drop all zero coefficients
    
    cv_covars <- attributes(cv_covars_temp)$names # Store all the non-zero hypotheses
    
    # Print the covariates in the best-fitting cross-validated model
    if (Output == TRUE) {
      print(paste0("CV covariates:"))
      print(cv_covars)
      print("")
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
    CV_res_temp[i] <- ifelse(setequal(cv_covars, target_covars) == TRUE, 1, 0)
    
    # Store if 'crit1', 'int1' and both in final model
    CV_res_temp_crit1[i] <- ifelse("crit1" %in% cv_covars == TRUE, 1, 0)
    CV_res_temp_int1[i] <- ifelse("int1" %in% cv_covars == TRUE, 1, 0)
    CV_res_temp_crit1int1[i] <- ifelse("crit1" %in% cv_covars == TRUE & "int1" %in% cv_covars == TRUE, 1, 0)
    CV_res_temp_crit1int1extra[i] <- ifelse(CV_res_temp_crit1int1[i] == 1 & CV_res_temp[i] == 0, 1, 0)
    
  }
  
  # Store the summaries of these results to transfer to the main results table
  res <- data.frame(LR_Nworked = sum(!is.na(LR_res_temp)),
                    LR_Ncorrect = sum(LR_res_temp, na.rm = TRUE),
                    LR_propcorrect = round(sum(LR_res_temp, na.rm = TRUE) / sum(!is.na(LR_res_temp)) * 100, 2),
                    LR_crit1correct = round(sum(LR_res_temp_crit1, na.rm = TRUE) / sum(!is.na(LR_res_temp)) * 100, 2),
                    LR_int1correct = round(sum(LR_res_temp_int1, na.rm = TRUE) / sum(!is.na(LR_res_temp)) * 100, 2),
                    LR_crit1int1correct = round(sum(LR_res_temp_crit1int1, na.rm = TRUE) / 
                                                  sum(!is.na(LR_res_temp)) * 100, 2),
                    LR_crit1int1extra = round(sum(LR_res_temp_crit1int1extra, na.rm = TRUE) / 
                                                sum(!is.na(LR_res_temp)) * 100, 2),
                    CV_propcorrect = round(sum(CV_res_temp) / n_sims * 100, 2),
                    CV_crit1correct = round(sum(CV_res_temp_crit1) / n_sims * 100, 2),
                    CV_int1correct = round(sum(CV_res_temp_int1) / n_sims * 100, 2),
                    CV_crit1int1correct = round(sum(CV_res_temp_crit1int1) / n_sims * 100, 2),
                    CV_crit1int1extra = round(sum(CV_res_temp_crit1int1extra) / n_sims * 100, 2))
  return(res)
  
}


## Next, the number of simulations per combination of parameters (1,000), the target hypotheses we simulated are the true model, and set up a data frame to store the results in
n_sims <- 1000
set.seed(6789)

target_covars <- c("high_sep", "crit1", "int1")

results_reduced <- data.frame(sampleSize = rep(c(1000, 10000), 16),
                      Exposure = rep(c(rep("Binary", 8), rep("Cont", 8)), 2),
                      Centered = rep(c("No", "No", "Yes", "Yes"), 8),
                      Collinear = rep(c(rep("Low", 4), rep("High", 4)), 4),
                      Outcome = c(rep("Cont", 16), rep("Binary", 16)),
                      LR_Nworked = rep(NA, nrow(results)),
                      LR_Ncorrect = rep(NA, nrow(results)),
                      LR_propcorrect = rep(NA, nrow(results)),
                      LR_crit1correct = rep(NA, nrow(results)),
                      LR_int1correct = rep(NA, nrow(results)),
                      LR_crit1int1correct = rep(NA, nrow(results)),
                      LR_crit1int1extra = rep(NA, nrow(results)),
                      CV_propcorrect = rep(NA, nrow(results)),
                      CV_crit1correct = rep(NA, nrow(results)),
                      CV_int1correct = rep(NA, nrow(results)),
                      CV_crit1int1correct = rep(NA, nrow(results)),
                      CV_crit1int1extra = rep(NA, nrow(results)))

results_reduced


# Also calculate time taken to run script
start_time <- Sys.time()


### First simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", 
                         Collinear = "Low",  Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 1
results_reduced[k, 6:17] <- res


### Second simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 2
results_reduced[k, 6:17] <- res


### Third simulation: Sample size = 1000; binary exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 3
results_reduced[k, 6:17] <- res


### Fourth simulation: Sample size = 10000; binary exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "Low",  Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 4
results_reduced[k, 6:17] <- res


### Fifth simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 5
results_reduced[k, 6:17] <- res


### Sixth simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 6
results_reduced[k, 6:17] <- res


### Seventh simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 7
results_reduced[k, 6:17] <- res


### Eighth simulation: Sample size = 10000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 8
results_reduced[k, 6:17] <- res


### Ninth simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 9
results_reduced[k, 6:17] <- res


### Tenth simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 10
results_reduced[k, 6:17] <- res


### Eleventh simulation: Sample size = 1000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 11
results_reduced[k, 6:17] <- res


### Twelfth simulation: Sample size = 10000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 12
results_reduced[k, 6:17] <- res


### 13th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 13
results_reduced[k, 6:17] <- res


### 14th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 14
results_reduced[k, 6:17] <- res


### 15th simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 15
results_reduced[k, 6:17] <- res


### 16th simulation: Sample size = 10000; continuous exposure; centered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 16
results_reduced[k, 6:17] <- res


### 17th simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 17
results_reduced[k, 6:17] <- res


### 18th simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 18
results_reduced[k, 6:17] <- res


### 19th simulation: Sample size = 1000; binary exposure; centered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 19
results_reduced[k, 6:17] <- res


### 20th simulation: Sample size = 10000; binary exposure; centered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 20
results_reduced[k, 6:17] <- res


### 21st simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 21
results_reduced[k, 6:17] <- res


### 22nd simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 22
results_reduced[k, 6:17] <- res


### 23rd simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 23
results_reduced[k, 6:17] <- res


### 24th simulation: Sample size = 10000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 24
results_reduced[k, 6:17] <- res


### 25th simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 25
results_reduced[k, 6:17] <- res


### 26th simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 26
results_reduced[k, 6:17] <- res


### 27th simulation: Sample size = 1000; continuous exposure; centered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 27
results_reduced[k, 6:17] <- res


### 28th simulation: Sample size = 10000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes",
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 28
results_reduced[k, 6:17] <- res


### 29th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 29
results_reduced[k, 6:17] <- res


### 30th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 30
results_reduced[k, 6:17] <- res


### 31st simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 31
results_reduced[k, 6:17] <- res


### 32nd simulation: Sample size = 10000; continuous exposure; centered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 32
results_reduced[k, 6:17] <- res


# Time taken to run script (took about 40 mins for 100 simulations, so should be about 6 hours for the full 1,000 simulations)
end_time <- Sys.time()

end_time - start_time


## Save the results table
results_reduced

write.csv(results_reduced, file = "simulationResults_reduced.csv", row.names = FALSE)


### Get some results from this
#results_reduced <- read.csv("simulationResults_reduced.csv")
#head(results_reduced)

## Plotting the percentage correct for each model
plot(results_reduced$LR_propcorrect, results_reduced$CV_propcorrect, pch = 16, xlim = c(0, 100), ylim = c(0, 100), 
     xlab = "% of correct models (LR test)", ylab = "% of correct models (1SE CV lasso)")

pdf(file = "LRbyCV_simulationResults_noChangeVars.pdf", height = 5, width = 7)
plot(results_reduced$LR_propcorrect, results_reduced$CV_propcorrect, pch = 16, xlim = c(0, 100), ylim = c(0, 100), 
     xlab = "% of correct models (LR test)", ylab = "% of correct models (1SE CV lasso)")
dev.off()

## Top performing models
# LR method
results_reduced[order(-results_reduced$LR_propcorrect), c(1:10)]

# 1SE method
results_reduced[order(-results_reduced$CV_propcorrect), c(1:5, 13:15)]

## Worst performing methods
# LR method
results_reduced[order(results_reduced$LR_propcorrect), c(1:10)]

# 1SE method
results_reduced[order(results_reduced$CV_propcorrect), c(1:5, 13:15)]


### Summary stats, split by factors we've varied

## Sample size
# LR test
summary(results_reduced[results_reduced$sampleSize == 1000, 8:12])
summary(results_reduced[results_reduced$sampleSize == 10000, 8:12])

# 1SE method
summary(results_reduced[results_reduced$sampleSize == 1000, 13:17])
summary(results_reduced[results_reduced$sampleSize == 10000, 13:17])

## Binary vs continuous exposure
# LR test
summary(results_reduced[results_reduced$Exposure == "Binary", 8:12])
summary(results_reduced[results_reduced$Exposure == "Cont", 8:12])

# 1SE method
summary(results_reduced[results_reduced$Exposure == "Binary", 13:17])
summary(results_reduced[results_reduced$Exposure == "Cont", 13:17])

## Exposures centered or not
# LR test
summary(results_reduced[results_reduced$Centered == "No", 8:12])
summary(results_reduced[results_reduced$Centered == "Yes", 8:12])

# 1SE method
summary(results_reduced[results_reduced$Centered == "No", 13:17])
summary(results_reduced[results_reduced$Centered == "Yes", 13:17])

## Differences in centering by whether exposure is continuous or binary? Slight improvement if center binary variables, but larger improvement if center continuous variables (unsurprisingly!)
# LR test
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Centered == "No", 8:12])
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Centered == "Yes", 8:12])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Centered == "No", 8:12])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Centered == "Yes", 8:12])

# 1SE method
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Centered == "No", 13:17])
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Centered == "Yes", 13:17])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Centered == "No", 13:17])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Centered == "Yes", 13:17])

## Exposures collinear or not
# LR test
summary(results_reduced[results_reduced$Collinear == "Low", 8:12])
summary(results_reduced[results_reduced$Collinear == "High", 8:12])

# 1SE method
summary(results_reduced[results_reduced$Collinear == "Low", 13:17])
summary(results_reduced[results_reduced$Collinear == "High", 13:17])

## Differences in collinearity by whether exposure is continuous or binary? Collinearity slightly worse if a binary exposure, but this could be due to the way variables were coded, potentially (as hard to say whether collinearity between binary and continuous variables is equivalent)
# LR test
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Collinear == "Low", 8:12])
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Collinear == "High", 8:12])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Collinear == "Low", 8:12])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Collinear == "High", 8:12])

# 1SE method
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Collinear == "Low", 13:17])
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Collinear == "High", 13:17])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Collinear == "Low", 13:17])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Collinear == "High", 13:17])

## Binary or continuous outcome
# LR test
summary(results_reduced[results_reduced$Outcome == "Binary", 8:12])
summary(results_reduced[results_reduced$Outcome == "Cont", 8:12])

# 1SE method
summary(results_reduced[results_reduced$Outcome == "Binary", 13:17])
summary(results_reduced[results_reduced$Outcome == "Cont", 13:17])

## Combination of binary or continuous outcome and binary or continuous exposure? Continuous exposure and binary outcome performs much worse than all other methods.
# LR test
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Outcome == "Binary", 8:12])
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Outcome == "Cont", 8:12])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Outcome == "Binary", 8:12])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Outcome == "Cont", 8:12])

# 1SE method
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Outcome == "Binary", 13:17])
summary(results_reduced[results_reduced$Exposure == "Binary" & results_reduced$Outcome == "Cont", 13:17])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Outcome == "Binary", 13:17])
summary(results_reduced[results_reduced$Exposure == "Cont" & results_reduced$Outcome == "Cont", 13:17])


## Comparing overall performance of simulations with vs without the 'change' variables. Results are better (although still lots of errors)
# LR test
summary(results[, 8:12])
summary(results_reduced[, 8:12])

# 1SE method
summary(results[, 13:17])
summary(results_reduced[, 13:17])
