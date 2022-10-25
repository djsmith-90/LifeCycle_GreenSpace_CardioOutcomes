### Simulation study script for structured lifecourse models with confounders and interactions

### Created 17/3/2022 by Dan Major-Smith
### R version 4.0.4

### Aim: To test the power of this structured life course approach with interactions to select the 'true' interaction term under varying conditions. Will explore all combinations of the various parameters:
#   - Sample sizes: 1,000 vs 10,000
#   - Exposures: Binary (access to green space; yes/no) vs continuous (distance to green space), and centered vs uncentered. Also want to vary the correlation between exposures to see how collinearity impacts the power of the lasso to detect the true model
#   - Outcome: Binary (overweight/obese) vs continuous (BMI)

## For these simulations, will use the same set-up as in the example simulation script, with 'access to green space' as the exposure measured at three time points, cardiometabolic health as the outcome (BMI/obesity), and SEP as a confounder/interaction term. SEP causes access to green space and the outcome (lower BMI/obesity if higher SEP), while the interaction between SEP and the most recent green space time-point also causes the outcome (access to green space causes lower BMI/obesity, but only in those from lower SEP backgrounds).

## In addition to this scenario, we will also vary the strength of the interaction term, to explore how this impacts the power to detect the interaction, as well as varying the specific life course interaction (i.e., the main model will explore an interaction with first first critical period, while other simulations will examine interactions with accumulation and change, to see whether this impacts conclusions). Given the number of simulations to run, and that each simulation takes about 6 hours, processing time on a standard laptop may be a bit prohibitive, so will use this script to set-up and test the simulation/code, and then run the actual simulations using University of Bristol's High Performance Computing suite.


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
                SEP_Green3_int [pos = "1.33,0.5"]
                
                SEP -> Green1
                SEP -> Green2
                SEP -> Green3
                SEP -> BMI
                SEP -> SEP_Green3_int
                Green3 -> SEP_Green3_int
                Green3 -> BMI
                SEP_Green3_int -> BMI
                Green1 -> Green2
                Green2 -> Green3
                }')
plot(dag)




#######################################################################################################
#### First, need to set the parameters for the simulation study. 

## This will require a bit of trial-and-error to get the correlations between the exposures right - Will aim for r = 0.3 for the 'less correlated' simulation and r = 0.9 for the more correlated simulation (with this set of values, the correlations between the critical period interaction terms should be even higher, as anything >0.9 could be considered a high degree of collinearity). To fix these values, will use a single large simulated dataset with 1 million observations, to remove most of the random variability.


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

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green3 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green3 here.
bmi <- 25 + (-4 * high_sep) + (-2 * green3) + (2 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi)
sd(bmi)

## descriptive stats split by possible combinations of SEP and green1 to check simulation worked
summary(bmi[high_sep == 1 & green3 == 1]) # High SEP and recent access to green space - Lowest BMI
summary(bmi[high_sep == 1 & green3 == 0]) # High SEP and no recent access to green space - Lowest BMI
summary(bmi[high_sep == 0 & green3 == 1]) # Low SEP and recent access to green space - Middle BMI
summary(bmi[high_sep == 0 & green3 == 0]) # Low SEP and no recent access to green space - Highest BMI

# Check the 'true' model, which is SEP as confounder, interaction with green 3, and main effect of green 3. green1 and green1 should also be pretty much null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))

AIC(lm(bmi ~ high_sep + green1 + green2 + green3)); AIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
BIC(lm(bmi ~ high_sep + green1 + green2 + green3)); BIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_noInt <- 25 + (-4 * high_sep) + (-2 * green3) + (0 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_noInt)
sd(bmi); sd(bmi_noInt)

# Very small interaction
bmi_vSmallInt <- 25 + (-4 * high_sep) + (-2 * green3) + (0.5 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_vSmallInt)
sd(bmi); sd(bmi_vSmallInt)

# Small interaction
bmi_smallInt <- 25 + (-4 * high_sep) + (-2 * green3) + (1 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_smallInt)
sd(bmi); sd(bmi_smallInt)

# Large interaction
bmi_largeInt <- 25 + (-4 * high_sep) + (-2 * green3) + (3 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_largeInt)
sd(bmi); sd(bmi_largeInt)

# very large interaction
bmi_vLargeInt <- 25 + (-4 * high_sep) + (-2 * green3) + (4 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_vLargeInt)
sd(bmi); sd(bmi_vLargeInt)


### Encode the life course hypotheses (here just using the original BMI values)

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


## Also want to create BMI outcome using different life course interaction to see if results differ. First, do interaction with accumulation (using the same model specs as above, so that access to green space only lowers BMI if from low SEP, but no difference if high SEP). As accumulation is continuous (from 0 to 3), will take two-thirds off the previous binary critical period coefficients so are approximately on the same scale. Make sure the BMI value is similar to the BMI value above for the equivalent critical period model
bmi_accum <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (0.67 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_accum)
sd(bmi); sd(bmi_accum)

# Check the 'true' model, which is SEP as confounder, interaction with accumulation, and main effect of accumulation. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi_accum ~ high_sep + accumulation))
summary(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))
anova(lm(bmi_accum ~ high_sep + accumulation), lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))

AIC(lm(bmi_accum ~ high_sep + accumulation)); AIC(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))
BIC(lm(bmi_accum ~ high_sep + accumulation)); BIC(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_accum_noInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (0 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_noInt)
sd(bmi_accum); sd(bmi_accum_noInt)

# Very small interaction
bmi_accum_vSmallInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (0.167 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_vSmallInt)
sd(bmi_accum); sd(bmi_accum_vSmallInt)

# Small interaction
bmi_accum_smallInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (0.33 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_smallInt)
sd(bmi_accum); sd(bmi_accum_smallInt)

# Large interaction
bmi_accum_largeInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (1 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_largeInt)
sd(bmi_accum); sd(bmi_accum_largeInt)

# very large interaction
bmi_accum_vLargeInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (1.33 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_vLargeInt)
sd(bmi_accum); sd(bmi_accum_vLargeInt)


## Next, construct a BMI outcome with one of the 'change' hypotheses as the main effect and interacting variable. Will say that an increase in green space from time 2 to time 3 lowers BMI, but only in those from lower SEP. As change is binary, will use the same parameters as from the critical period interaction model
bmi_change <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (2 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change)

# Check the 'true' model, which is SEP as confounder, interaction with green_inc23, and main effect of green_inc23. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi_change ~ high_sep + green_inc23))
summary(lm(bmi_change ~ high_sep + green_inc23 + high_sep*green_inc23))
anova(lm(bmi_change ~ high_sep + green_inc23), lm(bmi_change ~ high_sep + green_inc23 + high_sep*green_inc23))

AIC(lm(bmi_change ~ high_sep + green_inc23)); AIC(lm(bmi_change ~ high_sep + green_inc23 + high_sep*green_inc23))
BIC(lm(bmi_change ~ high_sep + green_inc23)); BIC(lm(bmi_change ~ high_sep + green_inc23 + high_sep*green_inc23))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_change_noInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (0 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_noInt)
sd(bmi_change); sd(bmi_change_noInt)

# Very small interaction
bmi_change_vSmallInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (0.5 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_vSmallInt)
sd(bmi_change); sd(bmi_change_vSmallInt)

# Small interaction
bmi_change_smallInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (1 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_smallInt)
sd(bmi_change); sd(bmi_change_smallInt)

# Large interaction
bmi_change_largeInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (3 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_largeInt)
sd(bmi_change); sd(bmi_change_largeInt)

# very large interaction
bmi_change_vLargeInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (4 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_vLargeInt)
sd(bmi_change); sd(bmi_change_vLargeInt)



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

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green3 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green3 here.
bmi <- 25 + (-4 * high_sep) + (-2 * green3) + (2 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi)

## descriptive stats split by possible combinations of SEP and green1 to check simulation worked
summary(bmi[high_sep == 1 & green3 == 1]) # High SEP and recent access to green space - Lowest BMI
summary(bmi[high_sep == 1 & green3 == 0]) # High SEP and no recent access to green space - Lowest BMI
summary(bmi[high_sep == 0 & green3 == 1]) # Low SEP and recent access to green space - Middle BMI
summary(bmi[high_sep == 0 & green3 == 0]) # Low SEP and no recent access to green space - Highest BMI

# Check the 'true' model, which is SEP as confounder, interaction with green 3, and main effect of green 3. green1 and green2 should also be pretty much null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))

AIC(lm(bmi ~ high_sep + green1 + green2 + green3)); AIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
BIC(lm(bmi ~ high_sep + green1 + green2 + green3)); BIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_noInt <- 25 + (-4 * high_sep) + (-2 * green3) + (0 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_noInt)
sd(bmi); sd(bmi_noInt)

# Very small interaction
bmi_vSmallInt <- 25 + (-4 * high_sep) + (-2 * green3) + (0.5 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_vSmallInt)
sd(bmi); sd(bmi_vSmallInt)

# Small interaction
bmi_smallInt <- 25 + (-4 * high_sep) + (-2 * green3) + (1 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_smallInt)
sd(bmi); sd(bmi_smallInt)

# Large interaction
bmi_largeInt <- 25 + (-4 * high_sep) + (-2 * green3) + (3 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_largeInt)
sd(bmi); sd(bmi_largeInt)

# very large interaction
bmi_vLargeInt <- 25 + (-4 * high_sep) + (-2 * green3) + (4 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_vLargeInt)
sd(bmi); sd(bmi_vLargeInt)



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


## Also want to create BMI outcome using different life course interaction to see if results differ. First, do interaction with accumulation (using the same model spces as above, so that access to green space only lowers BMI if from low SEP, but no difference if high SEP). As accumulation is continuous (from 0 to 3), will take two-thirds off the previous binary critical period coefficients so are approximately on the same scale
bmi_accum <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (0.67 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum)

# Check the 'true' model, which is SEP as confounder, interaction with accumulation, and main effect of accumulation. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi_accum ~ high_sep + accumulation))
summary(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))
anova(lm(bmi_accum ~ high_sep + accumulation), lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))

AIC(lm(bmi_accum ~ high_sep + accumulation)); AIC(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))
BIC(lm(bmi_accum ~ high_sep + accumulation)); BIC(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_accum_noInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (0 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_noInt)
sd(bmi_accum); sd(bmi_accum_noInt)

# Very small interaction
bmi_accum_vSmallInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (0.167 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_vSmallInt)
sd(bmi_accum); sd(bmi_accum_vSmallInt)

# Small interaction
bmi_accum_smallInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (0.33 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_smallInt)
sd(bmi_accum); sd(bmi_accum_smallInt)

# Large interaction
bmi_accum_largeInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (1 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_largeInt)
sd(bmi_accum); sd(bmi_accum_largeInt)

# very large interaction
bmi_accum_vLargeInt <- 25 + (-4 * high_sep) + (-0.67 * accumulation) + (1.33 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_vLargeInt)
sd(bmi_accum); sd(bmi_accum_vLargeInt)


## Next, construct a BMI outcome with one of the 'change' hypotheses as the main effect and interacting variable. Will say that an increase in green space from time 2 to time 3 lowers BMI, but only in those from lower SEP. As change is binary, will use the same parameters as from the critical period interaction model
bmi_change <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (2 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change)

# Check the 'true' model, which is SEP as confounder, interaction with green_inc23, and main effects of green_inc23. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi_change ~ high_sep + green_inc23))
summary(lm(bmi_change ~ high_sep + green_inc23 + high_sep*green_inc23))
anova(lm(bmi_change ~ high_sep + green_inc23), lm(bmi_change ~ high_sep + green_inc23 + high_sep*green_inc23))

AIC(lm(bmi_change ~ high_sep + green_inc23)); AIC(lm(bmi_change ~ high_sep + green_inc23 + high_sep*green_inc23))
BIC(lm(bmi_change ~ high_sep + green_inc23)); BIC(lm(bmi_change ~ high_sep + green_inc23 + high_sep*green_inc23))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_change_noInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (0 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_noInt)
sd(bmi_change); sd(bmi_change_noInt)

# Very small interaction
bmi_change_vSmallInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (0.5 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_vSmallInt)
sd(bmi_change); sd(bmi_change_vSmallInt)

# Small interaction
bmi_change_smallInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (1 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_smallInt)
sd(bmi_change); sd(bmi_change_smallInt)

# Large interaction
bmi_change_largeInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (3 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_largeInt)
sd(bmi_change); sd(bmi_change_largeInt)

# very large interaction
bmi_change_vLargeInt <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (4 * high_sep * green_inc23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_vLargeInt)
sd(bmi_change); sd(bmi_change_vLargeInt)



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

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green3 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green3 here. Have chosen 'green3' parameter to be 0.02 BMI per unit increase in green space, as the green3 standard deviation is ~50, and two times this should cover most of the variation in green space distance, making results broadly comparable to the binary green space effect of 2 BMI units (see: Gelman, A. (2008). Scaling regression inputs by dividing by two standard deviations. Statistics in medicine, 27(15), 2865-2873.)
bmi <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.02 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi)
sd(bmi)

# Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be pretty much null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))

AIC(lm(bmi ~ high_sep + green1 + green2 + green3)); AIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
BIC(lm(bmi ~ high_sep + green1 + green2 + green3)); BIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_noInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_noInt)
sd(bmi); sd(bmi_noInt)

# Very small interaction
bmi_vSmallInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.005 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_vSmallInt)
sd(bmi); sd(bmi_vSmallInt)

# Small interaction
bmi_smallInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.01 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_smallInt)
sd(bmi); sd(bmi_smallInt)

# Large interaction
bmi_largeInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.03 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_largeInt)
sd(bmi); sd(bmi_largeInt)

# very large interaction
bmi_vLargeInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.04 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_vLargeInt)
sd(bmi); sd(bmi_vLargeInt)


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
accumulation <- (green1 + green2 + green3) / 3

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



## Also want to create BMI outcome using different life course interaction to see if results differ. First, do interaction with accumulation (using the same model specs as above, so that access to green space only lowers BMI if from low SEP, but no difference if high SEP).
bmi_accum <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.02 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum)
sd(bmi_accum)

# Check the 'true' model, which is SEP as confounder, interaction with accumulation, and main effect of accumulation. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi_accum ~ high_sep + accumulation))
summary(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))
anova(lm(bmi_accum ~ high_sep + accumulation), lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))

AIC(lm(bmi_accum ~ high_sep + accumulation)); AIC(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))
BIC(lm(bmi_accum ~ high_sep + accumulation)); BIC(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_accum_noInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_noInt)
sd(bmi_accum); sd(bmi_accum_noInt)

# Very small interaction
bmi_accum_vSmallInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.005 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_vSmallInt)
sd(bmi_accum); sd(bmi_accum_vSmallInt)

# Small interaction
bmi_accum_smallInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.01 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_smallInt)
sd(bmi_accum); sd(bmi_accum_smallInt)

# Large interaction
bmi_accum_largeInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.03 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_largeInt)
sd(bmi_accum); sd(bmi_accum_largeInt)

# very large interaction
bmi_accum_vLargeInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.04 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_vLargeInt)
sd(bmi_accum); sd(bmi_accum_vLargeInt)


## Next, construct a BMI outcome with one of the 'change' hypotheses as the main effect and interacting variable. Will say that an increase in green space from time 2 to time 3 lowers BMI, but only in those from lower SEP.
bmi_change <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.02 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change)
sd(bmi_change)

# Check the 'true' model, which is SEP as confounder, interaction with green_inc23, and main effects of green_inc23. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi_change ~ high_sep + green_ch23))
summary(lm(bmi_change ~ high_sep + green_ch23 + high_sep*green_ch23))
anova(lm(bmi_change ~ high_sep + green_ch23), lm(bmi_change ~ high_sep + green_ch23 + high_sep*green_ch23))

AIC(lm(bmi_change ~ high_sep + green_ch23)); AIC(lm(bmi_change ~ high_sep + green_ch23 + high_sep*green_ch23))
BIC(lm(bmi_change ~ high_sep + green_ch23)); BIC(lm(bmi_change ~ high_sep + green_ch23 + high_sep*green_ch23))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_change_noInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_noInt)
sd(bmi_change); sd(bmi_change_noInt)

# Very small interaction
bmi_change_vSmallInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.005 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_vSmallInt)
sd(bmi_change); sd(bmi_change_vSmallInt)

# Small interaction
bmi_change_smallInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.01 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_smallInt)
sd(bmi_change); sd(bmi_change_smallInt)

# Large interaction
bmi_change_largeInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.03 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_largeInt)
sd(bmi_change); sd(bmi_change_largeInt)

# very large interaction
bmi_change_vLargeInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.04 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_vLargeInt)
sd(bmi_change); sd(bmi_change_vLargeInt)



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
accumulation <- (green1 + green2 + green3) / 3
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

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green3 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green3 here.
bmi <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.02 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi)
sd(bmi)

# Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be pretty much null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))

AIC(lm(bmi ~ high_sep + green1 + green2 + green3)); AIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
BIC(lm(bmi ~ high_sep + green1 + green2 + green3)); BIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_noInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_noInt)
sd(bmi); sd(bmi_noInt)

# Very small interaction
bmi_vSmallInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.005 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_vSmallInt)
sd(bmi); sd(bmi_vSmallInt)

# Small interaction
bmi_smallInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.01 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_smallInt)
sd(bmi); sd(bmi_smallInt)

# Large interaction
bmi_largeInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.03 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_largeInt)
sd(bmi); sd(bmi_largeInt)

# very large interaction
bmi_vLargeInt <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.04 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi); summary(bmi_vLargeInt)
sd(bmi); sd(bmi_vLargeInt)


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
accumulation <- (green1 + green2 + green3) / 3

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


## Also want to create BMI outcome using different life course interaction to see if results differ. First, do interaction with accumulation (using the same model specs as above, so that access to green space only lowers BMI if from low SEP, but no difference if high SEP).
bmi_accum <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.02 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum)

# Check the 'true' model, which is SEP as confounder, interaction with accumulation, and main effect of accumulation. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi_accum ~ high_sep + accumulation))
summary(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))
anova(lm(bmi_accum ~ high_sep + accumulation), lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))

AIC(lm(bmi_accum ~ high_sep + accumulation)); AIC(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))
BIC(lm(bmi_accum ~ high_sep + accumulation)); BIC(lm(bmi_accum ~ high_sep + accumulation + high_sep*accumulation))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_accum_noInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_noInt)
sd(bmi_accum); sd(bmi_accum_noInt)

# Very small interaction
bmi_accum_vSmallInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.005 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_vSmallInt)
sd(bmi_accum); sd(bmi_accum_vSmallInt)

# Small interaction
bmi_accum_smallInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.01 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_smallInt)
sd(bmi_accum); sd(bmi_accum_smallInt)

# Large interaction
bmi_accum_largeInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.03 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_largeInt)
sd(bmi_accum); sd(bmi_accum_largeInt)

# very large interaction
bmi_accum_vLargeInt <- 28 + (-4 * high_sep) + (-0.02 * accumulation) + (0.04 * high_sep * accumulation) + rnorm(n, 0, 3)
summary(bmi_accum); summary(bmi_accum_vLargeInt)
sd(bmi_accum); sd(bmi_accum_vLargeInt)



## Next, construct a BMI outcome with one of the 'change' hypotheses as the main effect and interacting variable. Will say that an increase in green space from time 2 to time 3 lowers BMI, but only in those from lower SEP.
bmi_change <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.02 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change)

# Check the 'true' model, which is SEP as confounder, interaction with green_inc23, and main effects of green_inc23. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi_change ~ high_sep + green_ch23))
summary(lm(bmi_change ~ high_sep + green_ch23 + high_sep*green_ch23))
anova(lm(bmi_change ~ high_sep + green_ch23), lm(bmi_change ~ high_sep + green_ch23 + high_sep*green_ch23))

AIC(lm(bmi_change ~ high_sep + green_ch23)); AIC(lm(bmi_change ~ high_sep + green_ch23 + high_sep*green_ch23))
BIC(lm(bmi_change ~ high_sep + green_ch23)); BIC(lm(bmi_change ~ high_sep + green_ch23 + high_sep*green_ch23))


## Test different interaction strengths to ensure data looks sensible

# No interaction
bmi_change_noInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_noInt)
sd(bmi_change); sd(bmi_change_noInt)

# Very small interaction
bmi_change_vSmallInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.005 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_vSmallInt)
sd(bmi_change); sd(bmi_change_vSmallInt)

# Small interaction
bmi_change_smallInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.01 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_smallInt)
sd(bmi_change); sd(bmi_change_smallInt)

# Large interaction
bmi_change_largeInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.03 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_largeInt)
sd(bmi_change); sd(bmi_change_largeInt)

# very large interaction
bmi_change_vLargeInt <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.04 * high_sep * green_ch23) + rnorm(n, 0, 3)
summary(bmi_change); summary(bmi_change_vLargeInt)
sd(bmi_change); sd(bmi_change_vLargeInt)



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
accumulation <- (green1 + green2 + green3) / 3
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
  AIC_res_temp <- rep(NA, n_sims)
  AIC_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
  AIC_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
  AIC_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
  AIC_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
  BIC_res_temp <- rep(NA, n_sims)
  BIC_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
  BIC_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
  BIC_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
  BIC_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
  CV_1SE_res_temp <- rep(NA, n_sims)
  CV_1SE_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
  CV_1SE_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
  CV_1SE_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
  CV_1SE_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
  CV_minMSE_res_temp <- rep(NA, n_sims)
  CV_minMSE_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
  CV_minMSE_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
  CV_minMSE_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
  CV_minMSE_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
  
  for (i in 1:n_sims) {
    
    if (Output == TRUE) {
      print(paste0("Simulation number: ", i, " / ", n_sims))
      print(paste0("Sample size = ", sampleSize, "; Exposure = ", Exposure, "; Centered = ", Centered, 
                   "; Collinear = ", Collinear, "; Outcome = ", Outcome))
      print("")
    }
    
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
    }
    
    ## If exposure is continuous and collinearity low
    if (Exposure == "Cont" & Collinear == "Low") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40) # First green space distance exposure in pregnancy
      green2 <- 240 + (-50 * high_sep) + (0.3 * green1) + rnorm(n, 0, 40) # Second green space distance exposure
      green3 <- 240 + (-50 * high_sep) + (0.3 * green2) + rnorm(n, 0, 40) # Third green space distance exposure
    }
    
    ## If exposure is continuous and collinearity high
    if (Exposure == "Cont" & Collinear == "High") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40) # First green space distance exposure in pregnancy
      green2 <- 50 + (-50 * high_sep) + (0.9 * green1) + rnorm(n, 0, 20) # Second green space distance exposure
      green3 <- 50 + (-50 * high_sep) + (0.9 * green2) + rnorm(n, 0, 20) # Third green space distance exposure
    }

    
    ## Encode the life course hypotheses
    
    # Critical periods (same for binary and continuous exposures)
    crit1 <- green1 # Critical period at first time point only
    int1 <- crit1 * high_sep # Interaction between SEP and first time point
    
    crit2 <- green2 # Critical period at second time point only
    int2 <- crit2 * high_sep # Interaction between SEP and second time point
    
    crit3 <- green3 # Critical period at third time point only
    int3 <- crit3 * high_sep # Interaction between SEP and third time point
    
    # Accumulation (different for binary and continuous exposures)
    if (Exposure == "Binary") {
      accumulation <- green1 + green2 + green3 # Linear accumulation of all exposures
      int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
    }
    
    if (Exposure == "Cont") {
      accumulation <- (green1 + green2 + green3) / 3 # Linear accumulation of all exposures
      int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
    }
    
    # Change (different for binary and continuous exposures)
    if (Exposure == "Binary") {
      green_inc12 <- (1 - green1) * green2 # Increase from time 1 to time 2
      green_inc12_int <- green_inc12 * high_sep # Increase from time 1 to time 2, with an interaction with SEP
      
      green_dec12 <- (1 - green2) * green1 # Decrease from time 1 to time 2
      green_dec12_int <- green_dec12 * high_sep # Decrease from time 1 to time 2, with an interaction with SEP
      
      green_inc23 <- (1 - green2) * green3 # Increase from time 2 to time 3
      green_inc23_int <- green_inc23 * high_sep # Increase from time 2 to time 3, with an interaction with SEP
      
      green_dec23 <- (1 - green3) * green2 # Decrease from time 2 to time 3
      green_dec23_int <- green_dec23 * high_sep # Decrease from time 2 to time 3, with an interaction with SEP
    }
    
    if (Exposure == "Cont") {
      green_ch12 <- green2 - green1 # Change from time 1 to time 2
      green_ch12_int <- green_ch12 * high_sep # Change from time 1 to time 2, with an interaction with SEP
      
      green_ch23 <- green3 - green2 # Change from time 2 to time 3
      green_ch23_int <- green_ch23 * high_sep # Change from time 2 to time 3, with an interaction with SEP
    }
    
    
    ## Create the outcomes (differs depending on whether exposures are binary or continuous)
    if (Exposure == "Binary") {
      bmi <- 25 + (-4 * high_sep) + (-2 * green3) + (2 * high_sep * green3) + rnorm(n, 0, 3)
    }
    
    if (Exposure == "Cont") {
      bmi <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.02 * high_sep * green3) + rnorm(n, 0, 3)
    }
    
    # If outcome is binary
    if (Outcome == "Binary") {
      overweight <- ifelse(bmi > 25, 1, 0) # Code BMI to binary 'overweight' variable if > 25
    }
    
    
    ## Center the exposures, if specified
    if (Centered == "Yes") {
      crit1 <- crit1 - mean(crit1) # Center critical period at time 1
      int1 <- crit1 * high_sep # Interaction between SEP and first time point
      
      crit2 <- crit2 - mean(crit2) # Center critical period at time 2
      int2 <- crit2 * high_sep # Interaction between SEP and second time point
      
      crit3 <- crit3 - mean(crit3) # Center critical period at time 3
      int3 <- crit3 * high_sep # Interaction between SEP and third time point
      
      accumulation <- accumulation - mean(accumulation) # Center the accumulation variable
      int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
      
      # Change (different for binary and continuous exposures)
      if (Exposure == "Binary") {
        green_inc12 <- green_inc12 - mean(green_inc12) # Center the increase from time 1 to time 2 variable
        green_inc12_int <- green_inc12 * high_sep # Increase from time 1 to time 2, with an interaction with SEP
        
        green_dec12 <- green_dec12 - mean(green_dec12) # Center the decrease from time 1 to time 2 variable
        green_dec12_int <- green_dec12 * high_sep # Decrease from time 1 to time 2, with an interaction with SEP
        
        green_inc23 <- green_inc23 - mean(green_inc23) # Center the increase from time 2 to time 3 variable
        green_inc23_int <- green_inc23 * high_sep # Increase from time 2 to time 3, with an interaction with SEP
        
        green_dec23 <- green_dec23 - mean(green_dec23) # Center the decrease from time 2 to time 3 variable
        green_dec23_int <- green_dec23 * high_sep # Decrease from time 2 to time 3, with an interaction with SEP
      }
      if (Exposure == "Cont") {
        green_ch12 <- green_ch12 - mean(green_ch12) # Center the change from time 1 to time 2 variable
        green_ch12_int <- green_ch12 * high_sep # Change from time 1 to time 2, with an interaction with SEP
        
        green_ch23 <- green_ch23 - mean(green_ch23) # Center the change from time 2 to time 3 variable
        green_ch23_int <- green_ch23 * high_sep # Change from time 2 to time 3, with an interaction with SEP
      }
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
    
    ## Loop over each step of the lasso, running a standard LM/GLM model (depending on outcome) whenever variables change. Store the AIC and BIC values and select the model with the lowest values
    
    # Make a dataframe of all the models identified by the lasso
    old_covars <- ""
    old_deviance <- 0
    old_varNum <- 0
    old_lambda <- NA
    df <- data.frame(matrix(ncol = 8, nrow = 0))
    #df
    
    for (j in 1:length(mod$df)) {
      #print(j)
      new_covars <- attributes(which(mod$beta[, j] != 0))$names
      new_deviance <- mod$dev.ratio[j]
      new_varNum <- mod$df[j]
      new_lambda <- mod$lambda[j]
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
    
    # Make a var to show number of steps where variables added, and rename the high_sep variable to blank (as is included by default)
    df$steps <- 1:nrow(df)
    df$Variables[df$steps == 1] <- ""
    #df
    
    
    ## Now run all combinations of the model variables in a standard LM/GLM, and store AIC and BIC values
    for (j in 1:nrow(df)) {
      
      vars_temp <- strsplit(df$model_vars[j], " ")[[1]] # Split the variables at each stage of the lasso
      x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp] # Matrix of the covariates at each step of lasso
      
      # If continuous outcome
      if (Outcome == "Cont") {
        mod_ic <- lm(bmi ~ x_hypos_new) # Run the model
        #print(summary(mod_ic))
      }
      
      # If binary outcome
      if (Outcome == "Binary") {
        mod_ic <- glm(overweight ~ x_hypos_new, family = "binomial") # Run the model
        #print(summary(mod_ic))
      }
      
      # Store the AIC values
      df$aic[j] <- AIC(mod_ic)
      df$bic[j] <- BIC(mod_ic)
      
    }
    
    # Select the models with the lowest AIC and BIC values and store these
    #df
    
    df$model_vars[which.min(df$aic)]
    df$model_vars[which.min(df$bic)]
    
    
    # Print the covariates in the best-fitting models if want to print this
    if (Output == TRUE) {
      print(paste0("Covariates in best-fitting AIC model:"))
      print(df$model_vars[which.min(df$aic)])
      print("")
      
      print(paste0("Covariates in best-fitting BIC model:"))
      print(df$model_vars[which.min(df$bic)])
      print("")
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
    vars_split_aic <- strsplit(df$model_vars[which.min(df$aic)], " ")[[1]] # Split the variables
    AIC_res_temp[i] <- ifelse(setequal(vars_split_aic, target_covars) == TRUE, 1, 0)
    
    vars_split_bic <- strsplit(df$model_vars[which.min(df$bic)], " ")[[1]] # Split the variables
    BIC_res_temp[i] <- ifelse(setequal(vars_split_bic, target_covars) == TRUE, 1, 0)
    
    # Store if 'crit3' in final model
    AIC_res_temp_crit3[i] <- ifelse("crit3" %in% vars_split_aic == TRUE, 1, 0)
    BIC_res_temp_crit3[i] <- ifelse("crit3" %in% vars_split_bic == TRUE, 1, 0)
    
    # Store if 'int3' in final model
    AIC_res_temp_int3[i] <- ifelse("int3" %in% vars_split_aic == TRUE, 1, 0)
    BIC_res_temp_int3[i] <- ifelse("int3" %in% vars_split_bic == TRUE, 1, 0)
    
    # Store if 'crit3' and 'int3' both in final model
    AIC_res_temp_crit3int3[i] <- ifelse("crit3" %in% vars_split_aic == TRUE & "int3" %in% vars_split_aic == TRUE, 1, 0)
    BIC_res_temp_crit3int3[i] <- ifelse("crit3" %in% vars_split_bic == TRUE & "int3" %in% vars_split_bic == TRUE, 1, 0)
    
    # Store if 'crit3' and 'int3' both in final model, plus other variables
    AIC_res_temp_crit3int3extra[i] <- ifelse(AIC_res_temp_crit3int3[i] == 1 & AIC_res_temp[i] == 0, 1, 0)
    BIC_res_temp_crit3int3extra[i] <- ifelse(BIC_res_temp_crit3int3[i] == 1 & BIC_res_temp[i] == 0, 1, 0)
    
    
    ### Next, want to summarise the 1SE and minimum MSE cross-validated lasso model, and see whether that corresponds to the correct model or not
    
    # If outcome is continuous
    if (Outcome == "Cont") {
      mod.cv <- cv.glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    # If outcome is binary
    if (Outcome == "Binary") {
      mod.cv <- cv.glmnet(x_hypos, overweight, alpha = 1, family = "binomial", 
                          penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    
    ## Start with the 1 SE model
    # Extract the non-zero covariates into a vector
    cv_1SE_covars_temp <- as.matrix(coef(mod.cv, s = mod.cv$lambda.1se))
    cv_1SE_covars_temp <- as.matrix(cv_1SE_covars_temp[!rownames(cv_1SE_covars_temp) %in% "(Intercept)", ]) # Drop the intercept
    cv_1SE_covars_temp <- cv_1SE_covars_temp[cv_1SE_covars_temp != 0, ] # Drop all zero coefficients
    
    cv_1SE_covars <- attributes(cv_1SE_covars_temp)$names # Store all the non-zero hypotheses
    
    # Print the covariates in the best-fitting cross-validated model
    if (Output == TRUE) {
      print(paste0("CV 1 SE covariates:"))
      print(cv_1SE_covars)
      print("")
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
    CV_1SE_res_temp[i] <- ifelse(setequal(cv_1SE_covars, target_covars) == TRUE, 1, 0)
    
    # Store if 'crit3', 'int3', both, and both with additional variables in final model
    CV_1SE_res_temp_crit3[i] <- ifelse("crit3" %in% cv_1SE_covars == TRUE, 1, 0)
    CV_1SE_res_temp_int3[i] <- ifelse("int3" %in% cv_1SE_covars == TRUE, 1, 0)
    CV_1SE_res_temp_crit3int3[i] <- ifelse("crit3" %in% cv_1SE_covars == TRUE & "int3" %in% cv_1SE_covars == TRUE, 1, 0)
    CV_1SE_res_temp_crit3int3extra[i] <- ifelse(CV_1SE_res_temp_crit3int3[i] == 1 & CV_1SE_res_temp[i] == 0, 1, 0)
    
    
    ## Next to the minimum MSE model
    # Extract the non-zero covariates into a vector
    cv_minMSE_covars_temp <- as.matrix(coef(mod.cv, s = mod.cv$lambda.min))
    cv_minMSE_covars_temp <- as.matrix(cv_minMSE_covars_temp[!rownames(cv_minMSE_covars_temp) %in% "(Intercept)", ]) # Drop the intercept
    cv_minMSE_covars_temp <- cv_minMSE_covars_temp[cv_minMSE_covars_temp != 0, ] # Drop all zero coefficients
    
    cv_minMSE_covars <- attributes(cv_minMSE_covars_temp)$names # Store all the non-zero hypotheses
    
    # Print the covariates in the best-fitting cross-validated model
    if (Output == TRUE) {
      print(paste0("CV min MSE covariates:"))
      print(cv_minMSE_covars)
      print("")
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
    CV_minMSE_res_temp[i] <- ifelse(setequal(cv_minMSE_covars, target_covars) == TRUE, 1, 0)
    
    # Store if 'crit3', 'int3', both, and both with additional variables in final model
    CV_minMSE_res_temp_crit3[i] <- ifelse("crit3" %in% cv_minMSE_covars == TRUE, 1, 0)
    CV_minMSE_res_temp_int3[i] <- ifelse("int3" %in% cv_minMSE_covars == TRUE, 1, 0)
    CV_minMSE_res_temp_crit3int3[i] <- ifelse("crit3" %in% cv_minMSE_covars == TRUE & "int3" %in% cv_minMSE_covars == TRUE, 1, 0)
    CV_minMSE_res_temp_crit3int3extra[i] <- ifelse(CV_minMSE_res_temp_crit3int3[i] == 1 & CV_minMSE_res_temp[i] == 0, 1, 0)
    
  }
  
  # Store the summaries of these results to transfer to the main results table
  res <- data.frame(AIC_propcorrect = round(sum(AIC_res_temp) / n_sims * 100, 2),
                    AIC_crit3correct = round(sum(AIC_res_temp_crit3) / n_sims * 100, 2),
                    AIC_int3correct = round(sum(AIC_res_temp_int3) / n_sims * 100, 2),
                    AIC_crit3int3correct = round(sum(AIC_res_temp_crit3int3) / n_sims * 100, 2),
                    AIC_crit3int3extra = round(sum(AIC_res_temp_crit3int3extra) / n_sims * 100, 2),
                    BIC_propcorrect = round(sum(BIC_res_temp) / n_sims * 100, 2),
                    BIC_crit3correct = round(sum(BIC_res_temp_crit3) / n_sims * 100, 2),
                    BIC_int3correct = round(sum(BIC_res_temp_int3) / n_sims * 100, 2),
                    BIC_crit3int3correct = round(sum(BIC_res_temp_crit3int3) / n_sims * 100, 2),
                    BIC_crit3int3extra = round(sum(BIC_res_temp_crit3int3extra) / n_sims * 100, 2),
                    CV_1SE_propcorrect = round(sum(CV_1SE_res_temp) / n_sims * 100, 2),
                    CV_1SE_crit3correct = round(sum(CV_1SE_res_temp_crit3) / n_sims * 100, 2),
                    CV_1SE_int3correct = round(sum(CV_1SE_res_temp_int3) / n_sims * 100, 2),
                    CV_1SE_crit3int3correct = round(sum(CV_1SE_res_temp_crit3int3) / n_sims * 100, 2),
                    CV_1SE_crit3int3extra = round(sum(CV_1SE_res_temp_crit3int3extra) / n_sims * 100, 2),
                    CV_minMSE_propcorrect = round(sum(CV_minMSE_res_temp) / n_sims * 100, 2),
                    CV_minMSE_crit3correct = round(sum(CV_minMSE_res_temp_crit3) / n_sims * 100, 2),
                    CV_minMSE_int3correct = round(sum(CV_minMSE_res_temp_int3) / n_sims * 100, 2),
                    CV_minMSE_crit3int3correct = round(sum(CV_minMSE_res_temp_crit3int3) / n_sims * 100, 2),
                    CV_minMSE_crit3int3extra = round(sum(CV_minMSE_res_temp_crit3int3extra) / n_sims * 100, 2))
  return(res)
  
}


## Next, the number of simulations per combination of parameters (1,000), the target hypotheses we simulated are the true model, and set up a data frame to store the results in
n_sims <- 10
set.seed(9876)

target_covars <- c("high_sep", "crit3", "int3")

results <- data.frame(sampleSize = rep(c(1000, 10000), 16),
                      Exposure = rep(c(rep("Binary", 8), rep("Cont", 8)), 2),
                      Centered = rep(c("No", "No", "Yes", "Yes"), 8),
                      Collinear = rep(c(rep("Low", 4), rep("High", 4)), 4),
                      Outcome = c(rep("Cont", 16), rep("Binary", 16)),
                      AIC_propcorrect = rep(NA, 32),
                      AIC_crit3correct = rep(NA, 32),
                      AIC_int3correct = rep(NA, 32),
                      AIC_crit3int3correct = rep(NA, 32),
                      AIC_crit3int3extra = rep(NA, 32),
                      BIC_propcorrect = rep(NA, 32),
                      BIC_crit3correct = rep(NA, 32),
                      BIC_int3correct = rep(NA, 32),
                      BIC_crit3int3correct = rep(NA, 32),
                      BIC_crit3int3extra = rep(NA, 32),
                      CV_1SE_propcorrect = rep(NA, 32),
                      CV_1SE_crit3correct = rep(NA, 32),
                      CV_1SE_int3correct = rep(NA, 32),
                      CV_1SE_crit3int3correct = rep(NA, 32),
                      CV_1SE_crit3int3extra = rep(NA, 32),
                      CV_minMSE_propcorrect = rep(NA, 32),
                      CV_minMSE_crit3correct = rep(NA, 32),
                      CV_minMSE_int3correct = rep(NA, 32),
                      CV_minMSE_crit3int3correct = rep(NA, 32),
                      CV_minMSE_crit3int3extra = rep(NA, 32))

results


# Also calculate time taken to run script
start_time <- Sys.time()


### First simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 1
results[k, 6:25] <- res


### Second simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 2
results[k, 6:25] <- res


### Third simulation: Sample size = 1000; binary exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 3
results[k, 6:25] <- res


### Fourth simulation: Sample size = 10000; binary exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 4
results[k, 6:25] <- res


### Fifth simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 5
results[k, 6:25] <- res


### Sixth simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 6
results[k, 6:25] <- res


### Seventh simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 7
results[k, 6:25] <- res


### Eighth simulation: Sample size = 10000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 8
results[k, 6:25] <- res


### Ninth simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 9
results[k, 6:25] <- res


### Tenth simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 10
results[k, 6:25] <- res


### Eleventh simulation: Sample size = 1000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 11
results[k, 6:25] <- res


### Twelfth simulation: Sample size = 10000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 12
results[k, 6:25] <- res


### 13th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 13
results[k, 6:25] <- res


### 14th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 14
results[k, 6:25] <- res


### 15th simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 15
results[k, 6:25] <- res


### 16th simulation: Sample size = 10000; continuous exposure; centered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 16
results[k, 6:25] <- res


### 17th simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 17
results[k, 6:25] <- res


### 18th simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 18
results[k, 6:25] <- res


### 19th simulation: Sample size = 1000; binary exposure; centered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 19
results[k, 6:25] <- res


### 20th simulation: Sample size = 10000; binary exposure; centered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 20
results[k, 6:25] <- res


### 21st simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 21
results[k, 6:25] <- res


### 22nd simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 22
results[k, 6:25] <- res


### 23rd simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 23
results[k, 6:25] <- res


### 24th simulation: Sample size = 10000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 24
results[k, 6:25] <- res


### 25th simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 25
results[k, 6:25] <- res


### 26th simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 26
results[k, 6:25] <- res


### 27th simulation: Sample size = 1000; continuous exposure; centered; low collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 27
results[k, 6:25] <- res


### 28th simulation: Sample size = 10000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 28
results[k, 6:25] <- res


### 29th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 29
results[k, 6:25] <- res


### 30th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 30
results[k, 6:25] <- res


### 31st simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 31
results[k, 6:25] <- res


### 32nd simulation: Sample size = 10000; continuous exposure; centered; High collinearity; binary outcome
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 32
results[k, 6:25] <- res


# Time taken to run script (took about 40 mins for 100 simulations, so should be about 6 hours for the full 1,000 simulations)
end_time <- Sys.time()

end_time - start_time


## Save the results table
results

write.csv(results, file = "simulationResults_test.csv", row.names = FALSE)



### Testing parameter combination 2 (Sample size = 10000; binary exposure; uncentered; low collinearity; continuous outcome) to explore why it correctly detects crit3 and int3 in all models, but erroneously includes other variables as well

# Initiate a vector to save results from this simulation to (i.e., whether method identified correct model or not)
set.seed(1234)
n_sims <- 1

AIC_res_temp <- rep(NA, n_sims)
AIC_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
AIC_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
AIC_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
AIC_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
BIC_res_temp <- rep(NA, n_sims)
BIC_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
BIC_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
BIC_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
BIC_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
CV_1SE_res_temp <- rep(NA, n_sims)
CV_1SE_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
CV_1SE_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
CV_1SE_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
CV_1SE_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
CV_minMSE_res_temp <- rep(NA, n_sims)
CV_minMSE_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
CV_minMSE_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
CV_minMSE_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
CV_minMSE_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model

for (i in 1:n_sims) {
  
  print(paste0("Simulation number: ", i))
  
  ## Start by simulating the data
  n <- 10000
  
  high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
  green1_p <- plogis(log(0.6) + (log(3) * high_sep)) # First green space exposure in pregnancy
  green1 <- rbinom(n, 1, green1_p)
  green2_p <- plogis(log(0.3) + (log(3) * high_sep) + (log(3) * green1)) # Second green space exposure
  green2 <- rbinom(n, 1, green2_p)
  green3_p <- plogis(log(0.3) + (log(3) * high_sep) + (log(3) * green2)) # Third green space exposure
  green3 <- rbinom(n, 1, green3_p)
  bmi <- 25 + (-4 * high_sep) + (-2 * green3) + (2 * high_sep * green3) + rnorm(n, 0, 3) # Cont. BMI outcome
  
  ## Encode the life course hypotheses
  crit1 <- green1 # Critical period at first time point only
  #crit1 <- crit1 - mean(crit1) # Center this variable
  int1 <- crit1 * high_sep # Interaction between SEP and first time point
  crit2 <- green2 # Critical period at second time point only
  #crit2 <- crit2 - mean(crit2) # Center this variable
  int2 <- crit2 * high_sep # Interaction between SEP and second time point
  crit3 <- green3 # Critical period at third time point only
  #crit3 <- crit3 - mean(crit3) # Center this variable
  int3 <- crit3 * high_sep # Interaction between SEP and third time point
  accumulation <- green1 + green2 + green3 # Linear accumulation of all exposures
  #accumulation <- accumulation - mean(accumulation) # Center this variable
  int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
  green_inc12 <- (1 - green1) * green2 # Increase from time 1 to time 2
  #green_inc12 <- green_inc12 - mean(green_inc12) # Center this variable
  green_inc12_int <- green_inc12 * high_sep # Increase from time 1 to time 2, with an interaction with SEP
  green_dec12 <- (1 - green2) * green1 # Decrease from time 1 to time 2
  #green_dec12 <- green_dec12 - mean(green_dec12) # Center this variable
  green_dec12_int <- green_dec12 * high_sep # Decrease from time 1 to time 2, with an interaction with SEP
  green_inc23 <- (1 - green2) * green3 # Increase from time 2 to time 3
  #green_inc23 <- green_inc23 - mean(green_inc23) # Center this variable
  green_inc23_int <- green_inc23 * high_sep # Increase from time 2 to time 3, with an interaction with SEP
  green_dec23 <- (1 - green3) * green2 # Decrease from time 2 to time 3
  #green_dec23 <- green_dec23 - mean(green_dec23) # Center this variable
  green_dec23_int <- green_dec23 * high_sep # Decrease from time 2 to time 3, with an interaction with SEP
  
  ## Combine all these into one matrix, including SEP as a confounder
  x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_inc12, 
                   green_inc12_int, green_dec12, green_dec12_int, green_inc23, green_inc23_int, green_dec23,
                   green_dec23_int)
  
  ## Run the model using glmnet
  mod <- glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
  
  ## Loop over each step of the lasso, running a standard LM/GLM model (depending on outcome) whenever variables change. Store the AIC and BIC values and select the model with the lowest values
  
  # Make a dataframe of all the models identified by the lasso
  old_covars <- ""
  old_deviance <- 0
  old_varNum <- 0
  old_lambda <- NA
  df <- data.frame(matrix(ncol = 8, nrow = 0))
  #df
  
  for (j in 1:length(mod$df)) {
    #print(j)
    new_covars <- attributes(which(mod$beta[, j] != 0))$names
    new_deviance <- mod$dev.ratio[j]
    new_varNum <- mod$df[j]
    new_lambda <- mod$lambda[j]
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
  
  # Make a var to show number of steps where variables added, and rename the high_sep variable to blank (as is included by default)
  df$steps <- 1:nrow(df)
  df$Variables[df$steps == 1] <- ""
  #df
  
  
  ## Now run all combinations of the model variables in a standard LM/GLM, and store AIC and BIC values
  for (j in 1:nrow(df)) {
    
    vars_temp <- strsplit(df$model_vars[j], " ")[[1]] # Split the variables at each stage of the lasso
    x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp] # Matrix of the covariates at each step of lasso
    
    mod_ic <- lm(bmi ~ x_hypos_new) # Run the model
    #print(summary(mod_ic))

    # Store the AIC values
    df$aic[j] <- AIC(mod_ic)
    df$bic[j] <- BIC(mod_ic)
    
  }
  
  # Print the covariates in the best-fitting models
  print(paste0("Covariates in best-fitting AIC model:"))
  print(df$model_vars[which.min(df$aic)])
  print("")
    
  print(paste0("Covariates in best-fitting BIC model:"))
  print(df$model_vars[which.min(df$bic)])
  print("")
  
  ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
  vars_split_aic <- strsplit(df$model_vars[which.min(df$aic)], " ")[[1]] # Split the variables
  AIC_res_temp[i] <- ifelse(setequal(vars_split_aic, target_covars) == TRUE, 1, 0)
  
  vars_split_bic <- strsplit(df$model_vars[which.min(df$bic)], " ")[[1]] # Split the variables
  BIC_res_temp[i] <- ifelse(setequal(vars_split_bic, target_covars) == TRUE, 1, 0)
  
  # Store if 'crit3' in final model
  AIC_res_temp_crit3[i] <- ifelse("crit3" %in% vars_split_aic == TRUE, 1, 0)
  BIC_res_temp_crit3[i] <- ifelse("crit3" %in% vars_split_bic == TRUE, 1, 0)
  
  # Store if 'int3' in final model
  AIC_res_temp_int3[i] <- ifelse("int3" %in% vars_split_aic == TRUE, 1, 0)
  BIC_res_temp_int3[i] <- ifelse("int3" %in% vars_split_bic == TRUE, 1, 0)
  
  # Store if 'crit3' and 'int3' both in final model
  AIC_res_temp_crit3int3[i] <- ifelse("crit3" %in% vars_split_aic == TRUE & "int3" %in% vars_split_aic == TRUE, 1, 0)
  BIC_res_temp_crit3int3[i] <- ifelse("crit3" %in% vars_split_bic == TRUE & "int3" %in% vars_split_bic == TRUE, 1, 0)
  
  # Store if 'crit3' and 'int3' both in final model, plus other variables
  AIC_res_temp_crit3int3extra[i] <- ifelse(AIC_res_temp_crit3int3[i] == 1 & AIC_res_temp[i] == 0, 1, 0)
  BIC_res_temp_crit3int3extra[i] <- ifelse(BIC_res_temp_crit3int3[i] == 1 & BIC_res_temp[i] == 0, 1, 0)
  
  
  ### Next, want to summarise the 1SE and minimum MSE cross-validated lasso model, and see whether that corresponds to the correct model or not
  mod.cv <- cv.glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
  
  ## Start with the 1 SE model
  # Extract the non-zero covariates into a vector
  cv_1SE_covars_temp <- as.matrix(coef(mod.cv, s = mod.cv$lambda.1se))
  cv_1SE_covars_temp <- as.matrix(cv_1SE_covars_temp[!rownames(cv_1SE_covars_temp) %in% "(Intercept)", ]) # Drop the intercept
  cv_1SE_covars_temp <- cv_1SE_covars_temp[cv_1SE_covars_temp != 0, ] # Drop all zero coefficients
  
  cv_1SE_covars <- attributes(cv_1SE_covars_temp)$names # Store all the non-zero hypotheses
  
  # Print the covariates in the best-fitting cross-validated model
  print(paste0("CV 1 SE covariates:"))
  print(cv_1SE_covars)
  print("")
  
  ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
  CV_1SE_res_temp[i] <- ifelse(setequal(cv_1SE_covars, target_covars) == TRUE, 1, 0)
  
  # Store if 'crit3', 'int3', both, and both with additional variables in final model
  CV_1SE_res_temp_crit3[i] <- ifelse("crit3" %in% cv_1SE_covars == TRUE, 1, 0)
  CV_1SE_res_temp_int3[i] <- ifelse("int3" %in% cv_1SE_covars == TRUE, 1, 0)
  CV_1SE_res_temp_crit3int3[i] <- ifelse("crit3" %in% cv_1SE_covars == TRUE & "int3" %in% cv_1SE_covars == TRUE, 1, 0)
  CV_1SE_res_temp_crit3int3extra[i] <- ifelse(CV_1SE_res_temp_crit3int3[i] == 1 & CV_1SE_res_temp[i] == 0, 1, 0)
  
  
  ## Next to the minimum MSE model
  # Extract the non-zero covariates into a vector
  cv_minMSE_covars_temp <- as.matrix(coef(mod.cv, s = mod.cv$lambda.min))
  cv_minMSE_covars_temp <- as.matrix(cv_minMSE_covars_temp[!rownames(cv_minMSE_covars_temp) %in% "(Intercept)", ]) # Drop the intercept
  cv_minMSE_covars_temp <- cv_minMSE_covars_temp[cv_minMSE_covars_temp != 0, ] # Drop all zero coefficients
  
  cv_minMSE_covars <- attributes(cv_minMSE_covars_temp)$names # Store all the non-zero hypotheses
  
  # Print the covariates in the best-fitting cross-validated model
  print(paste0("CV min MSE covariates:"))
  print(cv_minMSE_covars)
  print("")
  
  ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
  CV_minMSE_res_temp[i] <- ifelse(setequal(cv_minMSE_covars, target_covars) == TRUE, 1, 0)
  
  # Store if 'crit3', 'int3', both, and both with additional variables in final model
  CV_minMSE_res_temp_crit3[i] <- ifelse("crit3" %in% cv_minMSE_covars == TRUE, 1, 0)
  CV_minMSE_res_temp_int3[i] <- ifelse("int3" %in% cv_minMSE_covars == TRUE, 1, 0)
  CV_minMSE_res_temp_crit3int3[i] <- ifelse("crit3" %in% cv_minMSE_covars == TRUE & "int3" %in% cv_minMSE_covars == TRUE, 1, 0)
  CV_minMSE_res_temp_crit3int3extra[i] <- ifelse(CV_minMSE_res_temp_crit3int3[i] == 1 & CV_minMSE_res_temp[i] == 0, 1, 0)
  
}


### Look at each step of the lasso
df

## In this case, the first two models added are 'crit3' and 'int3', but there are some models with slightly lower AIC values which contain additional variables. It is also possible that in some cases the final model also includes additional variables because said variable gets added after crit3 but before int3 (where said variable *is* associated with the outcome), so does not get removed before 'int3' gets added to the model; but once 'int3' gets included, there is no association between said variable and the outcome. but gets included in the final model regardless.



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
  AIC_res_temp <- rep(NA, n_sims)
  AIC_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
  AIC_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
  AIC_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
  AIC_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
  BIC_res_temp <- rep(NA, n_sims)
  BIC_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
  BIC_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
  BIC_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
  BIC_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
  CV_1SE_res_temp <- rep(NA, n_sims)
  CV_1SE_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
  CV_1SE_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
  CV_1SE_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
  CV_1SE_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
  CV_minMSE_res_temp <- rep(NA, n_sims)
  CV_minMSE_res_temp_crit3 <- rep(NA, n_sims) # Store if 'crit3' main effect in final model
  CV_minMSE_res_temp_int3 <- rep(NA, n_sims) # Store if 'int3' interaction effect in final model
  CV_minMSE_res_temp_crit3int3 <- rep(NA, n_sims) # Store if 'crit3' and 'int3' in final model
  CV_minMSE_res_temp_crit3int3extra <- rep(NA, n_sims) # Store if 'crit3' and 'int3', plus extra vars, in final model
  
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
    }
    
    ## If exposure is continuous and collinearity low
    if (Exposure == "Cont" & Collinear == "Low") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40) # First green space distance exposure in pregnancy
      green2 <- 240 + (-50 * high_sep) + (0.3 * green1) + rnorm(n, 0, 40) # Second green space distance exposure
      green3 <- 240 + (-50 * high_sep) + (0.3 * green2) + rnorm(n, 0, 40) # Third green space distance exposure
    }
    
    ## If exposure is continuous and collinearity high
    if (Exposure == "Cont" & Collinear == "High") {
      high_sep <- rbinom(n, 1, 0.5) # SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
      green1 <- 325 + (-50 * high_sep) + rnorm(n, 0, 40) # First green space distance exposure in pregnancy
      green2 <- 50 + (-50 * high_sep) + (0.9 * green1) + rnorm(n, 0, 20) # Second green space distance exposure
      green3 <- 50 + (-50 * high_sep) + (0.9 * green2) + rnorm(n, 0, 20) # Third green space distance exposure
    }

    
    ## Encode the life course hypotheses
    
    # Critical periods (same for binary and continuous exposures)
    crit1 <- green1 # Critical period at first time point only
    int1 <- crit1 * high_sep # Interaction between SEP and first time point
    
    crit2 <- green2 # Critical period at second time point only
    int2 <- crit2 * high_sep # Interaction between SEP and second time point
    
    crit3 <- green3 # Critical period at third time point only
    int3 <- crit3 * high_sep # Interaction between SEP and third time point
    
    # Accumulation (different for binary and continuous exposures)
    if (Exposure == "Binary") {
      accumulation <- green1 + green2 + green3 # Linear accumulation of all exposures
      int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
    }
    
    if (Exposure == "Cont") {
      accumulation <- (green1 + green2 + green3) / 3 # Linear accumulation of all exposures
      int_accum <- high_sep * accumulation # Interaction between SEP and cumulative exposure
    }
    
    ## Create the outcomes (differs depending on whether exposures are binary or continuous)
    if (Exposure == "Binary") {
      bmi <- 25 + (-4 * high_sep) + (-2 * green3) + (2 * high_sep * green3) + rnorm(n, 0, 3)
    }
    
    if (Exposure == "Cont") {
      bmi <- 28 + (-4 * high_sep) + (-0.02 * green3) + (0.02 * high_sep * green3) + rnorm(n, 0, 3)
    }
    
    # If outcome is binary
    if (Outcome == "Binary") {
      overweight <- ifelse(bmi > 25, 1, 0) # Code BMI to binary 'overweight' variable if > 25
    }
    
    
    ## Center the exposures, if specified
    if (Centered == "Yes") {
      crit1 <- crit1 - mean(crit1) # Center critical period at time 1
      int1 <- crit1 * high_sep # Interaction between SEP and first time point
      
      crit2 <- crit2 - mean(crit2) # Center critical period at time 2
      int2 <- crit2 * high_sep # Interaction between SEP and second time point
      
      crit3 <- crit3 - mean(crit3) # Center critical period at time 3
      int3 <- crit3 * high_sep # Interaction between SEP and third time point
      
      accumulation <- accumulation - mean(accumulation) # Center the accumulation variable
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
    
    ## Loop over each step of the lasso, running a standard LM/GLM model (depending on outcome) whenever variables change. Store the AIC and BIC values and select the model with the lowest values
    
    # Make a dataframe of all the models identified by the lasso
    old_covars <- ""
    old_deviance <- 0
    old_varNum <- 0
    old_lambda <- NA
    df <- data.frame(matrix(ncol = 8, nrow = 0))
    #df
    
    for (j in 1:length(mod$df)) {
      #print(j)
      new_covars <- attributes(which(mod$beta[, j] != 0))$names
      new_deviance <- mod$dev.ratio[j]
      new_varNum <- mod$df[j]
      new_lambda <- mod$lambda[j]
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
    
    # Make a var to show number of steps where variables added, and rename the high_sep variable to blank (as is included by default)
    df$steps <- 1:nrow(df)
    df$Variables[df$steps == 1] <- ""
    #df
    
    
    ## Now run all combinations of the model variables in a standard LM/GLM, and store AIC and BIC values
    for (j in 1:nrow(df)) {
      
      vars_temp <- strsplit(df$model_vars[j], " ")[[1]] # Split the variables at each stage of the lasso
      x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp] # Matrix of the covariates at each step of lasso
      
      # If continuous outcome
      if (Outcome == "Cont") {
        mod_ic <- lm(bmi ~ x_hypos_new) # Run the model
        #print(summary(mod_ic))
      }
      
      # If binary outcome
      if (Outcome == "Binary") {
        mod_ic <- glm(overweight ~ x_hypos_new, family = "binomial") # Run the model
        #print(summary(mod_ic))
      }
      
      # Store the AIC values
      df$aic[j] <- AIC(mod_ic)
      df$bic[j] <- BIC(mod_ic)
      
    }
    
    # Select the models with the lowest AIC and BIC values and store these
    #df
    
    df$model_vars[which.min(df$aic)]
    df$model_vars[which.min(df$bic)]
    
    
    # Print the covariates in the best-fitting models if want to print this
    if (Output == TRUE) {
      print(paste0("Covariates in best-fitting AIC model:"))
      print(df$model_vars[which.min(df$aic)])
      
      print(paste0("Covariates in best-fitting BIC model:"))
      print(df$model_vars[which.min(df$bic)])
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
    vars_split_aic <- strsplit(df$model_vars[which.min(df$aic)], " ")[[1]] # Split the variables
    AIC_res_temp[i] <- ifelse(setequal(vars_split_aic, target_covars) == TRUE, 1, 0)
    
    vars_split_bic <- strsplit(df$model_vars[which.min(df$bic)], " ")[[1]] # Split the variables
    BIC_res_temp[i] <- ifelse(setequal(vars_split_bic, target_covars) == TRUE, 1, 0)
    
    # Store if 'crit3' in final model
    AIC_res_temp_crit3[i] <- ifelse("crit3" %in% vars_split_aic == TRUE, 1, 0)
    BIC_res_temp_crit3[i] <- ifelse("crit3" %in% vars_split_bic == TRUE, 1, 0)
    
    # Store if 'int3' in final model
    AIC_res_temp_int3[i] <- ifelse("int3" %in% vars_split_aic == TRUE, 1, 0)
    BIC_res_temp_int3[i] <- ifelse("int3" %in% vars_split_bic == TRUE, 1, 0)
    
    # Store if 'crit3' and 'int3' both in final model
    AIC_res_temp_crit3int3[i] <- ifelse("crit3" %in% vars_split_aic == TRUE & "int3" %in% vars_split_aic == TRUE, 1, 0)
    BIC_res_temp_crit3int3[i] <- ifelse("crit3" %in% vars_split_bic == TRUE & "int3" %in% vars_split_bic == TRUE, 1, 0)
    
    # Store if 'crit3' and 'int3' both in final model, plus other variables
    AIC_res_temp_crit3int3extra[i] <- ifelse(AIC_res_temp_crit3int3[i] == 1 & AIC_res_temp[i] == 0, 1, 0)
    BIC_res_temp_crit3int3extra[i] <- ifelse(BIC_res_temp_crit3int3[i] == 1 & BIC_res_temp[i] == 0, 1, 0)
    
    
    ### Next, want to summarise the 1SE and minimum MSE cross-validated lasso model, and see whether that corresponds to the correct model or not
    
    # If outcome is continuous
    if (Outcome == "Cont") {
      mod.cv <- cv.glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    # If outcome is binary
    if (Outcome == "Binary") {
      mod.cv <- cv.glmnet(x_hypos, overweight, alpha = 1, family = "binomial", 
                          penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
    }
    
    
    ## Start with the 1 SE model
    # Extract the non-zero covariates into a vector
    cv_1SE_covars_temp <- as.matrix(coef(mod.cv, s = mod.cv$lambda.1se))
    cv_1SE_covars_temp <- as.matrix(cv_1SE_covars_temp[!rownames(cv_1SE_covars_temp) %in% "(Intercept)", ]) # Drop the intercept
    cv_1SE_covars_temp <- cv_1SE_covars_temp[cv_1SE_covars_temp != 0, ] # Drop all zero coefficients
    
    cv_1SE_covars <- attributes(cv_1SE_covars_temp)$names # Store all the non-zero hypotheses
    
    # Print the covariates in the best-fitting cross-validated model
    if (Output == TRUE) {
      print(paste0("CV 1 SE covariates:"))
      print(cv_1SE_covars)
      print("")
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
    CV_1SE_res_temp[i] <- ifelse(setequal(cv_1SE_covars, target_covars) == TRUE, 1, 0)
    
    # Store if 'crit3', 'int3', both, and both with additional variables in final model
    CV_1SE_res_temp_crit3[i] <- ifelse("crit3" %in% cv_1SE_covars == TRUE, 1, 0)
    CV_1SE_res_temp_int3[i] <- ifelse("int3" %in% cv_1SE_covars == TRUE, 1, 0)
    CV_1SE_res_temp_crit3int3[i] <- ifelse("crit3" %in% cv_1SE_covars == TRUE & "int3" %in% cv_1SE_covars == TRUE, 1, 0)
    CV_1SE_res_temp_crit3int3extra[i] <- ifelse(CV_1SE_res_temp_crit3int3[i] == 1 & CV_1SE_res_temp[i] == 0, 1, 0)
    
    
    ## Next to the minimum MSE model
    # Extract the non-zero covariates into a vector
    cv_minMSE_covars_temp <- as.matrix(coef(mod.cv, s = mod.cv$lambda.min))
    cv_minMSE_covars_temp <- as.matrix(cv_minMSE_covars_temp[!rownames(cv_minMSE_covars_temp) %in% "(Intercept)", ]) # Drop the intercept
    cv_minMSE_covars_temp <- cv_minMSE_covars_temp[cv_minMSE_covars_temp != 0, ] # Drop all zero coefficients
    
    cv_minMSE_covars <- attributes(cv_minMSE_covars_temp)$names # Store all the non-zero hypotheses
    
    # Print the covariates in the best-fitting cross-validated model
    if (Output == TRUE) {
      print(paste0("CV min MSE covariates:"))
      print(cv_minMSE_covars)
      print("")
    }
    
    ## See whether this matches the 'true' model, and code as 0 if not and 1 if so
    CV_minMSE_res_temp[i] <- ifelse(setequal(cv_minMSE_covars, target_covars) == TRUE, 1, 0)
    
    # Store if 'crit3', 'int3', both, and both with additional variables in final model
    CV_minMSE_res_temp_crit3[i] <- ifelse("crit3" %in% cv_minMSE_covars == TRUE, 1, 0)
    CV_minMSE_res_temp_int3[i] <- ifelse("int3" %in% cv_minMSE_covars == TRUE, 1, 0)
    CV_minMSE_res_temp_crit3int3[i] <- ifelse("crit3" %in% cv_minMSE_covars == TRUE & "int3" %in% cv_minMSE_covars == TRUE, 1, 0)
    CV_minMSE_res_temp_crit3int3extra[i] <- ifelse(CV_minMSE_res_temp_crit3int3[i] == 1 & CV_minMSE_res_temp[i] == 0, 1, 0)
    
  }
  
  # Store the summaries of these results to transfer to the main results table
  res <- data.frame(AIC_propcorrect = round(sum(AIC_res_temp) / n_sims * 100, 2),
                    AIC_crit3correct = round(sum(AIC_res_temp_crit3) / n_sims * 100, 2),
                    AIC_int3correct = round(sum(AIC_res_temp_int3) / n_sims * 100, 2),
                    AIC_crit3int3correct = round(sum(AIC_res_temp_crit3int3) / n_sims * 100, 2),
                    AIC_crit3int3extra = round(sum(AIC_res_temp_crit3int3extra) / n_sims * 100, 2),
                    BIC_propcorrect = round(sum(BIC_res_temp) / n_sims * 100, 2),
                    BIC_crit3correct = round(sum(BIC_res_temp_crit3) / n_sims * 100, 2),
                    BIC_int3correct = round(sum(BIC_res_temp_int3) / n_sims * 100, 2),
                    BIC_crit3int3correct = round(sum(BIC_res_temp_crit3int3) / n_sims * 100, 2),
                    BIC_crit3int3extra = round(sum(BIC_res_temp_crit3int3extra) / n_sims * 100, 2),
                    CV_1SE_propcorrect = round(sum(CV_1SE_res_temp) / n_sims * 100, 2),
                    CV_1SE_crit3correct = round(sum(CV_1SE_res_temp_crit3) / n_sims * 100, 2),
                    CV_1SE_int3correct = round(sum(CV_1SE_res_temp_int3) / n_sims * 100, 2),
                    CV_1SE_crit3int3correct = round(sum(CV_1SE_res_temp_crit3int3) / n_sims * 100, 2),
                    CV_1SE_crit3int3extra = round(sum(CV_1SE_res_temp_crit3int3extra) / n_sims * 100, 2),
                    CV_minMSE_propcorrect = round(sum(CV_minMSE_res_temp) / n_sims * 100, 2),
                    CV_minMSE_crit3correct = round(sum(CV_minMSE_res_temp_crit3) / n_sims * 100, 2),
                    CV_minMSE_int3correct = round(sum(CV_minMSE_res_temp_int3) / n_sims * 100, 2),
                    CV_minMSE_crit3int3correct = round(sum(CV_minMSE_res_temp_crit3int3) / n_sims * 100, 2),
                    CV_minMSE_crit3int3extra = round(sum(CV_minMSE_res_temp_crit3int3extra) / n_sims * 100, 2))
  return(res)
  
}


## Next, the number of simulations per combination of parameters (1,000), the target hypotheses we simulated are the true model, and set up a data frame to store the results in
n_sims <- 10
set.seed(6789)

target_covars <- c("high_sep", "crit3", "int3")

results_reduced <- data.frame(sampleSize = rep(c(1000, 10000), 16),
                      Exposure = rep(c(rep("Binary", 8), rep("Cont", 8)), 2),
                      Centered = rep(c("No", "No", "Yes", "Yes"), 8),
                      Collinear = rep(c(rep("Low", 4), rep("High", 4)), 4),
                      Outcome = c(rep("Cont", 16), rep("Binary", 16)),
                      AIC_propcorrect = rep(NA, 32),
                      AIC_crit3correct = rep(NA, 32),
                      AIC_int3correct = rep(NA, 32),
                      AIC_crit3int3correct = rep(NA, 32),
                      AIC_crit3int3extra = rep(NA, 32),
                      BIC_propcorrect = rep(NA, 32),
                      BIC_crit3correct = rep(NA, 32),
                      BIC_int3correct = rep(NA, 32),
                      BIC_crit3int3correct = rep(NA, 32),
                      BIC_crit3int3extra = rep(NA, 32),
                      CV_1SE_propcorrect = rep(NA, 32),
                      CV_1SE_crit3correct = rep(NA, 32),
                      CV_1SE_int3correct = rep(NA, 32),
                      CV_1SE_crit3int3correct = rep(NA, 32),
                      CV_1SE_crit3int3extra = rep(NA, 32),
                      CV_minMSE_propcorrect = rep(NA, 32),
                      CV_minMSE_crit3correct = rep(NA, 32),
                      CV_minMSE_int3correct = rep(NA, 32),
                      CV_minMSE_crit3int3correct = rep(NA, 32),
                      CV_minMSE_crit3int3extra = rep(NA, 32))

results_reduced


# Also calculate time taken to run script
start_time <- Sys.time()


### First simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", 
                         Collinear = "Low",  Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 1
results_reduced[k, 6:25] <- res


### Second simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 2
results_reduced[k, 6:25] <- res


### Third simulation: Sample size = 1000; binary exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 3
results_reduced[k, 6:25] <- res


### Fourth simulation: Sample size = 10000; binary exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "Low",  Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 4
results_reduced[k, 6:25] <- res


### Fifth simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 5
results_reduced[k, 6:25] <- res


### Sixth simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 6
results_reduced[k, 6:25] <- res


### Seventh simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 7
results_reduced[k, 6:25] <- res


### Eighth simulation: Sample size = 10000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 8
results_reduced[k, 6:25] <- res


### Ninth simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 9
results_reduced[k, 6:25] <- res


### Tenth simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 10
results_reduced[k, 6:25] <- res


### Eleventh simulation: Sample size = 1000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 11
results_reduced[k, 6:25] <- res


### Twelfth simulation: Sample size = 10000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 12
results_reduced[k, 6:25] <- res


### 13th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 13
results_reduced[k, 6:25] <- res


### 14th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 14
results_reduced[k, 6:25] <- res


### 15th simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 15
results_reduced[k, 6:25] <- res


### 16th simulation: Sample size = 10000; continuous exposure; centered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "High", Outcome = "Cont", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 16
results_reduced[k, 6:25] <- res


### 17th simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 17
results_reduced[k, 6:25] <- res


### 18th simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 18
results_reduced[k, 6:25] <- res


### 19th simulation: Sample size = 1000; binary exposure; centered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 19
results_reduced[k, 6:25] <- res


### 20th simulation: Sample size = 10000; binary exposure; centered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 20
results_reduced[k, 6:25] <- res


### 21st simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 21
results_reduced[k, 6:25] <- res


### 22nd simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 22
results_reduced[k, 6:25] <- res


### 23rd simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 23
results_reduced[k, 6:25] <- res


### 24th simulation: Sample size = 10000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 24
results_reduced[k, 6:25] <- res


### 25th simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 25
results_reduced[k, 6:25] <- res


### 26th simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 26
results_reduced[k, 6:25] <- res


### 27th simulation: Sample size = 1000; continuous exposure; centered; low collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 27
results_reduced[k, 6:25] <- res


### 28th simulation: Sample size = 10000; continuous exposure; centered; low collinearity; continuous outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes",
                         Collinear = "Low", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 28
results_reduced[k, 6:25] <- res


### 29th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 29
results_reduced[k, 6:25] <- res


### 30th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 30
results_reduced[k, 6:25] <- res


### 31st simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 31
results_reduced[k, 6:25] <- res


### 32nd simulation: Sample size = 10000; continuous exposure; centered; High collinearity; binary outcome
res <- lasso_sim_reduced(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", 
                         Collinear = "High", Outcome = "Binary", Output = TRUE)

# Store the summaries of these results in the main results table
k <- 32
results_reduced[k, 6:25] <- res


# Time taken to run script (took about 40 mins for 100 simulations, so should be about 6 hours for the full 1,000 simulations)
end_time <- Sys.time()

end_time - start_time


## Save the results table
results_reduced

write.csv(results_reduced, file = "simulationResults_reduced_test.csv", row.names = FALSE)


### Get some results from this
#results_reduced <- read.csv("simulationResults_reduced_test.csv")
#head(results_reduced)


