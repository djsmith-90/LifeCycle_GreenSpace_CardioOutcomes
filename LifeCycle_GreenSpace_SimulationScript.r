### Example simulation script for structured lifecourse models with confounders and interactions

### Created 17/2/2022 by Dan Major-Smith
### R version 4.0.4


### Previous work has detailed a structured approach to life course modelling for both binary exposures (Smith et al., 2015: https://journals.lww.com/epidem/Fulltext/2015/09000/Model_Selection_of_the_Effect_of_Binary_Exposures.15.aspx) and continuous exposures with confounding (Smith et al., 2016: https://academic.oup.com/ije/article/45/4/1271/2951966?login=true) using LARS/lasso methods. Building on these approaches, we aim to demonstrate here how interactions can be incorporated into these models. 

## However, we will not be using the LARS/covariance test method here, as: 1) there are issues with the covariance test and it is no longer recommended; and 2) Using the LARS method, it is possible to 'force' continuous confounders to be included in the model, however, this does not appear possible for binary/categorical variables, especially if we want to include them as interaction terms. So here we will focus on using ordinary lasso models via glmnet, rather than the LARS method. Rather than being a stepwise procedure like LARS, this approach gradually increases the lambda value and lets more variables into the model; this can make it difficult to assess when adding a new covariate does or does not improve model fit. Interpretation is therefore more subjective, based on inspection of the improvement in the deviance ratio (a measure of goodness-of-fit similar to R2).

## To interpret these results, we will focus on three approaches:
##  - 1) Using a 'subjective' approach looking at the order in which hypotheses were entered into the model, combined with a plot of deviance/variance explained when each predictor was added, and making judgement based on these sources of information
##  - 2) Taking each variable(s) entered in turn in the lasso model, using a likelihood ratio test to compare standard linear/logistic models with/without the next predictor in. This will provide a formal test as to whether the hypothesis predicts the outcome.
##  - 3) Using cross-validated lasso and selecting the model within 1 SE of the best-fitting model (this method is similar to ordinary lasso, but uses k-fold cross-validation to identify the best-fitting model while minimising overfitting). Using the model within 1 SE of the best-fitting model, rather than the best-fitting model itself, should also avoid potential overfitting, leaving only variables more strongly associated with the outcome in the final model.

## If all these methods give a similar answer, then we can have more confidence in the conclusions drawn.



###################################################################################################################
#### Clear workspace and install/load packages

rm(list = ls())

#install.packages("glmnet")
#install.packages("dagitty")

library(glmnet)
library(dagitty)


## Working directory
setwd("C:\\Users\\ds16565\\OneDrive - University of Bristol\\MyFiles-Migrated\\Documents\\Projects\\Lifecycle\\SimResults")



####################################################################################################################
#### Simulating life course data with interaction terms

### Here, will assume exposures and covariates/interaction terms are binary, with a continuous outcome.
## Exposure: Access to green space (within 300 metres) in pregnancy, age 4 and age 7 (0 = no access; 1 = access).
## Confounder/interaction term: Socioeconomic position (0 = low; 1 = high)
## Outcome: BMI age 7 (continuous)

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

# Save this DAG
pdf("GreenSpaceDAG_simulation.pdf", height = 6, width = 8)
plot(dag)
dev.off()


### Generate the data (although note also that the simulated data here is purely to illustrate the logic and application of these structure life-course methods, and should not in any way be taken as a reflection of the real-world patterns or effect sizes of these variables)

# Sample size of ~10,000 (approximate ALSPAC participation early in study)
set.seed(4321)
n <- 10000

# SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
high_sep <- rbinom(n, 1, 0.5)
table(high_sep, useNA = "ifany")

# First green space exposure in pregnancy - Caused by SEP (higher SEP = more exposure to green space; three times the odds) - Approx. 50% near green space
green1_logit <- log(0.6) + (log(3) * high_sep)
green1_p <- exp(green1_logit) / (1 + exp(green1_logit))

green1 <- rbinom(n, 1, green1_p)
table(green1, useNA = "ifany")

table(high_sep, green1)

rm(green1_logit)
rm(green1_p)

# Second green space exposure at age 4 - Caused by Green1 and SEP
green2_logit <- log(0.3) + (log(3) * high_sep) + (log(5) * green1)
green2_p <- exp(green2_logit) / (1 + exp(green2_logit))

green2 <- rbinom(n, 1, green2_p)
table(green2, useNA = "ifany")

table(high_sep, green2)
table(green1, green2)

rm(green2_logit)
rm(green2_p)

# Third green space exposure at age 7 - Caused by Green2 and SEP
green3_logit <- log(0.3) + (log(3) * high_sep) + (log(5) * green2)
green3_p <- exp(green3_logit) / (1 + exp(green3_logit))

green3 <- rbinom(n, 1, green3_p)
table(green3, useNA = "ifany")

table(high_sep, green3)
table(green1, green3)
table(green2, green3)

rm(green3_logit)
rm(green3_p)

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green1 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green1 here.
bmi <- 25 + (-4 * high_sep) + (-2 * green1) + (2 * high_sep * green1) + rnorm(n, 0, 3)
summary(bmi)

## descriptive stats split by possible combinations of SEP and green1 to check simulation worked
summary(bmi[high_sep == 1 & green1 == 1]) # High SEP and access to green space in preg - Lowest BMI
summary(bmi[high_sep == 1 & green1 == 0]) # High SEP and no access to green space in preg - Lowest BMI
summary(bmi[high_sep == 0 & green1 == 1]) # Low SEP and access to green space in preg - Middle BMI
summary(bmi[high_sep == 0 & green1 == 0]) # Low SEP and no access to green space in preg - Highest BMI

# Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green1))


### Encode the life course hypotheses

# Critical period at first time point only
crit1 <- green1

# Interaction between SEP and first time point
int1 <- high_sep * green1

# Critical period at second time point only
crit2 <- green2

# Interaction between SEP and second time point
int2 <- high_sep * green2

# Critical period at third time point only
crit3 <- green3

# Interaction between SEP and third time point
int3 <- high_sep * green3

# Linear accumulation of all exposures
accumulation <- green1 + green2 + green3

# Interaction between SEP and cumulative exposure
int_accum <- high_sep * accumulation

# Increase in access to green space from time 1 to time 2
green_inc12 <- (1 - green1) * green2

# Increase in access to green space from time 1 to time 2, with an interaction with SEP
green_inc12_int <- (1 - green1) * green2 * high_sep

# Decrease in access to green space from time 1 to time 2
green_dec12 <- (1 - green2) * green1

# Decrease in access to green space from time 1 to time 2, with an interaction with SEP
green_dec12_int <- (1 - green2) * green1 * high_sep

# Increase in access to green space from time 2 to time 3
green_inc23 <- (1 - green2) * green3

# Increase in access to green space from time 2 to time 3, with an interaction with SEP
green_inc23_int <- (1 - green2) * green3 * high_sep

# Decrease in access to green space from time 2 to time 3
green_dec23 <- (1 - green3) * green2

# Decrease in access to green space from time 2 to time 3, with an interaction with SEP
green_dec23_int <- (1 - green3) * green2 * high_sep


## Combine all these into one matrix, including SEP as a confounder
x_hypos <- cbind(high_sep, crit1, int1, crit2, int2, crit3, int3, accumulation, int_accum, green_inc12, 
                 green_inc12_int, green_dec12, green_dec12_int, green_inc23, green_inc23_int, green_dec23,
                 green_dec23_int)
head(x_hypos)


## Check correlation matrix of all these hypotheses (>0.9 would be a cause for concern) - correlation between 'int_accum' and 'int2' is 0.91, but that should be fine here.
dim(x_hypos)
cor(x_hypos[,2:17])
cor(x_hypos[,2:17]) > 0.9


### Run the model using glmnet. 
# alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives 'high_sep' a weighting of '0', so is always included in the model (default is 1)
mod <- glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod

# Plot these results
plot(mod)

# Look at the variables included at each step (first model is the baseline SEP-only model)
coef(mod, s = max(mod$lambda[mod$df == 1])); min(mod$dev.ratio[mod$df == 1])
coef(mod, s = max(mod$lambda[mod$df == 2])); min(mod$dev.ratio[mod$df == 2])
coef(mod, s = max(mod$lambda[mod$df == 3])); min(mod$dev.ratio[mod$df == 3])
coef(mod, s = max(mod$lambda[mod$df == 4])); min(mod$dev.ratio[mod$df == 4])
coef(mod, s = max(mod$lambda[mod$df == 5])); min(mod$dev.ratio[mod$df == 5])

# The first variable included was green1, and then the second variable included was the interaction term between SEP and green1, which is correct, as this is how we coded the data.


#### Will explore a few methods to check whether we get the right result and do not include any other incorrect parameters in the 'best' model

### 1) Visual inspection

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 5, nrow = 0))
df

for (i in 1:length(mod$df)) {
  #print(i)
  new_covars <- attributes(which(mod$beta[, i] != 0))$names
  new_deviance <- mod$dev.ratio[i]
  new_varNum <- mod$df[i]
  new_lambda <- mod$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if covars added
  if (new_varNum > old_varNum) {
    change <- setdiff(new_covars, old_covars) # Find the new covariate(s) added
    change <- paste(change, collapse = " ") # Combine variable together, if > 1
    change <- paste0(change, " (+)") # Append a "(+)" sign
    dev_diff <- round((new_deviance - old_deviance) * 100, 3) # Diff in deviance between current and previous lambda
    new_dev <- round(new_deviance * 100, 3) # Current deviance value
    temp <- cbind(change, new_dev, dev_diff, new_varNum, new_lambda) # Combine values together
    df <- rbind(df, temp) # Merge with template data frame
  }
  
  # See if covars removed
  if (new_varNum < old_varNum) {
    change <- setdiff(old_covars, new_covars) # Find the covariate(s) removed
    change <- paste(change, collapse = " ") # Combine variable together, if > 1
    change <- paste0(change, " (-)") # Append a "(-)" sign
    dev_diff <- round((new_deviance - old_deviance) * 100, 3) # Diff in deviance between current and previous lambda
    new_dev <- round(new_deviance * 100, 3) # Current deviance value
    temp <- cbind(change, new_dev, dev_diff, new_varNum, new_lambda) # Combine values together
    df <- rbind(df, temp) # Merge with template data frame
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevRatio", "DevDiff", "VarNum", "Lambda")
df

# Make a var to show number of steps where variables added, and rename the covars
df$steps <- 1:nrow(df)
df$Variables[df$steps == 1] <- "covars"
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. This makes it clearer when different variables were added (not the neatest plot, though, as variables bunch up and overlap as lambda approaches 0...). Also, need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod$lambda, mod$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod$lambda)), ylim = c(0, max(mod$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod$log_lambda <- log(mod$lambda)
mod

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, as now much clearer that 'crit1' and 'int1' were entered first. However, one potentially misleading aspect is that lasso works by initialising the lambda value just above the threshold where no variables are entered (excluding covariates restrained to be in the model). This means that there will always be a 'first' variable entered early in the model, making it seem like this is the best fit to the data. However, because there always *has* to be one variable entered first, it doesn't mean that this is actually predictive of the outcome. So to interpret this we need to look at the deviance ratio. Here, once 'crit1' is added the deviance ratio increases by about 1.5%, after which 'int1' is added, and the deviance ratio increases by around another 2% until the next variable is entered in the model (green_inc12_int), after which the deviance ratio does not increase at all. This suggests that 'crit1' and 'int1' in combination explain most of the variation in the outcome BMI attributable to the life-course hypotheses, and that these variables are associated with the outcome (again, just as we simulated).
plot(mod$log_lambda, mod$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod$log_lambda)), ylim = c(0.2, max(mod$dev.ratio)))
text(df$log_lambda, 0.2, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_simulation.pdf", height = 7, width = 11)
plot(mod$log_lambda, mod$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod$log_lambda)), ylim = c(0.2, max(mod$dev.ratio)))
text(df$log_lambda, 0.2, labels = df$Variables, srt = 90, adj = 0)
dev.off()



### 2) Likelihood ratio tests at inclusion of each new parameter (if only one parameter is added at a time we don't really we don't need to do an LR test, as for single parameters the p-value from the model will be identical to those of the LR test, but we're doing this LR test because sometimes multiple terms get added to the lasso at one time-point)

## First, compare the baseline (confounder only) model to the model including the first term added (crit1)
coef(mod, s = max(mod$lambda[mod$df == 1])); min(mod$dev.ratio[mod$df == 1])
coef(mod, s = max(mod$lambda[mod$df == 2])); min(mod$dev.ratio[mod$df == 2])

base <- lm(bmi ~ x_hypos[, "high_sep"])
summary(base)
param1 <- lm(bmi ~ x_hypos[, "high_sep"] + x_hypos[, "crit1"])
summary(param1)

# Find strong support for more complex model which includes crit1 main effect
anova(base, param1)


## Next, compare the model including the next hypotheses added (int1)
coef(mod, s = max(mod$lambda[mod$df == 3])); min(mod$dev.ratio[mod$df == 3])

param2 <- lm(bmi ~ x_hypos[, "high_sep"] + x_hypos[, "crit1"] + x_hypos[, "int1"])
summary(param2)

# Adding the SEP by crit1 interaction term again improves model fit
anova(param1, param2)


## Does adding the next parameter (green_inc12_int) improve model fit?
coef(mod, s = max(mod$lambda[mod$df == 4])); min(mod$dev.ratio[mod$df == 4])

param3 <- lm(bmi ~ x_hypos[, "high_sep"] + x_hypos[, "crit1"] + x_hypos[, "int1"] + x_hypos[, "green_inc12_int"])
summary(param3)

# No improvement in model fit (as expected)
anova(param2, param3)



### 3) Cross-validated lasso to find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model. Will use the default number of k-folds (10).
set.seed(3930)
mod.cv <- cv.glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
mod.cv

plot(mod.cv) # Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters runs along the top of the plot)

# The 1SE model contains 3 parameters (SEP confounder, crit1 and int1 term), so identifies the correct model
coef(mod.cv, s = mod.cv$lambda.1se)

# The optimal' model contains 5 parameters, so does not correct identify the correct model (but as this optimal lasso is mainly for prediction, this increase in complexity is to be expected as the model hooks on to random noise in the data)
coef(mod.cv, s = mod.cv$lambda.min)

# Save this plot
pdf(file = "CVLasso_simulation.pdf", height = 7, width = 11)
plot(mod.cv)
dev.off()


#### Summary: All of the methods seem to match up well, give consistent answers, and identify the true model correctly.



###################################################################################################################
###################################################################################################################
#### Example life course model with confounders and interactions, but with a binary outcome

### Create the binary 'overweight' outcome - Will use the same model as for BMI - Overweight caused by SEP (higher SEP = lower probability of being overweight), plus interaction with green1 (lower SEP and access to green space = lower chances of being overweight compared to lower SEP and no access to green space). Assuming no main effect of green1 here.
set.seed(3930)

overweight <- NA
overweight[high_sep == 1 & green1 == 1] <- rbinom(sum(high_sep == 1 & green1 == 1), 1, 0.1)
overweight[high_sep == 1 & green1 == 0] <- rbinom(sum(high_sep == 1 & green1 == 0), 1, 0.1)
overweight[high_sep == 0 & green1 == 1] <- rbinom(sum(high_sep == 0 & green1 == 1), 1, 0.3)
overweight[high_sep == 0 & green1 == 0] <- rbinom(sum(high_sep == 0 & green1 == 0), 1, 0.5)
table(overweight, useNA = "ifany")

table(overweight, green1)
table(overweight, high_sep)

## descriptive stats split by possible combinations of SEP and green1 to check simulation worked
mean(overweight[high_sep == 1 & green1 == 1]) # High SEP and access to green space in preg - Low prob of overweight
mean(overweight[high_sep == 1 & green1 == 0]) # High SEP and no access to green space in preg - Low prob of overweight
mean(overweight[high_sep == 0 & green1 == 1]) # Low SEP and access to green space in preg - Middle prob of overweight
mean(overweight[high_sep == 0 & green1 == 0]) # Low SEP and no access to green space in preg - High prob of overweight

# Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(glm(overweight ~ high_sep + green1 + green2 + green3, family = "binomial"))
summary(glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green1, family = "binomial"))
exp(cbind(coef(summary(glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green1, family = "binomial")))[, 1], confint(glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green1, family = "binomial"))))
anova(glm(overweight ~ high_sep + green1 + green2 + green3, family = "binomial"), 
      glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green1, family = "binomial"), test = "Chisq")


### Everything else is the same as above, just using a logistic rather than linear lasso model

### Run the model using glmnet
# alpha = 1 specifies L1 regularisation (lasso model), and the penalty factor option gives 'high_sep' a weighting of '0', so is always included in the model (default is 1)
mod.binary <- glmnet(x_hypos, overweight, alpha = 1, family = "binomial", 
              penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))

# List the results of this model - Number of vars included, % deviance explained, and lambda value
mod.binary

# Plot these results
plot(mod.binary)

# Look at the variables included at each step (first model is the baseline SEP-only model)
coef(mod.binary, s = max(mod.binary$lambda[mod.binary$df == 1])); min(mod.binary$dev.ratio[mod.binary$df == 1])
coef(mod.binary, s = max(mod.binary$lambda[mod.binary$df == 2])); min(mod.binary$dev.ratio[mod.binary$df == 2])
coef(mod.binary, s = max(mod.binary$lambda[mod.binary$df == 3])); min(mod.binary$dev.ratio[mod.binary$df == 3])
coef(mod.binary, s = max(mod.binary$lambda[mod.binary$df == 4])); min(mod.binary$dev.ratio[mod.binary$df == 4])
coef(mod.binary, s = max(mod.binary$lambda[mod.binary$df == 5])); min(mod.binary$dev.ratio[mod.binary$df == 5])

# The first variable included was green1, and then the second variable included was the interaction term between SEP and green1, which is correct, as this is how we coded the data.


#### Will explore a few methods to check whether we get the right result and do not include any other incorrect parameters in the 'best' model

### 1) Visual inspection

# First, make a loop to pick out the changes in variables and the increment in deviance ratio
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 5, nrow = 0))
df

for (i in 1:length(mod.binary$df)) {
  #print(i)
  new_covars <- attributes(which(mod.binary$beta[, i] != 0))$names
  new_deviance <- mod.binary$dev.ratio[i]
  new_varNum <- mod.binary$df[i]
  new_lambda <- mod.binary$lambda[i]
  #print(new_covars); print(new_deviance); print(new_varNum)
  
  # See if covars added
  if (new_varNum > old_varNum) {
    change <- setdiff(new_covars, old_covars) # Find the new covariate(s) added
    change <- paste(change, collapse = " ") # Combine variable together, if > 1
    change <- paste0(change, " (+)") # Append a "(+)" sign
    dev_diff <- round((new_deviance - old_deviance) * 100, 3) # Diff in deviance between current and previous lambda
    new_dev <- round(new_deviance * 100, 3) # Current deviance value
    temp <- cbind(change, new_dev, dev_diff, new_varNum, new_lambda) # Combine values together
    df <- rbind(df, temp) # Merge with template data frame
  }
  
  # See if covars removed
  if (new_varNum < old_varNum) {
    change <- setdiff(old_covars, new_covars) # Find the covariate(s) removed
    change <- paste(change, collapse = " ") # Combine variable together, if > 1
    change <- paste0(change, " (-)") # Append a "(-)" sign
    dev_diff <- round((new_deviance - old_deviance) * 100, 3) # Diff in deviance between current and previous lambda
    new_dev <- round(new_deviance * 100, 3) # Current deviance value
    temp <- cbind(change, new_dev, dev_diff, new_varNum, new_lambda) # Combine values together
    df <- rbind(df, temp) # Merge with template data frame
  }
  
  # Rename the old covars, deviance and variable number
  old_covars <- new_covars
  old_deviance <- new_deviance
  old_varNum <- new_varNum
}

colnames(df) <- c("Variables", "DevRatio", "DevDiff", "VarNum", "Lambda")
df

# Make a var to show number of steps where variables added, and rename the covars
df$steps <- 1:nrow(df)
df$Variables[df$steps == 1] <- "covars"
df


# Make a plot of deviance ratio by lambda value, to show the time variables were entered and the improvement in model fit. This makes it clearer when different variables were added (not the neatest plot, though, as variables bunch up and overlap as lambda approaches 0...). Also, need to be careful about scale of y-axis, as improvement in model fit may not that large.
plot(mod.binary$lambda, mod.binary$dev.ratio, type = "l",
     xlab = "Lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod.binary$lambda)), ylim = c(0, max(mod.binary$dev.ratio)))
text(df$Lambda, 0, labels = df$Variables, srt = 90, adj = 0)


## Put lambda on the log scale, else differences between early models inflated, as lambda decreases on log scale. This does make the plot slightly more readable, and groups early-included variables closer together
mod.binary$log_lambda <- log(mod.binary$lambda)
mod

df$log_lambda <- log(as.numeric(df$Lambda))
df

# This looks better, as now clearer that 'crit1' and 'int1' were entered first. However, need to be aware of the deviance ratio scale (as mentioned above). Here, once 'crit1' is added the deviance ratio increases a little (~1%), after which 'int1' is added, and the deviance ratio increases by around another ~0.3% until the next variable is entered in the model (green_dec12), after which the deviance ratio does not increase by very much. This suggests that 'crit1' and 'int1' in combination explain most of the variation in the outcome BMI attributable to the life-course hypotheses, and that these variables are associated with the outcome (again, just as we simulated). Although unlike with the continuous outcome example above, the deviance ratios/variance explained by these 'crit1' and 'int1' variables is quite a bit weaker here with this binary outcome.
plot(mod.binary$log_lambda, mod.binary$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod.binary$log_lambda)), ylim = c(0.1, max(mod.binary$dev.ratio)))
text(df$log_lambda, 0.1, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_binaryOutcome_simulation.pdf", height = 7, width = 11)
plot(mod.binary$log_lambda, mod.binary$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod.binary$log_lambda)), ylim = c(0.1, max(mod.binary$dev.ratio)))
text(df$log_lambda, 0.1, labels = df$Variables, srt = 90, adj = 0)
dev.off()



### 2) Likelihood ratio tests at inclusion of each new parameter (if only one parameter is added at a time we don't really we don't need to do an LR test, as for single parameters the p-value from the model will be identical to those of the LR test, but we're doing this LR test because sometimes multiple terms get added to the lasso at one time-point)

## First, compare the baseline (confounder only) model to the model including the first term added (crit1)
coef(mod.binary, s = max(mod.binary$lambda[mod.binary$df == 1])); min(mod.binary$dev.ratio[mod.binary$df == 1])
coef(mod.binary, s = max(mod.binary$lambda[mod.binary$df == 2])); min(mod.binary$dev.ratio[mod.binary$df == 2])

base.binary <- glm(overweight ~ x_hypos[, "high_sep"], family = "binomial")
summary(base.binary)
param1.binary <- glm(overweight ~ x_hypos[, "high_sep"] + x_hypos[, "crit1"], family = "binomial")
summary(param1.binary)

# Find strong support for more complex model which includes crit1 main effect
anova(base.binary, param1.binary, test = "Chisq")


## Next, compare the model including the next hypotheses added (int1)
coef(mod.binary, s = max(mod.binary$lambda[mod.binary$df == 3])); min(mod.binary$dev.ratio[mod.binary$df == 3])

param2.binary <- glm(overweight ~ x_hypos[, "high_sep"] + x_hypos[, "crit1"] + x_hypos[, "int1"], family = "binomial")
summary(param2.binary)

# Adding the SEP by crit1 interaction term again improves model fit
anova(param1.binary, param2.binary, test = "Chisq")


## Does adding the next parameter (green_dec12) improve model fit?
coef(mod.binary, s = max(mod.binary$lambda[mod.binary$df == 4])); min(mod.binary$dev.ratio[mod.binary$df == 4])

param3.binary <- glm(overweight ~ x_hypos[, "high_sep"] + x_hypos[, "crit1"] + x_hypos[, "int1"] + 
                       x_hypos[, "green_dec12"], family = "binomial")
summary(param3.binary)

# No improvement in model fit (as expected)
anova(param2.binary, param3.binary, test = "Chisq")



### 3) Cross-validated lasso to find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model. Will use the default number of k-folds (10).
set.seed(3931)
mod.cv.binary <- cv.glmnet(x_hypos, overweight, alpha = 1, family = "binomial", 
                    penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
mod.cv.binary

plot(mod.cv.binary) # Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters runs along the top of the plot)

# Here, this 1SE model contains 3 parameters (SEP confounder, crit1 and int1 term), so identifies the correct model (however, note that these methods are quite sensitive, and if we change the set.seed number, either when simulating the data or before running the cross-validated model, we may get different results; for instance, if the seed number above the CV lasso is set to 3932 rather than 3931, the number of terms in the 1SE model reduces to 2 [high_sep and crit1 only]. The 1SE model may therefore be too conservative, especially when effect sizes are relatively small)
coef(mod.cv.binary, s = mod.cv.binary$lambda.1se)

# The optimal' model contains 5 parameters, so does not correct identify the correct model (but as this optimal lasso is mainly for prediction, this increase in complexity is to be expected as the model hooks on to random noise in the data)
coef(mod.cv.binary, s = mod.cv.binary$lambda.min)

# Save this plot
pdf(file = "CVLasso_binaryOutcome_simulation.pdf", height = 7, width = 11)
plot(mod.cv.binary)
dev.off()


#### Summary: All of the methods seem to match up well, give consistent answers, and identify the true model correctly. However, this may not always be the case, given that results do appear quite sensitive to random variation and different starting seeds/conditions (this is also exemplified in the accompanying Stata example script, where different methods - especially for cross-validated lasso and binary outcomes - do not always identify the true model); interpretation therefore should be made on a qualitative judgement taking information from all these methods into consideration, rather than just one.


