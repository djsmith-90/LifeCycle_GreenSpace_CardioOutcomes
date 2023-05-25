### Example simulation script for structured lifecourse models with confounders and interactions

### Created 17/2/2022 by Dan Major-Smith
### R version 4.0.4


### Previous work has detailed a structured approach to life course modelling for both binary exposures (Smith et al., 2015: https://journals.lww.com/epidem/Fulltext/2015/09000/Model_Selection_of_the_Effect_of_Binary_Exposures.15.aspx) and continuous exposures with confounding (Smith et al., 2016: https://academic.oup.com/ije/article/45/4/1271/2951966?login=true) using LARS/lasso methods. Building on these approaches, we aim to demonstrate here how interactions can be incorporated into these models. 

## However, we will not be using the LARS/covariance test method here, as: 1) there are issues with the covariance test and it is no longer recommended; and 2) Using the LARS method, it is possible to 'force' continuous confounders to be included in the model, however, this is more difficult for binary/categorical outcomes. So here we will focus on using ordinary lasso models via glmnet, rather than the LARS method. Rather than being a stepwise procedure like LARS, this approach gradually increases the lambda value and lets more variables into the model; this can make it difficult to assess when adding a new covariate does or does not improve model fit. Interpretation is therefore more subjective, based on inspection of the improvement in the deviance ratio (a measure of goodness-of-fit similar to R2).

## To interpret these results, we will focus on three approaches:
##  - 1) Using a 'subjective' approach looking at the order in which hypotheses were entered into the model, combined with a plot of deviance/variance explained when each predictor was added, and making judgement based on these sources of information
##  - 2) Using a 'relaxed lasso'-type approach, where use a standard LM/GLM on the model the lasso selects at each step, then comparing model fit of all these models to detect the best-fitting model. Will use both AIC and BIC as measures of model fit.
##  - 3) Using cross-validated lasso and selecting the model within 1 SE of the best-fitting model (this method is similar to ordinary lasso, but uses k-fold cross-validation to identify the best-fitting model while minimising overfitting). Using the model within 1 SE of the best-fitting model, rather than the best-fitting model itself, should also avoid potential overfitting, leaving only variables more strongly associated with the outcome in the final model. Will also explore model with lowest mean-squared error as well, as a comparison.

## If all these methods give a similar answer, then we can have more confidence in the conclusions drawn.



###################################################################################################################
#### Clear workspace and install/load packages

rm(list = ls())

#install.packages("glmnet")
#install.packages("dagitty")
#install.packages(lars)
#install.packages("selectiveInference")

library(glmnet)
library(dagitty)
library(lars)
library(selectiveInference)


## Working directory
setwd("C:\\Users\\ds16565\\OneDrive - University of Bristol\\MyFiles-Migrated\\Documents\\Projects\\Lifecycle\\SimResults")



####################################################################################################################
#### Simulating life course data with interaction terms

### Here, will assume exposures and covariates/interaction terms are binary, with a continuous outcome.
## Exposure: Access to green space (within 300 metres) in pregnancy, age 4 and age 7 (0 = no access; 1 = access).
## Confounder/interaction term: Socioeconomic position (0 = low; 1 = high)
## Outcome: BMI age 7 (continuous)

## Hypothetical DAG - The model here is that SEP is a confounder, but also interacts with most recent green space exposure to impact BMI at age 7
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

# Save this DAG
pdf("GreenSpaceDAG_simulation.pdf", height = 6, width = 8)
plot(dag)
dev.off()


### Generate the data (although note also that the simulated data here is purely to illustrate the logic and application of these structure life-course methods, and should not be taken as a reflection of the real-world patterns or effect sizes of these variables)

# Sample size of ~10,000 (approximate ALSPAC participation early in study)
set.seed(1234)
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

# Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green3 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green3 here if high SEP.
bmi <- 25 + (-4 * high_sep) + (-2 * green3) + (2 * high_sep * green3) + rnorm(n, 0, 3)
summary(bmi)

## descriptive stats split by possible combinations of SEP and green1 to check simulation worked
summary(bmi[high_sep == 1 & green3 == 1]) # High SEP and recent access to green space - Lowest BMI
summary(bmi[high_sep == 1 & green3 == 0]) # High SEP and no recent access to green space - Lowest BMI
summary(bmi[high_sep == 0 & green3 == 1]) # Low SEP and recent access to green space - Middle BMI
summary(bmi[high_sep == 0 & green3 == 0]) # Low SEP and no recent access to green space - Highest BMI

# Check the 'true' model, which is SEP as confounder, interaction with green 3, and main effect of green 3. green1 and green2 should also be null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(lm(bmi ~ high_sep + green1 + green2 + green3))
summary(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
anova(lm(bmi ~ high_sep + green1 + green2 + green3), lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))

AIC(lm(bmi ~ high_sep + green1 + green2 + green3)); AIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))
BIC(lm(bmi ~ high_sep + green1 + green2 + green3)); BIC(lm(bmi ~ high_sep + green1 + green2 + green3 + high_sep*green3))


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

# The first variable included was green3, and then the second variable included was the interaction term between SEP and green3, which is correct, as this is how we coded the data.


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
    change_drop <- setdiff(old_covars, new_covars) # Make sure no covariates dropped at same time
    if (!identical(change_drop, character(0))) { # If a covar also dropped, then combine with 'change'
      change_drop <- paste(change_drop, collapse = " ") # Combine variable together, if > 1
      change <- paste0(change, " ", change_drop, " (-)")
    }
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
    change_add <- setdiff(new_covars, old_covars) # Make sure no covariates added at same time
    if (!identical(change_add, character(0))) { # If a covar also dropped, then combine with 'change'
      change_add <- paste(change_add, collapse = " ") # Combine variable together, if > 1
      change <- paste0(change, " ", change_add, " (+)")
    }
    dev_diff <- round((new_deviance - old_deviance) * 100, 3) # Diff in deviance between current and previous lambda
    new_dev <- round(new_deviance * 100, 3) # Current deviance value
    temp <- cbind(change, new_dev, dev_diff, new_varNum, new_lambda) # Combine values together
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

# Make a var to show number of steps where variables added, and rename the high_sep variable to blank (as is included by default)
df$steps <- 1:nrow(df)
df$Variables[df$steps == 1] <- ""
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

# This looks better, as now much clearer that 'crit3' and 'int3' were entered first. However, one potentially misleading aspect is that lasso works by initialising the lambda value just above the threshold where no variables are entered (excluding covariates restrained to be in the model). This means that there will always be a 'first' variable entered early in the model, making it seem like this is the best fit to the data. However, because there always *has* to be one variable entered first, it doesn't mean that this is actually predictive of the outcome. So to interpret this we need to look at the deviance ratio. Here, once 'crit3' is added the deviance ratio increases by about 1.5%, after which 'int3' is added, and the deviance ratio increases by around another 1.5% until the next variable is entered in the model (green_dec23_int), after which the deviance ratio barely increases at all. This suggests that 'crit3' and 'int3' in combination explain most of the variation in the outcome BMI attributable to the life-course hypotheses, and that these variables are associated with the outcome (again, just as we simulated).
plot(mod$log_lambda, mod$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod$log_lambda)), ylim = c(0.2, max(mod$dev.ratio)))
text(df$log_lambda, 0.2, labels = df$Variables, srt = 90, adj = 0)

# save this plot
pdf(file = "LogLambdaPlot_simulation.pdf", height = 6, width = 10)
plot(mod$log_lambda, mod$dev.ratio, type = "l",
     xlab = "Log lambda value", ylab = "Deviance ratio", 
     xlim = rev(range(mod$log_lambda)), ylim = c(0.2, max(mod$dev.ratio)))
text(df$log_lambda, 0.2, labels = df$Variables, srt = 90, adj = 0)
dev.off()


### 2) Comparison of AIC and BIC over all different model combinations suggested by lasso, and select the one with the best model fit

## Loop over each step of the lasso, performing a standard linear model (or GLM if binary outcome) and storing the AIC and BIC values
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 8, nrow = 0))
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
df

# Make a var to show number of steps where variables added, and rename the high_sep variable to blank (as is included by default)
df$steps <- 1:nrow(df)
df$Variables[df$steps == 1] <- ""
df


## Now run all combinations of the model variables in a standard LM, and store AIC and BIC values
for (i in 1:nrow(df)) {
  
  vars_temp <- strsplit(df$model_vars[i], " ")[[1]] # Split the variables at each stage of the lasso
    x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp] # Matrix of the covariates at each step of lasso
  mod_temp <- lm(bmi ~ x_hypos_new) # Run the model
  #print(summary(mod_temp))
  
  df$aic[i] <- AIC(mod_temp)
  df$bic[i] <- BIC(mod_temp)
  
}

# Check the dataframe, and select the models with the lowest AIC and BIC values
df

df$model_vars[which.min(df$aic)]
df$model_vars[which.min(df$bic)]


### The best-fitting model, according to both AIC and BIC, is the true model - i.e., that with 'high_sep' as default, plus 'crit3' and 'int3'


## Run selective inference on the best-fitting lasso model (this is the same for both AIC and BIC)

# select the chosen lambda value
lambda <- as.numeric(df$Lambda[which.min(df$aic)])
lambda

# Extract the coefficients (have to exclude the intercept)
beta <- coef(mod, s = lambda, x = x_hypos, y = bmi, 
             penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))), exact = TRUE)[-1]

# Perform selective inference
si_res <- fixedLassoInf(x_hypos, bmi, beta, lambda, alpha = 0.05)
si_res
si_res$vars

# Save results in nicer data frame
si_res <- as.data.frame(cbind(si_res$vars, si_res$coef0, si_res$sd, si_res$pv, si_res$ci))
colnames(si_res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res[,-1] <- lapply(si_res[,-1], function(x) round(as.numeric(as.character(x)), 4))
si_res

## And compare against standard regression results - Results are practically identical, but p-values and CIs are usually slightly wider/corrected in the selective inference method as this accounts for the multiple comparisons (although here they are practically identical). Another benefit of the glmnet method is that we get to see the confounder/covariate results (which are partialled-out/implicitly adjusted for when using LARS and the FWL method).
summary(lm(bmi ~ high_sep + crit3 + high_sep:crit3))
confint(lm(bmi ~ high_sep + crit3 + high_sep:crit3))


### One benefit of this approach over LARS and the FWL method (in addition to being easier to implement and the regression parameters being on the original scale), is that it is possible to edit the penalty factor when running selective inference, so if a variable does not appear in the final lasso model we can constrain it to be in the reported model - particularly handy if lasso selects an interaction term but not a main effect, but we want to report a main effect in the reported model. This functionality is not possible with LARs and the FWL method, as the confounders/covariates are fixed throughout the process.

## Here's an example - Say our best-fitting model from the lasso algorithm was the fourth model, with high_sep, crit3, int3 and green_dec23_int (as potentially suggested by the visual inspection plot). In this example, we may want to force the 'green_dec23' main effect into the final model.
df

# select the chosen lambda value (Note: Sometimes the behaviour of this selective inference package is quite strange, as have to specify the lambda value manually to get the script to work - More on this below.)
lambda <- as.numeric(df$Lambda[df$VarNum == 4])
lambda
lambda <- 0.08160667

# Extract the coefficients (have to exclude the intercept) - Here we edit the penalty term here so that 'green_dec23' (the 16th variable in x_hypos) is forced into the model (also have to force 'green_dec23_int' - variable 17 - into the model here as well).
beta <- coef(mod, s = lambda, x = x_hypos, y = bmi, 
             penalty.factor = (c(0, rep(1, 14), 0, 0)), exact = TRUE)[-1]

# Perform selective inference. Note: There are severe issues with model fit here, any there are lots of warning messages from the 'fixedLassoInf' command, the confidence intervals contain infinite values and the p-values are non-sensical (either 0 or 1).
si_res <- fixedLassoInf(x_hypos, bmi, beta, lambda, alpha = 0.05)
si_res

## The output suggested lowering the glmnet tolerance factor to help improve convergence, so will do that here (Which seems to work and produces more sensible outputs)

# glmnet with different threshold
mod2 <- glmnet(x_hypos, bmi, alpha = 1, thresh = 1E-5, 
               penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))

# Check get same parameters as in previous model (we do)
coef(mod2, s = max(mod2$lambda[mod2$df == 4])); min(mod2$dev.ratio[mod2$df == 4])

# Store lambda value
mod2
lambda <- 0.09830

# Extract betas (have to specify all the parameters to keep here)
beta <- coef(mod2, s = lambda, x = x_hypos, y = bmi, 
              penalty.factor = (c(0, rep(1, 4), 0, 0, rep(1, 8), 0, 0)), exact = TRUE)[-1]

# Perform selective inference (now get more sensible results, although one side of the confidence intervals is much further away from the coefficient than the other side, which is quite odd)
si_res <- fixedLassoInf(x_hypos, bmi, beta, lambda, alpha = 0.05)
si_res
si_res$vars

# Save results in nicer data frame
si_res <- as.data.frame(cbind(si_res$vars, si_res$coef0, si_res$sd, si_res$pv, si_res$ci))
colnames(si_res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res[,-1] <- lapply(si_res[,-1], function(x) round(as.numeric(as.character(x)), 4))
si_res

## And compare against standard regression results - Results similar, but wider (and aymmetrical) confidence intervals for lasso model, so, does appear that none of the extra variables, other than high_sep, crit3 and int3 are associated with the outcome
summary(lm(bmi ~ high_sep + crit3 + high_sep:crit3 + green_dec23 + green_dec23_int))
confint(lm(bmi ~ high_sep + crit3 + high_sep:crit3 + green_dec23 + green_dec23_int))



### 3) Cross-validated lasso to find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model (i.e., lowest MSE). Will use the default number of k-folds (10).
set.seed(3930)
mod.cv <- cv.glmnet(x_hypos, bmi, alpha = 1, penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
mod.cv

plot(mod.cv) # Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters runs along the top of the plot)

# The 1SE model contains 3 parameters (SEP confounder, crit3 and int3 term), so identifies the correct model
coef(mod.cv, s = mod.cv$lambda.1se)

# The model with the lowest MSE contains 12 parameters, so does not correct identify the correct model (but as this optimal lasso is mainly for prediction, this increase in complexity is to be expected as the model hooks on to random noise in the data) - It does contain the two 'true' parameters - 'crit3' and 'int3' - however.
coef(mod.cv, s = mod.cv$lambda.min)

# Save this plot
pdf(file = "CVLasso_simulation.pdf", height = 5, width = 8)
plot(mod.cv)
dev.off()


### Perform selective inference on these models

## Minimum MSE/prediction error model first
# select the chosen lambda value
lambda <- mod.cv$lambda.min
lambda

# Extract the coefficients (have to exclude the intercept)
beta <- coef(mod.cv, s = lambda, x = x_hypos, y = bmi, 
             penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))), exact = TRUE)[-1]


### As an aside, am noticing some very strange behaviour here. If I select the lambda value from the cross-validated model then added this to the 'beta' command, I don't get an error if I exclude the 'penalty.factor' argument, and hence can't force some variables into the model. But if I assign the lambda value manually, the penalty.factor argument works correctly. As 'lambda' is just a standard number, I don't know why this would differ... Compare:

## Standard lambda assignment without the penalty.factor argument in 'beta' - Works fine (which it shouldn't do)
lambda <- mod.cv$lambda.min
beta <- coef(mod.cv, s = lambda, x = x_hypos, y = bmi, exact = TRUE)[-1]

# Manually assigning the lambda value - Now get an error message, and have to add back in the penalty.factor argument to get it to work. Very strange! To be on the safe side, it's probably best to assign lambda manually.
lambda
lambda <- 0.004157
beta <- coef(mod.cv, s = lambda, x = x_hypos, y = bmi, exact = TRUE)[-1]

beta <- coef(mod.cv, s = lambda, x = x_hypos, y = bmi, 
             penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))), exact = TRUE)[-1]


# Perform selective inference - There are some warnings here, and the p-values and confidence intervals have not been estimated (possibly because the model is too complex with 12 parameters). Unlike standard glmnet, can't add a 'thresh' parameter to try and get around these convergence issues.
si_res_CV_minMSE <- fixedLassoInf(x_hypos, bmi, beta, lambda, alpha = 0.05)
si_res_CV_minMSE
si_res_CV_minMSE$vars

# Save results in nicer data frame
si_res_CV_minMSE <- as.data.frame(cbind(si_res_CV_minMSE$vars, si_res_CV_minMSE$coef0, si_res_CV_minMSE$sd, 
                                        si_res_CV_minMSE$pv, si_res_CV_minMSE$ci))
colnames(si_res_CV_minMSE) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res_CV_minMSE[,-1] <- lapply(si_res_CV_minMSE[,-1], function(x) round(as.numeric(as.character(x)), 4))
si_res_CV_minMSE

## And compare against standard regression results (although impossible to really compare as no CIs or p-values for selective inference results. Still, does appear that none of the extra variables, other than high_sep, crit3 and int3 are associated with the outcome)
summary(lm(bmi ~ high_sep + crit3 + high_sep:crit3 + int1 + crit2 + green_inc12 + green_dec12 + green_dec12_int +
             green_inc23 + green_inc23_int + green_dec23))
confint(lm(bmi ~ high_sep + crit3 + high_sep:crit3 + int1 + crit2 + green_inc12 + green_dec12 + green_dec12_int +
             green_inc23 + green_inc23_int + green_dec23))


## Now the model within 1 SE of the minimum MSE model
# select the chosen lambda value
(lambda <- mod.cv$lambda.1se)
lambda <- 0.0983

# Extract the coefficients (have to exclude the intercept)
beta <- coef(mod.cv, s = lambda, x = x_hypos, y = bmi, 
             penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))), exact = TRUE)[-1]

# Perform selective inference - As this model identified crit3 and int3, this model is the same as that reported above for AIC/BIC
si_res_CV_1SE <- fixedLassoInf(x_hypos, bmi, beta, lambda, alpha = 0.05)
si_res_CV_1SE
si_res_CV_1SE$vars

# Save results in nicer data frame
si_res_CV_1SE <- as.data.frame(cbind(si_res_CV_1SE$vars, si_res_CV_1SE$coef0, si_res_CV_1SE$sd, 
                                     si_res_CV_1SE$pv, si_res_CV_1SE$ci))
colnames(si_res_CV_1SE) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res_CV_1SE[,-1] <- lapply(si_res_CV_1SE[,-1], function(x) round(as.numeric(as.character(x)), 4))
si_res_CV_1SE


#### Summary: All of the methods seem to match up well, give consistent answers, and identify the true model correctly (other than the cross-validated lasso with the lowest mean-squared error).



### Quick demonstration of the LARS method using Frish-Waugh-Lovell theorem approach to remove confounding, to show how this can be done (Mostly taken/adapted from this gitHub page: https://github.com/thedunnlab/SLCMA-pipeline/blob/master/SLCMA_completecase_analysis.R)

# Specify confounders, and convert to factor if necessary
covars <- x_hypos[, "high_sep"]
covars <- as.factor(covars)

# Drop SEP confounder from matrix of encoded hypotheses
x_hypos_lars <- x_hypos[, !colnames(x_hypos) == "high_sep"]

# Get the residuals from the exposures/encoded hypotheses, adjusting for covariates
x_resid <-lm(x_hypos_lars ~ covars)$resid
head(x_resid)
summary(x_resid)

# Now get residuals for outcome on covariates
bmi_resid <- lm(bmi ~ covars)$resid
head(bmi_resid)
summary(bmi_resid)

## Now apply the LARS method
lars <- lars(x_resid, bmi_resid)
lars

# Name the first few variables selected
attributes(lars$actions[[1]])$name
attributes(lars$actions[[2]])$name
attributes(lars$actions[[3]])$name
attributes(lars$actions[[4]])$name

## Make a nicer dataframe for these results
last.action <- length(lars$actions)
variables <- numeric(last.action)
selection <- vector("list", last.action)

for(action in 1:last.action) {
  variables[action] <- sum(lars$beta[as.character(action), ] != 0)
  selection[[action]] <- attributes(which(lars$beta[as.character(action), ] != 0))$names
}

additions <- character(last.action)
for(action in 2:last.action) {
  current_selection <- lars$beta[as.character(action-1), ] != 0
  new_selection <- lars$beta[as.character(action), ] != 0
  if(variables[action] > variables[action-1]) {
    additions[action] <-
      dimnames(x_resid)[[2]][new_selection != current_selection]
  }
}

additions[1] <- selection[1]

## 'additions' shows the new variables added at each time point, while 'selection' shows the full model selected at each time point (the 'additions' object doesn't seem to say if/when variables are dropped, though.)
additions
selection

# Elbow plot
par(mar=c(5,4,4,5)+0.1)
plot(c(0, variables), lars$R2, type='l', # lars$R2 - elbow plots model R2 value
     xlab="Variables selected", ylab="R-squared")
text(c(0, variables), rep(0, max(variables)+1),
     labels=c("", additions), srt=90, adj=0)

## Some of the encoded variables are repeated twice in the plot (as included, then dropped, then re-included again), but results do point to the same picture as with the glmnet/lasso approach above, as R2 increases when the first two variables are added - crit3 and int3 - after which there is no noticeable improvement.


## Post-selection inference

# Need to take the sum of squares and re-normalise the results, so they are on the original scale (as currently are using residuals via the FWL method)
sumsq <- lars$normx
x_normed <- scale(x_resid, scale = sumsq)
head(x_normed)
summary(x_normed)

# creating output
step <- 2 # This is the number of encoded variables we want in our final model (based on the elbow plot)
fli <- fixedLassoInf(x_normed, bmi_resid, lars$beta[step + 1,], lars$lambda[step + 1], alpha = 0.05)
scale <- sumsq[fli$vars]
fli$coef0 <- fli$coef0/scale
fli$sd <- fli$sd/scale
fli$ci <- fli$ci/cbind(scale,scale)
fli$vars <- names(fli$vars)

# Print the results on the original scale and save into dataframe
fli

fli.res <- as.data.frame(cbind(fli$vars, fli$coef0, fli$sd, fli$pv, fli$ci))
colnames(fli.res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

fli.res[,-1] <- lapply(fli.res[,-1], function(x) round(as.numeric(as.character(x)), 4))
fli.res

# And compare against standard regression results - Results are practically identical, but p-values and CIs are usually slightly wider/corrected in the selective inference method as this accounts for the multiple comparisons. Another benefit of the glmnet method is that we get to see the confounder/covariate results (which are 'partialled-out'/implicitly conditioned on when using LARS and the FWL method).
summary(lm(bmi ~ high_sep + crit3 + high_sep:crit3))
confint(lm(bmi ~ high_sep + crit3 + high_sep:crit3))




###################################################################################################################
###################################################################################################################
#### Example life course model with confounders and interactions, but with a binary outcome

### Create the binary 'overweight' outcome - For simplicity here, will just recode the 'bmi' variable into overweight or not, based on a BMI > 25 or not.
overweight <- ifelse(bmi > 25, 1, 0)
table(overweight)

# Check the 'true' model, which is SEP as confounder, interaction with green 3, and main effect of green 3. green1 and green2 should also be null. Yup, model works as expected and interaction model better fit than non-interaction model
summary(glm(overweight ~ high_sep + green1 + green2 + green3, family = "binomial"))
summary(glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green3, family = "binomial"))
exp(cbind(coef(summary(glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green3, family = "binomial")))[, 1], confint(glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green3, family = "binomial"))))
anova(glm(overweight ~ high_sep + green1 + green2 + green3, family = "binomial"), 
      glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green3, family = "binomial"), test = "Chisq")

AIC(glm(overweight ~ high_sep + green1 + green2 + green3, family = "binomial")); AIC(glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green3, family = "binomial"))
BIC(glm(overweight ~ high_sep + green1 + green2 + green3, family = "binomial")); BIC(glm(overweight ~ high_sep + green1 + green2 + green3 + high_sep*green3, family = "binomial"))



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

# The first variable included was green3, and then the second variable was 'green_inc23', while the third variable included was int3'. As 'green_inc23' does not cause the outcome, this is incorrect, meaning that this method does not appear to have identified the true model. This is likely because we are using a binary outcome, which reduces the power of the analysis to detect the correct model (see the more formal simulation study for additional details and verification).


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
    change_drop <- setdiff(old_covars, new_covars) # Make sure no covariates dropped at same time
    if (!identical(change_drop, character(0))) { # If covar also dropped, then combine with 'change'
      change_drop <- paste(change_drop, collapse = " ") # Combine variable together, if > 1
      change <- paste0(change, " ", change_drop, " (-)")
    }
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
    change_add <- setdiff(new_covars, old_covars) # Make sure no covariates added at same time
    if (!identical(change_add, character(0))) { # If covar also dropped, then combine with 'change'
      change_add <- paste(change_add, collapse = " ") # Combine variable together, if > 1
      change <- paste0(change, " ", change_add, " (+)")
    }
    dev_diff <- round((new_deviance - old_deviance) * 100, 3) # Diff in deviance between current and previous lambda
    new_dev <- round(new_deviance * 100, 3) # Current deviance value
    temp <- cbind(change, new_dev, dev_diff, new_varNum, new_lambda) # Combine values together
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

# Make a var to show number of steps where variables added, and rename the high_sep variable to remove the '+'
df$steps <- 1:nrow(df)
df$Variables[df$steps == 1] <- "high_sep"
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

# This looks better, as now clearer that 'crit3' was entered first, followed by 'green_inc23' and then 'int3' shortly after. Also need to be aware of the deviance ratio scale (as mentioned above). Here, once 'crit3' is added the deviance ratio increases a little (~1.5%), after which 'green_inc23' and 'int3' are added, and the deviance ratio increases by around another ~0.5% until the next variable is entered in the model (green_dec12), after which the deviance ratio does not increase by very much. This suggests that 'crit3', 'int3' and 'green_inc23' in combination explain most of the variation in the outcome BMI attributable to the life-course hypotheses, and that these variables are associated with the outcome. However, unlike with the continuous outcome example above, the deviance ratios/variance explained by these variables is quite a bit weaker here with this binary outcome, and the model does not identify the true combination of parameters that we simulated.
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



### 2) Comparison of AIC and BIC over all different model combinations suggested by lasso, and select the one with the best model fit

## Loop over each step of the lasso, performing a standard linear model (or GLM if binary outcome) and storing the AIC and BIC values
old_covars <- ""
old_deviance <- 0
old_varNum <- 0
old_lambda <- NA
df <- data.frame(matrix(ncol = 8, nrow = 0))
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
df

# Make a var to show number of steps where variables added, and rename the high_sep variable to blank (as is included by default)
df$steps <- 1:nrow(df)
df$Variables[df$steps == 1] <- ""
df


## Now run all combinations of the model variables in a standard GLM, and store AIC and BIC values
for (i in 1:nrow(df)) {
  
  vars_temp <- strsplit(df$model_vars[i], " ")[[1]] # Split the variables at each stage of the lasso
  x_hypos_new <- x_hypos[, colnames(x_hypos) %in% vars_temp] # Matrix of the covariates at each step of lasso
  mod_temp <- glm(overweight ~ x_hypos_new, family = "binomial") # Run the model
  #print(summary(mod_temp))
  
  df$aic[i] <- AIC(mod_temp)
  df$bic[i] <- BIC(mod_temp)
  
}

# Check the dataframe, and select the models with the lowest AIC and BIC values
df

df$model_vars[which.min(df$aic)]
df$model_vars[which.min(df$bic)]


### The best-fitting model, according to both AIC and BIC, is the not true model, as it includes 'high_sep' (as default), plus 'crit3' and 'int3' (as per the true model), but also includes the additional variable 'green_inc23' (as per the visual inspection method above). These models do still include the correct interaction term ('int3') in the final model, however.


## Run selective inference on the best-fitting lasso model (this is the same for both AIC and BIC)

# select the chosen lambda value
(lambda <- as.numeric(df$Lambda[which.min(df$aic)]))
lambda <- 0.010055

# Extract the coefficients (note: have to *include* the intercept for logistic models).
beta <- coef(mod.binary, s = lambda, x = x_hypos, y = overweight, 
             penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))), exact = TRUE)

# Perform selective inference
si_res <- fixedLassoInf(x_hypos, overweight, beta, lambda, alpha = 0.05, family = "binomial")
si_res
si_res$vars

# Save results in nicer data frame
si_res <- as.data.frame(cbind(si_res$vars, si_res$coef0, si_res$sd, si_res$pv, si_res$ci))
colnames(si_res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res[,-1] <- lapply(si_res[,-1], function(x) round(as.numeric(as.character(x)), 4))
si_res

## And compare against standard regression results - The coefficients are almost the same, but the confidence intervals are much wider, and the p-values much larger, in the selective inference model, to the extent that the interaction term between SEP and critical period 5 is now only 'borderline significant' (and the CIs cross the null slightly) in the selective inference model, although the overall pattern of results is comparable.
summary(glm(overweight ~ high_sep + crit3 + high_sep:crit3 + green_inc23, family = "binomial"))
confint(glm(overweight ~ high_sep + crit3 + high_sep:crit3 + green_inc23, family = "binomial"))



### 3) Cross-validated lasso to find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model. Will use the default number of k-folds (10).
set.seed(3931)
mod.cv.binary <- cv.glmnet(x_hypos, overweight, alpha = 1, family = "binomial", 
                    penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))))
mod.cv.binary

plot(mod.cv.binary) # Plot the log lambda by MSE to show both 'optimal' and '1SE' models (number of parameters runs along the top of the plot)

# Here, this 1SE model contains only 2 parameters (SEP confounder and crit3), so does not identify the correct model as it excludes the 'int3' interaction term. This is likely due to reduced power with using a binary outcome, plus potentially the 1SE cross-validated lasso model being too conservative (see the full simulation study for a more detailed exploration)
coef(mod.cv.binary, s = mod.cv.binary$lambda.1se)

# The 'optimal' model contains 9 parameters, so does not correct identify the correct model (but as this optimal lasso is mainly for prediction, this increase in complexity is to be expected as the model hooks on to random noise in the data). This does contain the true parameters, though.
coef(mod.cv.binary, s = mod.cv.binary$lambda.min)

# Save this plot
pdf(file = "CVLasso_binaryOutcome_simulation.pdf", height = 7, width = 11)
plot(mod.cv.binary)
dev.off()


### Run selective inference on the best-fitting cross-validated lasso models (this differs for minimum MSE and 1SE models)

## Minimum MSE/prediction error
# select the chosen lambda value
(lambda <- mod.cv.binary$lambda.min)
lambda <- 0.002269

# Extract the coefficients (note: have to *include* the intercept for logistic models)
beta <- coef(mod.cv.binary, s = lambda, x = x_hypos, y = overweight, 
             penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))), exact = TRUE)

# Perform selective inference (as lots of parameters, get error message and non-sensical results - Lots of 'infinite' values)
si_res_CV_minMSE <- fixedLassoInf(x_hypos, overweight, beta, lambda, alpha = 0.05, family = "binomial")
si_res_CV_minMSE
si_res_CV_minMSE$vars

# Save results in nicer data frame
si_res_CV_minMSE <- as.data.frame(cbind(si_res_CV_minMSE$vars, si_res_CV_minMSE$coef0, si_res_CV_minMSE$sd, 
                                        si_res_CV_minMSE$pv, si_res_CV_minMSE$ci))
colnames(si_res_CV_minMSE) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res_CV_minMSE[,-1] <- lapply(si_res_CV_minMSE[,-1], function(x) round(as.numeric(as.character(x)), 4))
si_res_CV_minMSE


## Model within 1 SE of the minimum MSE model
# select the chosen lambda value
(lambda <- mod.cv.binary$lambda.1se)
lambda <- 0.0337

# Extract the coefficients (note: have to *include* the intercept for logistic models)
beta <- coef(mod.cv.binary, s = lambda, x = x_hypos, y = overweight, 
             penalty.factor = (c(0, rep(1, ncol(x_hypos) - 1))), exact = TRUE)

# Perform selective inference
si_res_CV_1SE <- fixedLassoInf(x_hypos, overweight, beta, lambda, alpha = 0.05, family = "binomial")
si_res_CV_1SE
si_res_CV_1SE$vars

# Save results in nicer data frame
si_res_CV_1SE <- as.data.frame(cbind(si_res_CV_1SE$vars, si_res_CV_1SE$coef0, si_res_CV_1SE$sd, 
                                     si_res_CV_1SE$pv, si_res_CV_1SE$ci))
colnames(si_res_CV_1SE) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

si_res_CV_1SE[,-1] <- lapply(si_res_CV_1SE[,-1], function(x) round(as.numeric(as.character(x)), 4))
si_res_CV_1SE



#### Summary: When using a binary outcome, the methods give slightly less consistent answers compared to when using continuous outcomes, and do not always identify the true model. For instance, although the visual inspection and AIC/BIC methods both included the correct variables 'crit3' and 'int3', they erroneously included an additional variable. On the other hand, the 1SE cross-validated lasso model only identified the 'crit3' term in the final model; interpretation therefore should be made on a qualitative judgement taking information from all these methods into consideration, rather than just one. As with the continuous outcome, the cross-validated lasso with the lowest MSE selects too many variables.


