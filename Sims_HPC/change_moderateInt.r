### Simulation study script for structured lifecourse models with confounders and interactions

### Created 19/8/2022 by Dan Major-Smith
### R version 4.1.0

### Aim: To test the power of this structured life course approach with interactions to select the 'true' interaction term under varying conditions. Will explore all combinations of the various parameters:
#   - Sample sizes: 1,000 vs 10,000
#   - Exposures: Binary (access to green space; yes/no) vs continuous (distance to green space), and centered vs uncentered. Also want to vary the correlation between exposures to see how collinearity impacts the power of the lasso to detect the true model
#   - Outcome: Binary (overweight/obese) vs continuous (BMI)

## For these simulations, will use the same set-up as in the example simulation script, with 'access to green space' as the exposure measured at three time points, cardiometabolic health as the outcome (BMI/obesity), and SEP as a confounder/interaction term. SEP causes access to green space and the outcome (lower BMI/obesity if higher SEP), while the interaction between SEP and the first green space time-point also causes the outcome.

## In addition to this scenario, we will also vary the strength of the interaction term, to explore how this impacts the power to detect the interaction, as well as varying the specific life course interaction (i.e., the main model will explore an interaction with the first critical period, while other simulations will examine interactions with accumulation and change, to see whether this impacts conclusions).


### This script is for: BMI caused by change from time 2 to time 3 and interaction with this and SEP interaction size = moderate [2 for binary exposures; 0.02 for continuous exposures])


## Load packages (will assume packages already loaded on HPC)
library(glmnet)

## Working directory
setwd("/user/home/ds16565/LifeCycle/Results/")


###########################################################################################################
#### Actual simulations (to see how these simulation parameters were arrived at, see the 'LifeCycle_GreenSpace_SimulationStudyScript.r' script)

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
  AIC_res_temp_main <- rep(NA, n_sims) # Store if correct main effect in final model
  AIC_res_temp_int <- rep(NA, n_sims) # Store if correct interaction effect in final model
  AIC_res_temp_mainint <- rep(NA, n_sims) # Store if correct main effect and interaction in final model
  AIC_res_temp_mainintextra <- rep(NA, n_sims) # Store if main effect and int, plus extra vars, in final model
  BIC_res_temp <- rep(NA, n_sims)
  BIC_res_temp_main <- rep(NA, n_sims) # Store if correct main effect in final model
  BIC_res_temp_int <- rep(NA, n_sims) # Store if correct interaction effect in final model
  BIC_res_temp_mainint <- rep(NA, n_sims) # Store if correct main effect and interaction in final model
  BIC_res_temp_mainintextra <- rep(NA, n_sims) # Store if main effect and int, plus extra vars, in final model
  CV_1SE_res_temp <- rep(NA, n_sims)
  CV_1SE_res_temp_main <- rep(NA, n_sims) # Store if correct main effect in final model
  CV_1SE_res_temp_int <- rep(NA, n_sims) # Store if correct interaction effect in final model
  CV_1SE_res_temp_mainint <- rep(NA, n_sims) # Store if correct main effect and interaction in final model
  CV_1SE_res_temp_mainintextra <- rep(NA, n_sims) # Store if main effect and int, plus extra vars, in final model
  CV_minMSE_res_temp <- rep(NA, n_sims)
  CV_minMSE_res_temp_main <- rep(NA, n_sims) # Store if correct main effect in final model
  CV_minMSE_res_temp_int <- rep(NA, n_sims) # Store if correct interaction effect in final model
  CV_minMSE_res_temp_mainint <- rep(NA, n_sims) # Store if correct main effect and interaction in final model
  CV_minMSE_res_temp_mainintextra <- rep(NA, n_sims) # Store if main effect and int, plus extra vars, in final model
  
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
		bmi <- 25 + (-4 * high_sep) + (-2 * green_inc23) + (2 * high_sep * green_inc23) + rnorm(n, 0, 3)
	}
	
	if (Exposure == "Cont") {
		bmi <- 25 + (-4 * high_sep) + (-0.02 * green_ch23) + (0.02 * high_sep * green_ch23) + rnorm(n, 0, 3)
	}
	if (Output == TRUE) {
      print("BMI summary stats:")
      print(summary(bmi))
      print("")
    }
    
    # If outcome is binary
    if (Outcome == "Binary") {
      overweight <- ifelse(bmi > 25, 1, 0) # Code BMI to binary 'overweight' variable if > 25
      if (Output == TRUE) {
        print(table(overweight))
        print("")
      }
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
    
    # Store if correct main effect in final model - Split by whether binary or continuous exposures
	if (Exposure == "Binary") {
		AIC_res_temp_main[i] <- ifelse("green_inc23" %in% vars_split_aic == TRUE, 1, 0)
		BIC_res_temp_main[i] <- ifelse("green_inc23" %in% vars_split_bic == TRUE, 1, 0)
	}
    if (Exposure == "Cont") {
		AIC_res_temp_main[i] <- ifelse("green_ch23" %in% vars_split_aic == TRUE, 1, 0)
		BIC_res_temp_main[i] <- ifelse("green_ch23" %in% vars_split_bic == TRUE, 1, 0)
	}
    
    # Store if correct interaction in final model - Split by whether binary or continuous exposures
	if (Exposure == "Binary") {
		AIC_res_temp_int[i] <- ifelse("green_inc23_int" %in% vars_split_aic == TRUE, 1, 0)
		BIC_res_temp_int[i] <- ifelse("green_inc23_int" %in% vars_split_bic == TRUE, 1, 0)
	}
    if (Exposure == "Cont") {
		AIC_res_temp_int[i] <- ifelse("green_ch23_int" %in% vars_split_aic == TRUE, 1, 0)
		BIC_res_temp_int[i] <- ifelse("green_ch23_int" %in% vars_split_bic == TRUE, 1, 0)
	}
    
    # Store if correct main effect and interaction both in final model - Split by whether binary or continuous exposures
	if (Exposure == "Binary") {
		AIC_res_temp_mainint[i] <- ifelse("green_inc23" %in% vars_split_aic == TRUE & "green_inc23_int" %in% vars_split_aic == TRUE, 1, 0)
		BIC_res_temp_mainint[i] <- ifelse("green_inc23" %in% vars_split_bic == TRUE & "green_inc23_int" %in% vars_split_bic == TRUE, 1, 0)
	}
	if (Exposure == "Cont") {
		AIC_res_temp_mainint[i] <- ifelse("green_ch23" %in% vars_split_aic == TRUE & "green_ch23_int" %in% vars_split_aic == TRUE, 1, 0)
		BIC_res_temp_mainint[i] <- ifelse("green_ch23" %in% vars_split_bic == TRUE & "green_ch23_int" %in% vars_split_bic == TRUE, 1, 0)
	}
	
    # Store if correct main effect and interaction both in final model, plus other variables
    AIC_res_temp_mainintextra[i] <- ifelse(AIC_res_temp_mainint[i] == 1 & AIC_res_temp[i] == 0, 1, 0)
    BIC_res_temp_mainintextra[i] <- ifelse(BIC_res_temp_mainint[i] == 1 & BIC_res_temp[i] == 0, 1, 0)
    
    
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
    
    # Store if main effect, interaction, both, and both with additional variables in final model  - Split by whether binary or continuous exposures
	if (Exposure == "Binary") {
		CV_1SE_res_temp_main[i] <- ifelse("green_inc23" %in% cv_1SE_covars == TRUE, 1, 0)
		CV_1SE_res_temp_int[i] <- ifelse("green_inc23_int" %in% cv_1SE_covars == TRUE, 1, 0)
		CV_1SE_res_temp_mainint[i] <- ifelse("green_inc23" %in% cv_1SE_covars == TRUE & "green_inc23_int" %in% cv_1SE_covars == TRUE, 1, 0)
    }
	if (Exposure == "Cont") {
		CV_1SE_res_temp_main[i] <- ifelse("green_ch23" %in% cv_1SE_covars == TRUE, 1, 0)
		CV_1SE_res_temp_int[i] <- ifelse("green_ch23_int" %in% cv_1SE_covars == TRUE, 1, 0)
		CV_1SE_res_temp_mainint[i] <- ifelse("green_ch23" %in% cv_1SE_covars == TRUE & "green_ch23_int" %in% cv_1SE_covars == TRUE, 1, 0)
    }
	CV_1SE_res_temp_mainintextra[i] <- ifelse(CV_1SE_res_temp_mainint[i] == 1 & CV_1SE_res_temp[i] == 0, 1, 0)
	
	
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
    
    # Store if main effect, interaction, both, and both with additional variables in final model  - Split by whether binary or continuous exposures
	if (Exposure == "Binary") {
		CV_minMSE_res_temp_main[i] <- ifelse("green_inc23" %in% cv_minMSE_covars == TRUE, 1, 0)
		CV_minMSE_res_temp_int[i] <- ifelse("green_inc23_int" %in% cv_minMSE_covars == TRUE, 1, 0)
		CV_minMSE_res_temp_mainint[i] <- ifelse("green_inc23" %in% cv_minMSE_covars == TRUE & "green_inc23_int" %in% cv_minMSE_covars == TRUE, 1, 0)
    }
	if (Exposure == "Cont") {
		CV_minMSE_res_temp_main[i] <- ifelse("green_ch23" %in% cv_minMSE_covars == TRUE, 1, 0)
		CV_minMSE_res_temp_int[i] <- ifelse("green_ch23_int" %in% cv_minMSE_covars == TRUE, 1, 0)
		CV_minMSE_res_temp_mainint[i] <- ifelse("green_ch23" %in% cv_minMSE_covars == TRUE & "green_ch23_int" %in% cv_minMSE_covars == TRUE, 1, 0)
    }
	CV_minMSE_res_temp_mainintextra[i] <- ifelse(CV_minMSE_res_temp_mainint[i] == 1 & CV_minMSE_res_temp[i] == 0, 1, 0)
	
  }
  
  # Store the summaries of these results to transfer to the main results table
  res <- data.frame(AIC_propcorrect = round(sum(AIC_res_temp) / n_sims * 100, 2),
                    AIC_maincorrect = round(sum(AIC_res_temp_main) / n_sims * 100, 2),
                    AIC_intcorrect = round(sum(AIC_res_temp_int) / n_sims * 100, 2),
                    AIC_mainintcorrect = round(sum(AIC_res_temp_mainint) / n_sims * 100, 2),
                    AIC_mainintextra = round(sum(AIC_res_temp_mainintextra) / n_sims * 100, 2),
                    BIC_propcorrect = round(sum(BIC_res_temp) / n_sims * 100, 2),
                    BIC_maincorrect = round(sum(BIC_res_temp_main) / n_sims * 100, 2),
                    BIC_intcorrect = round(sum(BIC_res_temp_int) / n_sims * 100, 2),
                    BIC_mainintcorrect = round(sum(BIC_res_temp_mainint) / n_sims * 100, 2),
                    BIC_mainintextra = round(sum(BIC_res_temp_mainintextra) / n_sims * 100, 2),
                    CV_1SE_propcorrect = round(sum(CV_1SE_res_temp) / n_sims * 100, 2),
                    CV_1SE_maincorrect = round(sum(CV_1SE_res_temp_main) / n_sims * 100, 2),
                    CV_1SE_intcorrect = round(sum(CV_1SE_res_temp_int) / n_sims * 100, 2),
                    CV_1SE_mainintcorrect = round(sum(CV_1SE_res_temp_mainint) / n_sims * 100, 2),
                    CV_1SE_mainintextra = round(sum(CV_1SE_res_temp_mainintextra) / n_sims * 100, 2),
                    CV_minMSE_propcorrect = round(sum(CV_minMSE_res_temp) / n_sims * 100, 2),
                    CV_minMSE_maincorrect = round(sum(CV_minMSE_res_temp_main) / n_sims * 100, 2),
                    CV_minMSE_intcorrect = round(sum(CV_minMSE_res_temp_int) / n_sims * 100, 2),
                    CV_minMSE_mainintcorrect = round(sum(CV_minMSE_res_temp_mainint) / n_sims * 100, 2),
                    CV_minMSE_mainintextra = round(sum(CV_minMSE_res_temp_mainintextra) / n_sims * 100, 2))
  return(res)
  
}


## Next, the number of simulations per combination of parameters (1,000), the target hypotheses we simulated are the true model (as this differs in this simulation depending on whether the exposures are binary or continuous, will define this before each simulation run instead), and set up a data frame to store the results in
n_sims <- 1000
set.seed(32876)

#target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
#target_covars <- c("high_sep", "green_ch23", "green_ch23_int")

results <- data.frame(sampleSize = rep(c(1000, 10000), 16),
                      Exposure = rep(c(rep("Binary", 8), rep("Cont", 8)), 2),
                      Centered = rep(c("No", "No", "Yes", "Yes"), 8),
                      Collinear = rep(c(rep("Low", 4), rep("High", 4)), 4),
                      Outcome = c(rep("Cont", 16), rep("Binary", 16)),
                      AIC_propcorrect = rep(NA, 32),
                      AIC_maincorrect = rep(NA, 32),
                      AIC_intcorrect = rep(NA, 32),
                      AIC_mainintcorrect = rep(NA, 32),
                      AIC_mainintextra = rep(NA, 32),
                      BIC_propcorrect = rep(NA, 32),
                      BIC_maincorrect = rep(NA, 32),
                      BIC_intcorrect = rep(NA, 32),
                      BIC_mainintcorrect = rep(NA, 32),
                      BIC_mainintextra = rep(NA, 32),
                      CV_1SE_propcorrect = rep(NA, 32),
                      CV_1SE_maincorrect = rep(NA, 32),
                      CV_1SE_intcorrect = rep(NA, 32),
                      CV_1SE_mainintcorrect = rep(NA, 32),
                      CV_1SE_mainintextra = rep(NA, 32),
                      CV_minMSE_propcorrect = rep(NA, 32),
                      CV_minMSE_maincorrect = rep(NA, 32),
                      CV_minMSE_intcorrect = rep(NA, 32),
                      CV_minMSE_mainintcorrect = rep(NA, 32),
                      CV_minMSE_mainintextra = rep(NA, 32))

#results


### First simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; continuous outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 1
results[k, 6:25] <- res


### Second simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; continuous outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 2
results[k, 6:25] <- res


### Third simulation: Sample size = 1000; binary exposure; centered; low collinearity; continuous outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 3
results[k, 6:25] <- res


### Fourth simulation: Sample size = 10000; binary exposure; centered; low collinearity; continuous outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 4
results[k, 6:25] <- res


### Fifth simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 5
results[k, 6:25] <- res


### Sixth simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; continuous outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 6
results[k, 6:25] <- res


### Seventh simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 7
results[k, 6:25] <- res


### Eighth simulation: Sample size = 10000; binary exposure; centered; High collinearity; continuous outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 8
results[k, 6:25] <- res


### Ninth simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; continuous outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 9
results[k, 6:25] <- res


### Tenth simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; continuous outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 10
results[k, 6:25] <- res


### Eleventh simulation: Sample size = 1000; continuous exposure; centered; low collinearity; continuous outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 11
results[k, 6:25] <- res


### Twelfth simulation: Sample size = 10000; continuous exposure; centered; low collinearity; continuous outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 12
results[k, 6:25] <- res


### 13th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; continuous outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 13
results[k, 6:25] <- res


### 14th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; continuous outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 14
results[k, 6:25] <- res


### 15th simulation: Sample size = 1000; binary exposure; centered; High collinearity; continuous outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 15
results[k, 6:25] <- res


### 16th simulation: Sample size = 10000; continuous exposure; centered; High collinearity; continuous outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Cont", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 16
results[k, 6:25] <- res


### 17th simulation: Sample size = 1000; binary exposure; uncentered; low collinearity; binary outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 17
results[k, 6:25] <- res


### 18th simulation: Sample size = 10000; binary exposure; uncentered; low collinearity; binary outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 18
results[k, 6:25] <- res


### 19th simulation: Sample size = 1000; binary exposure; centered; low collinearity; binary outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 19
results[k, 6:25] <- res


### 20th simulation: Sample size = 10000; binary exposure; centered; low collinearity; binary outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 20
results[k, 6:25] <- res


### 21st simulation: Sample size = 1000; binary exposure; uncentered; High collinearity; continuous outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 21
results[k, 6:25] <- res


### 22nd simulation: Sample size = 10000; binary exposure; uncentered; High collinearity; binary outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 22
results[k, 6:25] <- res


### 23rd simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 23
results[k, 6:25] <- res


### 24th simulation: Sample size = 10000; binary exposure; centered; High collinearity; binary outcome
target_covars <- c("high_sep", "green_inc23", "green_inc23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Binary", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 24
results[k, 6:25] <- res


### 25th simulation: Sample size = 1000; continuous exposure; uncentered; low collinearity; binary outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 25
results[k, 6:25] <- res


### 26th simulation: Sample size = 10000; continuous exposure; uncentered; low collinearity; binary outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "Low",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 26
results[k, 6:25] <- res


### 27th simulation: Sample size = 1000; continuous exposure; centered; low collinearity; binary outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 27
results[k, 6:25] <- res


### 28th simulation: Sample size = 10000; continuous exposure; centered; low collinearity; binary outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "Low",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 28
results[k, 6:25] <- res


### 29th simulation: Sample size = 1000; continuous exposure; uncentered; High collinearity; binary outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 29
results[k, 6:25] <- res


### 30th simulation: Sample size = 10000; continuous exposure; uncentered; High collinearity; binary outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "No", Collinear = "High",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 30
results[k, 6:25] <- res


### 31st simulation: Sample size = 1000; binary exposure; centered; High collinearity; binary outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 1000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 31
results[k, 6:25] <- res


### 32nd simulation: Sample size = 10000; continuous exposure; centered; High collinearity; binary outcome
target_covars <- c("high_sep", "green_ch23", "green_ch23_int")
res <- lasso_sim(n_sims = n_sims, sampleSize = 10000, Exposure = "Cont", Centered = "Yes", Collinear = "High",
                 Outcome = "Binary", Output = FALSE)

# Store the summaries of these results in the main results table
k <- 32
results[k, 6:25] <- res


## Save the results table
results

write.csv(results, file = "simulationResults_change_moderateInt.csv", row.names = FALSE)

