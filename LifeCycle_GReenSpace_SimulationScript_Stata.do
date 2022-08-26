*** Example simulation script for structured lifecourse models with confounders and interactions

*** Created 17/2/2022 by Dan Major-Smith
*** Stata v.17.0


*** Previous work has detailed a structured approach to life course modelling for both binary exposures (Smith et al., 2015: https://journals.lww.com/epidem/Fulltext/2015/09000/Model_Selection_of_the_Effect_of_Binary_Exposures.15.aspx) and continuous exposures with confounding (Smith et al., 2016: https://academic.oup.com/ije/article/45/4/1271/2951966?login=true) using LARS/lasso methods. Building on these approaches, we aim to demonstrate here how interactions can be incorporated into these models. 

** However, we will not be using the LARS/covariance test method here, as: 1) there are issues with the covariance test and it is no longer recommended; and 2) Using the LARS method, it is possible to 'force' continuous confounders to be included in the model, however, this is more difficult for binary/categorical covariates and/or binary outcomes. So here we will focus on using ordinary lasso models via glmnet, rather than the LARS method. Rather than being a stepwise procedure like LARS, this approach gradually increases the lambda value and lets more variables into the model; this can make it difficult to assess when adding a new covariate does or does not improve model fit. Interpretation is therefore more subjective, based on inspection of the improvement in the deviance ratio (a measure of goodness-of-fit similar to R2).

** To interpret these results, we will focus on three approaches:
**  - 1) Using a 'subjective' approach looking at the order in which hypotheses were entered into the model, combined with a plot of deviance/variance explained when each predictor was added, and making judgement based on these sources of information
**  - 2) Using a 'relaxed lasso'-type approach, where use a standard LM/GLM on the model the lasso selects at each step, then comparing model fit of all these models to detect the best-fitting model. Will use both AIC and BIC as measures of model fit.
**  - 3) Using cross-validated lasso and selecting the model within 1 SE of the best-fitting model (this method is similar to ordinary lasso, but uses k-fold cross-validation to identify the best-fitting model while minimising overfitting). Using the model within 1 SE of the best-fitting model, rather than the best-fitting model itself, should also avoid potential overfitting, leaving only variables more strongly associated with the outcome in the final model. Will also explore model with lowest mean-squared error as well, as a comparison.

** If all these methods give a similar answer, then we can have more confidence in the conclusions drawn.


** Note that the Stata lasso command was only introduced in Stata version 16; the scripts below will not work for early versions of Stata (if using an earlier version of Stata, see the user-written packages 'lassopack': https://ideas.repec.org/c/boc/bocode/s458458.html)



****************************************************************************************
**** Simulating life course data with interaction terms

*** Here, will assume exposures and covariates/interaction terms are binary, with a continuous outcome.
** Exposure: Access to green space (within 300 metres) in pregnancy, age 4 and age 7 (0 = no access; 1 = access).
** Confounder/interaction term: Socioeconomic position (0 = low; 1 = high)
** Outcome: BMI age 7 (continuous)

** The model/DAG we're simulating here is that SEP is a confounder, but also interacts with first green space exposure in pregnancy to impact BMI at age 7


*** Generate the data (although note also that the simulated data here is purely to illustrate the logic and application of these structure life-course methods, and should not be taken as a reflection of the real-world patterns or effect sizes of these variables)

* Sample size of ~10,000 (approximate ALSPAC participation early in study)
clear
set seed 5678
set obs 10000

* SEP - Not caused by anything, so just do 50/50 split (1 = high SEP)
gen high_sep = rbinomial(1, 0.5)
tab high_sep

* First green space exposure in pregnancy - Caused by SEP (higher SEP = more exposure to green space; three times the odds) - Approx. 50% near green space
gen green1_p = invlogit(log(0.6) + (log(3) * high_sep))
sum green1_p

gen green1 = rbinomial(1, green1_p)
tab green1

tab high_sep green1

drop green1_p

* Second green space exposure at age 4 - Caused by Green1 and SEP
gen green2_p = invlogit(log(0.3) + (log(3) * high_sep) + (log(5) * green1))
sum green2_p

gen green2 = rbinomial(1, green2_p)
tab green2

tab high_sep green2
tab green1 green2

drop green2_p

* Third green space exposure at age 7 - Caused by Green2 and SEP
gen green3_p = invlogit(log(0.3) + (log(3) * high_sep) + (log(5) * green2))
sum green3_p

gen green3 = rbinomial(1, green3_p)
tab green3

tab high_sep green3
tab green1 green3
tab green2 green3

drop green3_p

* Now for continuous BMI outcome - Caused by SEP (higher SEP = lower BMI), plus interaction with green1 (lower SEP and access to green space = lower BMI compared to lower SEP and no access to green space). Assuming no main effect of green1 here.
gen bmi = 25 + (-4 * high_sep) + (-2 * green1) + (2 * high_sep * green1) + rnormal(0, 3)
sum bmi

** Descriptive stats split by possible combinations of SEP and green1 to check simulation worked
table high_sep green1, stat(mean bmi)

* Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be null. Yup, model works as expected and interaction model better fit than non-interaction model
regress bmi high_sep green1 green2 green3
estat ic
est store A
regress bmi high_sep##green1 green2 green3
estat ic
est store B

lrtest A B

estimates clear


*** Encode the life course hypotheses

* Critical period at first time point only
gen crit1 = green1

* Interaction between SEP and first time point
gen int1 = high_sep * green1

* Critical period at second time point only
gen crit2 = green2

* Interaction between SEP and second time point
gen int2 = high_sep * green2

* Critical period at third time point only
gen crit3 = green3

* Interaction between SEP and third time point
gen int3 = high_sep * green3

* Linear accumulation of all exposures
gen accumulation = green1 + green2 + green3

* Interaction between SEP and cumulative exposure
gen int_accum = high_sep * accumulation

* Increase in access to green space from time 1 to time 2
gen green_inc12 = (1 - green1) * green2

* Increase in access to green space from time 1 to time 2, with an interaction with SEP
gen green_inc12_int = (1 - green1) * green2 * high_sep

* Decrease in access to green space from time 1 to time 2
gen green_dec12 = (1 - green2) * green1

* Decrease in access to green space from time 1 to time 2, with an interaction with SEP
gen green_dec12_int = (1 - green2) * green1 * high_sep

* Increase in access to green space from time 2 to time 3
gen green_inc23 = (1 - green2) * green3

* Increase in access to green space from time 2 to time 3, with an interaction with SEP
gen green_inc23_int = (1 - green2) * green3 * high_sep

* Decrease in access to green space from time 2 to time 3
gen green_dec23 = (1 - green3) * green2

* Decrease in access to green space from time 2 to time 3, with an interaction with SEP
gen green_dec23_int = (1 - green3) * green2 * high_sep


** Check correlation matrix of all these hypotheses (>0.9 would be a cause for concern) - correlation between 'int_accum' and 'int2' is 0.90, but that should be fine here.
cor crit1-green_dec23_int


*** Run the model using standard lasso (i.e., no cross-validation; equivalent to generic 'glmnet' in R). Variables in brackets are constrained to always be included in the model
lasso linear bmi (high_sep) crit1-green_dec23_int, selection(none)
est store sim

** Summarise these results - The first variable included was green1, and then the second variable included was the interaction term between SEP and green1, which is correct, as this is how we coded the data, although 'green_dec12_int' was added shortly after 'int1', so not a clear separation between the two 'true' variables and all others.
lassoknots

lassoselect id = 2
lassocoef, display(coef, penalized) nolegend
lassogof

lassoselect id = 12
lassocoef, display(coef, penalized) nolegend
lassogof

lassoselect id = 15
lassocoef, display(coef, penalized) nolegend
lassogof

* Plot these results
coefpath



*** Will explore a few methods to check whether we get the right result and do not include any other incorrect parameters in the 'best' model

*** 1) Visual inspection

** The lassoknots command gives the steps of the lasso, along with the lambda value, the number of non-zero coefficients, R-squared value, and the variables added/dropped to the model (option 'all' displays all steps, while the default only gives steps when variables were added/removed).
lassoknots
lassoknots, all

* It should be possible to convert this to a plot of (log) lambda by R-squared, with the variables added at the correct position (as done in the R script), but the stored results do not include the variables added/removed, so to convert this to a format ready to create a plot will take a fair bit of manual work which I will not attempt here (plus, the 'lassoknots' tables gives pretty much the same information as the plot).
matrix res = r(table)
matrix list res


*** 2) Comparison of AIC and BIC over all different model combinations suggested by lasso, and select the one with the best model fit

** I am being a bit lazy here and not testing all of the lasso models (like in the R script), but am just showing the concepts.
lassoknots

* First, compare the baseline (confounder only) model to the model including the first term added (crit1)
regress bmi high_sep
estat ic
local mod1_aic = r(S)[1,"AIC"]
local mod1_bic = r(S)[1,"BIC"]

regress bmi high_sep crit1
estat ic
local mod2_aic = r(S)[1,"AIC"]
local mod2_bic = r(S)[1,"BIC"]

* Find strong support for more complex model which includes crit1 main effect
di "Model 1 AIC = " `mod1_aic' "; Model 2 AIC = " `mod2_aic'
di "Model 1 BIC = " `mod1_bic' "; Model 2 BIC = " `mod2_bic'


** Next, compare the model including the next hypotheses added (int1)
est restore sim
lassoknots

regress bmi high_sep crit1 int1
estat ic
local mod3_aic = r(S)[1,"AIC"]
local mod3_bic = r(S)[1,"BIC"]

* Adding the SEP by crit1 interaction term again improves model fit
di "Model 2 AIC = " `mod2_aic' "; Model 3 AIC = " `mod3_aic'
di "Model 2 BIC = " `mod2_bic' "; Model 3 BIC = " `mod3_bic'


** Does adding the next parameter (green_dec12_int) improve model fit?
est restore sim
lassoknots

regress bmi high_sep crit1 int1 green_dec23_int
estat ic
local mod4_aic = r(S)[1,"AIC"]
local mod4_bic = r(S)[1,"BIC"]

* No real improvement in model fit with BIC (as expected), but AIC shows minor improvement in model fit
di "Model 3 AIC = " `mod3_aic' "; Model 4 AIC = " `mod4_aic'
di "Model 3 BIC = " `mod3_bic' "; Model 4 BIC = " `mod4_bic'


** Does adding the next parameter (green_inc23_int) improve model fit?
est restore sim
lassoknots

regress bmi high_sep crit1 int1 green_dec23_int green_inc23_int
estat ic
local mod5_aic = r(S)[1,"AIC"]
local mod5_bic = r(S)[1,"BIC"]

* Now there's no improvement in AIC or BIC
di "Model 4 AIC = " `mod4_aic' "; Model 5 AIC = " `mod5_aic'
di "Model 4 BIC = " `mod4_bic' "; Model 5 BIC = " `mod5_bic'

estimates clear


*** 3) Cross-validated lasso to find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model. Will use the default number of k-folds (10).

** Specify the 1SE rule first
lasso linear bmi (high_sep) crit1-green_dec23_int, selection(cv, folds(10) serule) rseed(3930)

* Display results - Note that this method only includes one additional variable (crit1) on top of high_sep, so does not include 'int1' and therefore is too conservative and does not select the correct model
lassoknots
lassocoef, display(coef, penalized) nolegend
lassogof
lassoinfo

* Plot the cross-validated results (this gives both 1SE and CV min values)
cvplot

** Now for 'optimal' lambda value based on cross-validation function (CV mean prediction error)
lasso linear bmi (high_sep) crit1-green_dec23_int, selection(cv, folds(10)) rseed(3930)

* Display results - This method includes 9 variables in addition to high_sep, so is too liberal and also does not specify the correct model (although crit1 and int1 are included in the final model)
lassoknots
lassocoef, display(coef, penalized) nolegend
lassogof
lassoinfo


**** Summary: The methods seem match up relatively well and give broadly consistent answers, although in this example the 1SE cross-validated lasso was too conservative and did not include the 'int1' interaction term in the final model, while model comparison using AIC (but not BIC) found that inclusion of an additional term beyond 'crit1' and 'int1' marginal improved model fit. The cross-validated lasso with the lowest MSE contained the two true parameters, but also included numerous additional variables as well)



****************************************************************************************
****************************************************************************************
**** Example life course model with confounders and interactions, but with a binary outcome

*** Create the binary 'overweight' outcome - For simplicity here, will just recode the 'bmi' variable into overweight or not, based on a BMI > 25 or not.
gen overweight = 0
replace overweight = 1 if bmi > 25
tab overweight

* Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be null. Yup, model works as expected and interaction model better fit than non-interaction model
logistic overweight high_sep green1 green2 green3
estat ic
est store A

logistic overweight high_sep##green1 green2 green3
estat ic
est store B

lrtest A B


*** Everything else is the same as above, just using a logistic rather than linear lasso model


*** Run the model using standard lasso (i.e., no cross-validation; equivalent to generic 'glmnet' in R). Variables in brackets are constrained to always be included in the model
lasso logit overweight (high_sep) crit1-green_dec23_int, selection(none)
est store sim

** Summarise these results - The first variable included was crit1, and then the second variable included was the interaction term between SEP and green1 ('int1'; which is also correct). The next variable added after was 'green_dec12_int', although the improvement in model fit seems minimal. So it would appear that this method has identified the true model (although given the reduction in power when using binary outcomes, the correct model may not always be identified; see the full simulation results for more details)
lassoknots

lassoselect id = 2
lassocoef, display(coef, penalized eform) nolegend
lassogof

lassoselect id = 18
lassocoef, display(coef, penalized eform) nolegend
lassogof

lassoselect id = 28
lassocoef, display(coef, penalized eform) nolegend
lassogof

* Plot these results
coefpath


*** Will explore a few methods to check whether we get the right result and do not include any other incorrect parameters in the 'best' model

*** 1) Visual inspection

** The lassoknots command gives the steps of the lasso, along with the lambda value, the number of non-zero coefficients, in-sample deviance ratio, and the variables added/dropped to the model (option 'all' displays all steps, while the default only gives steps when variables were added/removed).
lassoknots
lassoknots, all


*** 2) Comparison of AIC and BIC over all different model combinations suggested by lasso, and select the one with the best model fit

** As above, I am being a bit lazy here and not testing all of the lasso models (like in the R script), but am just showing the concepts.
lassoknots

** First, compare the baseline (confounder only) model to the model including the first term added (crit1)
logistic overweight high_sep
estat ic
local mod1_aic = r(S)[1,"AIC"]
local mod1_bic = r(S)[1,"BIC"]

logistic overweight high_sep crit1
estat ic
local mod2_aic = r(S)[1,"AIC"]
local mod2_bic = r(S)[1,"BIC"]

* Find strong support for more complex model which includes crit1 main effect
di "Model 1 AIC = " `mod1_aic' "; Model 2 AIC = " `mod2_aic'
di "Model 1 BIC = " `mod1_bic' "; Model 2 BIC = " `mod2_bic'


** Next, compare the model including the next hypotheses added (int1)
est restore sim
lassoknots

logistic overweight high_sep crit1 int1
estat ic
local mod3_aic = r(S)[1,"AIC"]
local mod3_bic = r(S)[1,"BIC"]

* Adding the SEP by crit1 interaction term again improves model fit
di "Model 2 AIC = " `mod2_aic' "; Model 3 AIC = " `mod3_aic'
di "Model 2 BIC = " `mod2_bic' "; Model 3 BIC = " `mod3_bic'


** Does adding the next parameter (green_dec12_int) improve model fit?
est restore sim
lassoknots

logistic overweight high_sep crit1 int1 green_dec12_int
estat ic
local mod4_aic = r(S)[1,"AIC"]
local mod4_bic = r(S)[1,"BIC"]

* No real improvement in model fit with either AIC or BIC (as expected)
di "Model 3 AIC = " `mod3_aic' "; Model 4 AIC = " `mod4_aic'
di "Model 3 BIC = " `mod3_bic' "; Model 4 BIC = " `mod4_bic'

estimates clear



*** 3) Cross-validated lasso to find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model. Will use the default number of k-folds (10).

** Specify the 1SE rule first
lasso logit overweight (high_sep) crit1-green_dec23_int, selection(cv, folds(10) serule) rseed(3932)

* Display results - Note that this method only includes the variable 'crit1' variable in addition to high_sep, so does not include 'crit1' and therefore is too conservative and does not select the correct model.
lassoknots
lassocoef, display(coef, penalized eform) nolegend
lassogof
lassoinfo

* Plot the cross-validated results (this gives both 1SE and CV min values)
cvplot

** Now for 'optimal' lambda value based on cross-validation function (CV mean prediction error)
lasso logit overweight (high_sep) crit1-green_dec23_int, selection(cv, folds(10)) rseed(3932)

* Display results - This method includes 4 variables in addition to crit1 and int1, so is somewhat too liberal and also does not specify the correct model
lassoknots
lassocoef, display(coef, penalized eform) nolegend
lassogof
lassoinfo


**** Summary: When using a binary outcome, the methods give slightly less consistent answers compared to when using continuous outcomes, and do not always identify the true model (likely due to power issues when using binary outcomes; see for formal simulation study results for a more detailed exploration). For instance, the visual inspection and likelihood ratio test methods both included the correct variables 'crit1' and 'int1', while the 1SE cross-validated lasso model only identified the 'crit1' term in the final model; interpretation therefore should be made on a qualitative judgement taking information from all these methods into consideration, rather than just one. As with the continuous outcome, the cross-validated lasso with the lowest MSE selects too many variables.


