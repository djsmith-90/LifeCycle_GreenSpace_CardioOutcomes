*** Example simulation script for structured lifecourse models with confounders and interactions

*** Created 17/2/2022 by Dan Major-Smith
*** Stata v.17.0


*** Previous work has detailed a structured approach to life course modelling for both binary exposures (Smith et al., 2015: https://journals.lww.com/epidem/Fulltext/2015/09000/Model_Selection_of_the_Effect_of_Binary_Exposures.15.aspx) and continuous exposures with confounding (Smith et al., 2016: https://academic.oup.com/ije/article/45/4/1271/2951966?login=true) using LARS/lasso methods. Building on these approaches, we aim to demonstrate here how interactions can be incorporated into these models. 

** However, we will not be using the LARS/covariance test method here, as: 1) there are issues with the covariance test and it is no longer recommended; and 2) Using the LARS method, it is possible to 'force' continuous confounders to be included in the model, however, this does not appear possible for binary/categorical variables, especially if we want to include them as interaction terms. So here we will focus on using ordinary lasso models via glmnet, rather than the LARS method. Rather than being a stepwise procedure like LARS, this approach gradually increases the lambda value and lets more variables into the model; this can make it difficult to assess when adding a new covariate does or does not improve model fit. Interpretation is therefore more subjective, based on inspection of the improvement in the deviance ratio (a measure of goodness-of-fit similar to R2).

** To interpret these results, we will focus on three approaches:
**  - 1) Using a 'subjective' approach looking at the order in which hypotheses were entered into the model, combined with a plot of deviance/variance explained when each predictor was added, and making judgement based on these sources of information
**  - 2) Taking each variable(s) entered in turn in the lasso model, using a likelihood ratio test to compare standard linear/logistic models with/without the next predictor in. This will provide a formal test as to whether the hypothesis predicts the outcome.
**  - 3) Using cross-validated lasso and selecting the model within 1 SE of the best-fitting model (this method is similar to ordinary lasso, but uses k-fold cross-validation to identify the best-fitting model while minimising overfitting). Using the model within 1 SE of the best-fitting model, rather than the best-fitting model itself, should also avoid potential overfitting, leaving only variables more strongly associated with the outcome in the final model.

** If all these methods give a similar answer, then we can have more confidence in the conclusions drawn.


** Note that the Stata lasso command was only introduced in Stata version 16; the scripts below will not work for early versions of Stata (if using an earlier version of Stata, see the user-written packages 'lassopack': https://ideas.repec.org/c/boc/bocode/s458458.html)



****************************************************************************************
**** Simulating life course data with interaction terms

*** Here, will assume exposures and covariates/interaction terms are binary, with a continuous outcome.
** Exposure: Access to green space (within 300 metres) in pregnancy, age 4 and age 7 (0 = no access; 1 = access).
** Confounder/interaction term: Socioeconomic position (0 = low; 1 = high)
** Outcome: BMI age 7 (continuous)

** The model/DAG we're simulating here is that SEP is a confounder, but also interacts with first green space exposure in pregnancy to impact BMI at age 7


*** Generate the data (although note also that the simulated data here is purely to illustrate the logic and application of these structure life-course methods, and should not in any way be taken as a reflection of the real-world patterns or effect sizes of these variables)

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
est store A
regress bmi high_sep##green1 green2 green3
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


*** 2) Likelihood ratio tests at inclusion of each new parameter (if only one parameter is added at a time we don't really we don't need to do an LR test, as for single parameters the p-value from the model will be identical to those of the LR test, but we're doing this LR test because sometimes multiple terms get added to the lasso at one time-point)

** First, compare the baseline (confounder only) model to the model including the first term added (crit1)
lassoknots

regress bmi high_sep
est store base

regress bmi high_sep crit1
est store param1

* Find strong support for more complex model which includes crit1 main effect
lrtest base param1


** Next, compare the model including the next hypotheses added (int1)
est restore sim
lassoknots

regress bmi high_sep crit1 int1
est store param2

* Adding the SEP by crit1 interaction term again improves model fit
lrtest param1 param2


** Does adding the next parameter (green_dec12_int) improve model fit?
est restore sim
lassoknots

regress bmi high_sep crit1 int1 green_dec23_int
est store param3

* No real improvement in model fit (as expected)
lrtest param2 param3

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


**** Summary: The methods seem match up relatively well and give broadly consistent answers, although in this example the 1SE cross-validated lasso was too conservative and did not include the 'int1' interaction term in the final model.



****************************************************************************************
****************************************************************************************
**** Example life course model with confounders and interactions, but with a binary outcome

*** Create the binary 'overweight' outcome - Will use the same model as for BMI - Overweight caused by SEP (higher SEP = lower probability of being overweight), plus interaction with green1 (lower SEP and access to green space = lower chances of being overweight compared to lower SEP and no access to green space). Assuming no main effect of green1 here.
set seed 3930

gen overweight = .
replace overweight = rbinomial(1, 0.1) if high_sep == 1 & green1 == 1
replace overweight = rbinomial(1, 0.1) if high_sep == 1 & green1 == 0
replace overweight = rbinomial(1, 0.3) if high_sep == 0 & green1 == 1
replace overweight = rbinomial(1, 0.5) if high_sep == 0 & green1 == 0
tab overweight

tab overweight green1
tab overweight high_sep

** descriptive stats split by possible combinations of SEP and green1 to check simulation worked
table high_sep green1, stat(mean overweight)

* Check the 'true' model, which is SEP as confounder, interaction with green 1, and main effect of green 1. green2 and green3 should also be null. Yup, model works as expected and interaction model better fit than non-interaction model
logistic overweight high_sep green1 green2 green3
est store A

logistic overweight high_sep##green1 green2 green3
est store B

lrtest A B


*** Everything else is the same as above, just using a logistic rather than linear lasso model


*** Run the model using standard lasso (i.e., no cross-validation; equivalent to generic 'glmnet' in R). Variables in brackets are constrained to always be included in the model
lasso logit overweight (high_sep) crit1-green_dec23_int, selection(none)
est store sim

** Summarise these results - The first variable included was green1, but then the second variables included were the interaction term between SEP and green1 (which is correct), but also the 'assumulation' hypothesis (which is incorrect).
lassoknots

lassoselect id = 2
lassocoef, display(coef, penalized eform) nolegend
lassogof

lassoselect id = 15
lassocoef, display(coef, penalized eform) nolegend
lassogof

lassoselect id = 37
lassocoef, display(coef, penalized eform) nolegend
lassogof

* Plot these results
coefpath


*** Will explore a few methods to check whether we get the right result and do not include any other incorrect parameters in the 'best' model

*** 1) Visual inspection

** The lassoknots command gives the steps of the lasso, along with the lambda value, the number of non-zero coefficients, R-squared value, and the variables added/dropped to the model (option 'all' displays all steps, while the default only gives steps when variables were added/removed).
lassoknots
lassoknots, all


*** 2) Likelihood ratio tests at inclusion of each new parameter (if only one parameter is added at a time we don't really we don't need to do an LR test, as for single parameters the p-value from the model will be identical to those of the LR test, but we're doing this LR test because sometimes multiple terms get added to the lasso at one time-point)

** First, compare the baseline (confounder only) model to the model including the first term added (crit1)
lassoknots

logistic overweight high_sep
est store base

logistic overweight high_sep crit1
est store param1

* Find strong support for more complex model which includes crit1 main effect
lrtest base param1


** Next, compare the model including the next hypotheses added (int1 and accumulation)
est restore sim
lassoknots

logistic overweight high_sep crit1 int1 accumulation
est store param2

* Adding the SEP by crit1 interaction term and accumulation terms again improves model fit (however, the 'accumlation' term was not strongly associated with the outcome, while the 'int1' interaction term was)
lrtest param1 param2


** Does adding the next parameters (int_accum, green_inc12 and green_inc23) improve model fit?
est restore sim
lassoknots

logistic overweight high_sep crit1 int1 accumulation int_accum green_inc12 green_inc23
est store param3

* No real improvement in model fit
lrtest param2 param3

estimates clear


*** 3) Cross-validated lasso to find 'optimal' model for out-of-sample prediction. Will compare both 'optimal' and '1 SE' models, although the 1 SE model is probably better to avoid overfitting as it selects the best model within 1 SE of the 'optimal' model. Will use the default number of k-folds (10).

** Specify the 1SE rule first
lasso logit overweight (high_sep) crit1-green_dec23_int, selection(cv, folds(10) serule) rseed(3932)

* Display results - Note that this method does not include any additional variable on top of high_sep, so does not include 'crit1' or 'int1' and therefore is far too conservative and does not select the correct model.
lassoknots
lassocoef, display(coef, penalized eform) nolegend
lassogof
lassoinfo

* Plot the cross-validated results (this gives both 1SE and CV min values)
cvplot

** Now for 'optimal' lambda value based on cross-validation function (CV mean prediction error)
lasso logit overweight (high_sep) crit1-green_dec23_int, selection(cv, folds(10)) rseed(3932)

* Display results - This method includes 4 variables in addition to high_sep (crit1, int1 and accumulation), so is somewhat too liberal and also does not specify the correct model (although it is closer to the 'truth' than the 1SE model)
lassoknots
lassocoef, display(coef, penalized eform) nolegend
lassogof
lassoinfo


**** Summary: The methods seem to broadly give similar answers, although none of the methods identifies the true model correctly. These lassos incorrectly include both 'int1' and 'accumulation' at the same time, while the 1SE cross-validated lasso does not include any hypotheses beyond the 'high_sep' confounder/covariate. This highlights that we need to use multiple methods to 'triangulate' results which are consistent with each method, and that interpretation is rather qualitative rather than quantative/definitive.



