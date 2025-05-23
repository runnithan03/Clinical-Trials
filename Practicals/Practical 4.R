install.packages(c("HSAUR", "HSAUR3", "pROC", "survival", "ggsurvfit", "survminer", "car","ggplot2", "dpylr"))

library(HSAUR)
data("respiratory")
# Keep only months 0 and 4
resp_04 = respiratory[respiratory$month %in% c(0,4),]
# Use month 4 data to form basis of overall dataset
resp_df = respiratory[respiratory$month %in% c(4),]
# add month 0 statusses in as baseline, and remove month column
resp_df$status0 = resp_04$status[resp_04$month==0]
resp_df = resp_df[ ,-6]

# For each patient, we have the following baseline covariates: sex, age, 
# treatment centre (centre 1 or centre 2) and symptom status (poor or good).

# The outcome variable is whether the status of the patient’s symptoms are poor or good after four months of the trial. 
# Note that unusually we do have the baseline observation of this quantity. We will build up through some of the models we studied in Chapter 6.

# Our first step will be to fit a 95% Wilson confidence interval for the absolute risk difference, as we did in Section 6.2.1.3.
rC = sum((resp_df$status == "good")&(resp_df$treatment == "placebo"))
nC = sum(resp_df$treatment == "placebo")
pC = rC / nC
rT = sum((resp_df$status == "good")&(resp_df$treatment == "treatment"))
nT = sum(resp_df$treatment == "treatment")
pT = rT / nT

# To find the limits of the individual CIs, we need to rearrange the quadratic equation above into the form we’re used to.
# It is simple to define a function to solve quadratic equations in this form, for example
quad_fun = function(a,b,c){
  disc = b^2 - 4*a*c
  lower = (-b - sqrt(disc))/(2*a)
  upper = (-b + sqrt(disc))/(2*a)
  c(lower, upper)
}
# and we can now find the limits of each CI:
# Control group
lim_resp_C = quad_fun(
  a = 1 + qnorm(0.975)^2/nC,
  b = -2*pC - qnorm(0.975)^2/nC,
  c = pC^2
)
lim_resp_C

# Treatment group
lim_resp_T = quad_fun(
  a = 1 + qnorm(0.975)^2/nT,
  b = -2*pT - qnorm(0.975)^2/nT,
  c = pT^2
)
lim_resp_T

# Notice that because pC is quite close to 0.5, the Newcome interval for  
# πC is close to symmetric. The CI for πT is more skewed.
# Finally we can use these values to find the Wilson CI for the Absolute Risk Difference (ARD):
# Make objects match notation
lC = lim_resp_C[1]
uC = lim_resp_C[2]
lT = lim_resp_T[1]
uT = lim_resp_T[2]

l_newc = pT - pC - sqrt((pT - lT)^2 + (uC - pC)^2)
u_newc = pT - pC + sqrt((uT - pT)^2 + (pC - lC)^2)

c(l_newc, u_newc)

# Finally we can use this to find a 95% CI for the NNT by finding the reciprocals of the limits of the CI for ARD:
c(1/u_newc, 1/l_newc)

# This is a huge range and would make it difficult to make a clinical decision.
# With the quantities we’ve found in exercise D.1 we can also find an approximate 95% CI for the risk ratio (see Section 6.2.2.1).

# Find a 95% CI for the risk ratio  πT/πC and comment on your result.

# D.1.2 Logistic regression
# We next looked at logistic regression as way to account for the baseline covariates. 
# To fit a logistic regression in R we use the function glm. This fits all kinds of generalised linear models, 
# and so we specify that we want a logistic regression by choosing family = binomial(link='logit').
# The code is included in all examples in Sections 6.3 and 6.4, so if you aren’t sure how to do something they would be a good place to look.

# Exercise D.3 Fit a logistic regression model to the resp_df data set, being careful to check for any issues with the data.
model1 = glm(status ~ centre + treatment + gender + age + status0, 
             family = binomial(link='logit'), data=resp_df)
summary(model1)

# and we can inspect the variance inflation factor to check for multicollinearity:
library(car)
vif(model1)

# We can also try one with age squared, since the effect might not be linear (the age range is 11 to 68 which is quite a wide range). 
# While we’re at it we’ll remove sex from the model since it appears to have no effect
model2 = glm(status ~ centre + treatment + I(age^2) + status0, 
             family = binomial(link='logit'), data=resp_df)
summary(model2)
vif(model2)

# Because the age range actually goes quite low, we could try a model with a log transform
model3 = glm(status ~ centre + treatment + I(log(age)) + status0, 
             family = binomial(link='logit'), data=resp_df)
summary(model3)
vif(model3)

# Exercise D.4 Perform some diagnostics on the model you fitted in Exercise D.3, 
# focussing on the model’s ability to discriminate between outcomes, and also on how well the model calibrates to the observed data.
# Firstly we’ll use ROC analysis to assess the model in terms of how well it discriminates the participants in terms of their respiratory status at 4 months.
library(pROC)
fit_resp = fitted(model3)   # Fitted values from model3
out_resp = resp_df$status   # outcome values (1 or 2)
roc_resp_df = data.frame(fit = fit_resp, out = out_resp)
roc_resp = roc(data=roc_resp_df, response = out, predictor=fit)
roc_resp

# We can also plot the ROC curve
library(ggplot2)
ggroc(roc_resp, legacy.axes=T) + geom_abline(slope=1, intercept=0, type=2)

# Next we’ll assess the model in terms of calibration.
## Group the observations into age groups (I've chosen 10 year bins)
resp_df$age10 = round(resp_df$age, -1)
# find mean status (minus 1 because factor levels are 1 and 2) and number of obs
# for each combination of factor/age group levels
library(dplyr)
resp_sum = resp_df %>% 
  group_by(age10, centre, treatment, status0) %>%
  summarise(mean = mean(as.numeric(status))-1, n=length(age10), .groups="keep")

# Plot the observations, using facet_wrap to deal with two of the factors 

obs_plot = ggplot(data=resp_sum, aes(x=age10, col=treatment)) +
  geom_point(aes(y=mean, size=n), pch=16) + 
  theme_bw()+
  facet_wrap(c("centre", "status0")) + theme(legend.position = "bottom")

## To include the estimate and SD from the model, we use the original dataset with continuous age,
# fit model 3 including an estimate of SE, and use geom_line and geom_ribbon to add the fitted model with 95% intervals

fit_resp = predict(model3, newdata=resp_df, se.fit=T, type="response")

resp_df$fit = fit_resp$fit
resp_df$fit_se = fit_resp$se.fit

obs_plot + geom_line(data=resp_df, aes(x=age, y=fit)) +
  geom_ribbon(data=resp_df, aes(x=age, ymin = fit - 1.96*fit_se, ymax = fit + 1.96*fit_se, fill = treatment), alpha=0.3)

# Figure D.1: Observations and fitted model for combinations of centre (1: top and 2: bottom) for baseline status (poor: left and good: right)
# Overall, the model does not appear to be a great fit, but there is no systematic cause for concern.
# These plots also give us an idea of the model’s fit for different categories of patient. 
# For example, the probability of good respiratory symptoms at 4 months is much higher for a patient in the treatment group who had good symptoms at month zero, especially if they are young and belong to centre 2. 
# The model’s output is much more bleak for a patient in group C at centre 1 who had poor symptoms at month zero.

# Having performed some diagnostics, we can proceed to use our model to provide information about the effect of the treatment to improve respiratory symptoms.
# Exercise D.5: Compute a 95% CI for the odd ratio for this trial, using the model you built in Exercise D.3.
summary(model3)

# and find that the estimate of the log OR is 1.0257, with an SE of 0.4605, and therefore a 95% CI for the log OR is
est_logOR = summary(model3)$coefficients[3,1]
se_logOR = summary(model3)$coefficients[3,2]
logOR_CI = c(est_logOR - qnorm(0.975)*se_logOR, est_logOR + qnorm(0.975)*se_logOR)

# Finally we can use this to find a CI for the OR itself
exp(logOR_CI)

# This is further away from the null value (1 for the OR) than any of our previous confidence intervals, showing that including the baseline covariates has reduced our uncertainty in the treatment effect in a helpful way.

# D.2 Analysis for Survival data
library(survival)
library(ggsurvfit)
library(survminer)

# The first dataset we will use is aml - this dataset is from a trial investigating whether the standard course of chemotheraphy should be extended for some additional cycles (‘maintenance’) for patients with Acute Myelogenous Leukemia.
aml

# The first step is to combine the first two columns into a form we can use. We do this using the Surv function in the package survival, which creates a ‘survival object’ that we can then use in various other functions. 
# This object contains the times and the information about which observations are censored.
# Use the Surv function now on the aml data. The output will contain some notation you probably haven’t seen before - can you work out what it means?
surv_aml = Surv(aml$time, aml$status)

# This is a vector of the time values, in the same order as in aml. You’ll notice that some have a ‘+’ attached. This denotes the censored observations (the notation reflects the fact that the true time of death/the event will be greater than this).

#D.2.1 Fitting a survival curve
# The next thing we probably want to do is to estimate the survival function and plot the survival curve. The first method we’ll use is the Kaplan-Meier estimator.

# D.2.1.1 Kaplan-Meier
# To fit a Kaplan-Meier survival curve, we use the function survfit, which is specified using a formula, much like lm or glm. To fit a Kaplan-Meier estimate with a data frame split by treatment effect, the general form is
# survfit(<surv_obj> ~ <treatment var>, data = <dataframe>)
# We can then use summary to see the intermediary calculations at each step, and plot (for base plot) or ggsurvfit (from ggplot2 and ggsurvfit) to plot the curves.
# To fit the Kaplan-Meier estimator we use:
km_aml <- survfit(surv_aml ~ x, data = aml)

# We can then look at the summary table and plot the data by
summary(km_aml)

km_aml %>%  ggsurvfit() +
  add_censor_mark() +
  add_risktable() +
  add_confidence_interval()

# We can see from the table that the lower curve is the non-maintained arm - there is only one survivor of this group, and the data finish at t=45.
# The function add_risktable() adds information about the numbers at risk and numbers of events. The function add_confidence_interval() adds a confidence interval, and we see that with the aml data the uncertainty is huge.

# D.2.1.2 Fitting an exponential distribution
# Exercise D.9 Fit an exponential distribution for each treatment group to the aml data and plot the resulting estimated survival curves, along with the Kaplan Meier estimators from Exercise D.8 (for comparison).
mC_aml = sum((aml$status==1)&(aml$x=="Nonmaintained"))
mT_aml = sum((aml$status==1)&(aml$x=="Maintained"))
tsum_aml_C = sum(aml$time[aml$x=="Nonmaintained"])
tsum_aml_T = sum(aml$time[aml$x=="Maintained"])
lamhat_aml_C = mC_aml / tsum_aml_C
lamhat_aml_T = mT_aml / tsum_aml_T

# We can then plot the Survival curves using geom_function


# Define survival function for exponential density
exp_st = function(t, lambda){exp(-lambda*t)}
km_aml %>% ggsurvfit() + ylim(0,1) + theme_bw() +
  add_censor_mark() +
  geom_function(fun=exp_st, args=list(lambda = lamhat_aml_C), col="darkturquoise") +
  geom_function(fun=exp_st, args=list(lambda = lamhat_aml_T), col="red") 

# D.2.2 Comparing survival curves
# Having found the MLEs for the aml dataset, assuming an exponential distribution, we can now immediately conduct a likelihood ratio test.

# D.2.2.1 Likelihood ratio test
# Recall that our test statistic (which we found in Section 8.1) is ...
# which we then refer to a  Chi-Squared_1 distribution.
# Exercise D.10 Conduct a likelihood ratio test for the aml data, with the null hypothesis that the survival curves are the same for both treatment groups. Before you calculate the answer, think about what you expect to see.

# We already have the mX and t + X from Exercise D.9, and we can easily find t^+ and m from these.
m_aml = mT_aml + mC_aml 
tsum_aml = tsum_aml_T + tsum_aml_C #and we can use these to compute  λLR:
LRstat_aml =  2*(mC_aml*log(mC_aml/tsum_aml_C) + mT_aml*log(mT_aml/tsum_aml_T) - m_aml*log(m_aml/tsum_aml))
LRstat_aml
## [1] 4.061349
# Finally, we refer this to  Chi-Squared 1
1-pchisq(LRstat_aml, df=1)
## [1] 0.04387544
# We find that we have just enough evidence to reject  H0 at the 95% level.

# D.2.2.2 Log-rank test
# The log-rank test is most easily found using the function survdiff. This function has an argument rho that controls the type of test. If we set rho=0 then it performs a log-rank test.
# The general form is similar to survfit:
survdiff(surv_aml ~ x, data = aml, rho = 0)

# This is not that close to the result of the likelihood ratio test, probably reflecting the less-than-perfect fit of the exponential survival curve.

# D.2.2.3 Cox regression

# Exercise D.12 Fit a Cox proportional hazards model to the data. Do you think this is an appropriate model to use? How influential is the treatment?
# Solution. In the aml dataset there are no baseline covariates, so our only dependent variable is the treatment group variable x.
cox_aml = coxph(formula = Surv(time, status)~x, data=aml)
cox_aml
## Call:
## coxph(formula = Surv(time, status) ~ x, data = aml)
## 
##                  coef exp(coef) se(coef)     z      p
## xNonmaintained 0.9155    2.4981   0.5119 1.788 0.0737
## 
## Likelihood ratio test=3.38  on 1 df, p=0.06581
## n= 23, number of events= 18
# We’ve already seen that the survival curves don’t cross, so we can be reasonably comfortable fitting this model. To see how the model performs, we can use the summary function:
summary(cox_aml)
library(survival)
library(ggsurvfit)
library(survminer)

# To check this more carefully we can use the log-log plot
print(surv_aml)

# Log-Log plot - HOW TO PLOT THE LOG-LOG PLOT??!!
km_cox = survfit(Surv(time, status)~x, data=aml)
ggsurvplot(km_cox, fun = "cloglog")
# and we see that the curves are close to parallel.

# Although the coefficient of x isn’t quite significant at the 0.05 level, 
# it appears there is reasonable evidence that the treatment has an effect (though probably not enough to make a clinical decision).

# Exercise D.13 The dataset colon, also from the survival package, contains data from a trial of colon cancer patients, comparing three treatments: observation (obs), levamisole (Lev) and levamisole + 5-FU (Lev+5FU). To simplify things, we will restrict the data to those patients on Obs or Lev+5FU. The main report of this trial is Moertel et al. (1990).

colondf = colon[colon$rx!="Lev",]
colondf$rx = as.factor(as.character(colondf$rx)) # Removes Lev factor level

# For this data
# Look at the help file and make sure you understand what the columns mean
# Fit Kaplan-Meier estimators to the survival curves for the two groups.
# Perform some tests to see whether the treatment has a significant effect on the outcome. What do you find?

# From the help file we see that rx is the treatment group, stop is the time variable, status gives the censoring status.
str(colondf)

# We can therefore fit (and plot) the Kaplan-Meier estimator split by treatment group.
km_colon = survfit(Surv(time, status) ~ rx, data=colondf)
km_colon %>% ggsurvfit + ylim(0,1) +
  add_censor_mark() +
  add_confidence_interval()

# Next we can find the MLE for each treatment group
mC_colon = sum((colondf$status==1)&(colondf$rx=="Obs"))
mT_colon = sum((colondf$status==1)&(colondf$rx=="Lev+5FU"))
tsum_colon_C = sum(colondf$time[colondf$rx=="Obs"])
tsum_colon_T = sum(colondf$time[colondf$rx=="Lev+5FU"])
lamhat_colon_C = mC_colon / tsum_colon_C
lamhat_colon_T = mT_colon / tsum_colon_T

# We can then plot the Survival curves using geom_function
# Define survival function for exponential density
km_colon %>% ggsurvfit() +
  add_censor_mark() +
  add_risktable() +
  add_confidence_interval() + 
  ylim(0,1) + theme_bw() +
  geom_function(fun=exp_st, args=list(lambda = lamhat_colon_C), col="darkturquoise") +
  geom_function(fun=exp_st, args=list(lambda = lamhat_colon_T), col="red") 

# Figure D.2: Survival curves for groups C and T in colon study, fitted assuming an exponential distribution. Kaplan-Meier estimates also shown.
# We see that this fit is quite poor.
m_colon = mT_colon + mC_colon
tsum_colon = tsum_colon_T + tsum_colon_C

# and we can use these to compute  λLR:
LRstat_colon =  2*(mC_colon*log(mC_colon/tsum_colon_C) + mT_colon*log(mT_colon/tsum_colon_T) - m_colon*log(m_colon/tsum_colon))
LRstat_colon

# Finally, we refer this to Chi-Sqaured 1
1-pchisq(LRstat_colon, df=1)

# Highly significant, but because of the very poor fit in Figure D.2 not especially trustworthy.
# Next, we can perform a log-rank test:
survdiff(Surv(time, status)~rx, data=colondf, rho=0) # rx is the treatment group

# We see that the test statistic, although still very significant, is much lower than for the likelihood ratio test.
# Finally, we fit a Cox regression model. In the first instance we can do this with just the treatment group as a covariate:
coxph(formula = Surv(time, status)~rx, data=colondf)

# But we can also include other baseline covariates
coxph(formula = Surv(time, status)~rx + sex + age + obstruct + nodes, data=colondf)

# and we see that as well as the treatment arm, the number of lymph nodes with detectable cancer (given by nodes) is highly significant, and sex is also fairly significant.
# We can visualise this by further subsetting the Kaplan-Meier estimator
km_colon_bl = survfit(Surv(time, status) ~ rx + sex, data=colondf)
km_colon_bl %>% ggsurvfit() +
  add_censor_mark() 

# nodes is a numeric output, so one way to get a visual impression of its effect on the survival curve we can bin it. For example, we can choose  
# nodes≤4 and nodes > 4.
colondf$nodes4 = sapply(1:nrow(colondf), function(i){ifelse(colondf$nodes[i]>4, 1, 0)})
km_colon_bl = survfit(Surv(time, status) ~ rx + nodes4, data=colondf)
km_colon_bl %>% ggsurvfit()

# To check our model more carefully we can use the log-log plot for the categorical variables:
km_colon = survfit(Surv(time, status) ~ rx + sex, data=colondf)
ggsurvplot(km_colon, fun = "cloglog")

# and the Schoenfeld residuals for the continuous residuals (after removing the insignificant variables):
cox_colon = coxph(formula = Surv(time, status)~rx + sex + nodes, data=colondf)
ggcoxzph(cox.zph(cox_colon), var = "nodes")
# and we see that there is insufficient evidence to reject the null hypothesis (of proportional hazards).