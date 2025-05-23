# Simulating Ideal Sample Size
install.packages("ggplot2")

# Simulation to find the appropriate size for a trial with a binary primary outcome variable
binom_sim = function(
    nsim,    # Number of simulations
    N,       # Number of participants in each arm
    piC,     # Assumed proportion in control group
    piT      # minimum detectably different proportion in treatment group
){
  H0_reject_vec = rep(NA, nsim) # to store 1 if H0 rejected, 0 if fail to reject H0
  for (i in 1:nsim){
    rC = rbinom(n=1, size=N, prob=piC) # generate numbers of successes
    rT = rbinom(n=1, size=N, prob=piT)
    
    pC = rC/N # Find estimates of proportions
    pT = rT/N
    p_pooled = (rC+rT)/(2*N)
    
    z_stat = (pT - pC) / (sqrt(p_pooled*(1-p_pooled)*(1/N + 1/N)))
    
    p_val = 2*(1-pnorm(z_stat, mean=0, sd=1))
    
    H0_reject_vec[i] = ifelse(p_val<0.05, 1, 0)
  }
  
  power.est = mean(H0_reject_vec)
  power.est
}

exp_sim_LR = function(
    nsim,     # Number of simulations
    N,        # Number of participants per group
    rate_C,   # Hazard rate in control group (lambda_C)
    rate_T    # Hazard rate in treatment group (lambda_T)
){
  H0_reject_vec = rep(NA, nsim) # Store 1 if H0 rejected, 0 otherwise
  
  for (i in 1:nsim){
    # Simulate survival times from exponential distributions
    surv_C = rexp(n = N, rate = rate_C)  # Control group
    surv_T = rexp(n = N, rate = rate_T)  # Treatment group
    
    # Estimate means (just like proportions in binom_sim)
    mean_C = mean(surv_C)
    mean_T = mean(surv_T)
    
    # Pooled standard deviation estimate
    pooled_var = var(surv_C)/N + var(surv_T)/N
    z_stat = (mean_T - mean_C) / sqrt(pooled_var)
    
    # Two-sided p-value from standard normal
    p_val = 2 * (1 - pnorm(abs(z_stat)))
    
    H0_reject_vec[i] = ifelse(p_val < 0.05, 1, 0)
  }
  
  power.est = mean(H0_reject_vec)
  return(power.est)
}

# Define hazard rates
rate_C <- 103 / 274
rate_T <- 138 / 274

# Sample sizes to test
N_vals <- seq(240, 260, by = 1)
results <- data.frame(N = N_vals, power = NA)

# Loop through each N and estimate average power
for (j in 1:length(N_vals)) {
  sim_vec <- rep(NA, 100)
  for (i in 1:100) {
    sim_vec[i] <- exp_sim_LR(nsim = 100, N = N_vals[j],
                          rate_C = rate_C, rate_T = rate_T)
  }
  results$power[j] <- mean(sim_vec)
}

# View results
print(results)

sim_vec = rep(NA, 100)
for (i in 1:100){
  sim_vec[i] = exp_sim_LR(nsim=100, N=250, rate_C=rate_C, rate_T=rate_T)
}
library(ggplot2)
ggplot(mapping = aes(sim_vec)) + geom_histogram(bins=10)

participants_to_allocate <- read.csv("C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\CTs\\Summatives\\Summative 2\\participant_data.csv")
any(is.na(participants_to_allocate))
# all(is.na(participants_to_allocate))

# Allocation 
install.packages("Minirand")
library(Minirand)

colnames(participants_to_allocate)
table(participants_to_allocate$DiseaseLevel)
summary(participants_to_allocate$Age)

participants_to_allocate$AgeGroup <- ifelse(
  participants_to_allocate$Age < 70,
  "50-70",
  "70+"
)

nsample = nrow(participants_to_allocate)
covmat=as.matrix(participants_to_allocate[, c("DiseaseLevel", "AgeGroup")])
covwt=rep(1,1)/2
ratio=c(1,1)
ntrt=2
trtseq=c(0,1)
method="Range"
p = 0.5

res = rep(NA, nsample)
res[1] = sample(c(0,1), 1, replace = TRUE, prob = c(0.5,0.5)) 
# work through the remaining patients sequentially
for (j in 2:nsample){
  # get treatment assignment sequentially for all subjects
  # The vector res is updated and so all previous allocations are accounted for
  # covmat is the data frame of participant data - including only the covariates
  res[j] <- Minirand(
    covmat=covmat, j, covwt=covwt, ratio=ratio, ntrt=ntrt, trtseq=trtseq, method=method, result=res, p=p
  )
}

allocated_data <- read.csv("C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\CTs\\Summatives\\Summative 2\\participant_data.csv")
allocated_data$arm = res

## You can now investigate the balance of the design as usual
table(allocated_data$arm)

# covmat <- covmat %>% 
#   mutate(across(everything(), as.numeric))
# covmat

# Display the number of randomized subjects at covariate factors
balance1 <- randbalance(res, covmat=covmat, ntrt=2, trtseq=c(0,1)) 
balance1

# Calculate the total imbalance of the allocation
totimbal(trt = res, covmat = covmat, covwt = covwt, 
         ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = method)

allocated_data$arm <- ifelse(res == 0, "A", "B")
# View(allocated_data$arm)

table(allocated_data$arm)
colnames(allocated_data)

write.csv(allocated_data, 
          "C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\CTs\\Summatives\\Summative 2\\allocated_data.csv", 
          row.names = FALSE)

## Analysis
install.packages(c("HSAUR", "HSAUR3", "pROC", "survival", "ggsurvfit", "survminer", "car", "dplyr"))

library(HSAUR)
library(HSAUR3)
library(pROC)
library(survival)
library(ggsurvfit)
library(survminer)
library(car)
library(dplyr)

results <- read.csv("C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\CTs\\Summatives\\Summative 2\\trial_results.csv")
results

# D.2.1 Fitting a survival curve

# D.2.1.1 Kaplan-Meier
surv_obj = Surv(results$eventtime, results$status)
surv_obj

# You’ll notice that some have a ‘+’ attached. 
# This denotes the censored observations (the notation reflects the fact that the true time of death/the event will be greater than this).

# To fit a Kaplan-Meier survival curve, we use the function survfit, which is specified using a formula, much like lm or glm. 
# To fit a Kaplan-Meier estimate with a data frame split by treatment effect, the general form is:
km = survfit(surv_obj ~ DiseaseLevel, data=results)
summary(km)

km %>%  ggsurvfit() +
  add_censor_mark() +
  add_risktable() +
  add_confidence_interval()

# D.2.1.2 Fitting an exponential distribution

# To fit an exponential distribution, we need to estimate  λC and λT, using
mC = sum((results$status==1)&(results$DiseaseLevel=="Moderate"))
mT = sum((results$status==1)&(results$DiseaseLevel=="Severe"))
tsum_C = sum(results$eventtime[results$DiseaseLevel=="Moderate"])
tsum_T = sum(results$eventtime[results$DiseaseLevel=="Severe"])
lamhat_C = mC / tsum_C
lamhat_T = mT / tsum_T

# We can then plot the Survival curves using geom_function

# Define survival function for exponential density

exp_st = function(t, lambda){exp(-lambda*t)}

km %>% ggsurvfit() + ylim(0,1) + theme_bw() +
  add_censor_mark() +
  geom_function(fun=exp_st, args=list(lambda = lamhat_C), col="darkturquoise") +
  geom_function(fun=exp_st, args=list(lambda = lamhat_T), col="red") 

# D.2.2 Comparing survival curves

# Having found the MLEs for our dataset, assuming an exponential distribution, 
# we can now immediately conduct a likelihood ratio test.

# D.2.2.1 Likelihood ratio test

# We already have the mX and tX+ from Exercise D.9, and we can easily find t+ and m from these.
m = mT + mC
tsum = tsum_T + tsum_C

# and we can use these to compute λLR:
LRstat =  2*(mC*log(mC/tsum_C) + mT*log(mT/tsum_T) - m*log(m/tsum))
LRstat

# Finally, we refer this to Chi-Squared_1:
1-pchisq(LRstat, df=1)

# We find that we have easily enough evidence to reject H0 at the 95% level.

# D.2.2.2 Log-rank test
# The log-rank test is most easily found using the function survdiff. 
# This function has an argument rho that controls the type of test. 
# If we set rho=0 then it performs a log-rank test.

# The general form is similar to survfit:
survdiff(surv_obj ~ DiseaseLevel, data=results, rho=0)

# This is actually quite close to the result of the likelihood ratio test,
# in spite of the less-than-perfect fit of the exponential survival curve.

# D.2.2.3 Cox regression

# Finally, we will fit the Cox regression model. 
# This is done using the function coxph in the package survival.
initial_cox = coxph(formula = Surv(eventtime, status)~DiseaseLevel, data=results)
initial_cox

# There are baseline covariates here so use them in the Cox Regression model
cox = coxph(formula = Surv(eventtime, status)~DiseaseLevel + Age + arm, data=results)
cox
summary(cox)

# We can visualise this by further subsetting the Kaplan-Meier estimator
km_full = survfit(surv_obj ~ DiseaseLevel + arm, data=results)
km_full %>% ggsurvfit() +
  add_censor_mark() 

# eventtime is a numeric output, so one way to get a visual impression of its effect on the survival curve we can bin it. 
# For example, we can choose  eventtime ≤ 142 and eventtime > 142.
results$eventtime142 = sapply(1:nrow(results), function(i){ifelse(results$eventtime[i]>142, 1, 0)})
km_full = survfit(surv_obj ~ DiseaseLevel + eventtime142, data=results)
km_full %>% ggsurvfit()

# To check our model more carefully we can use the log-log plot for the categorical variables:
km_full = survfit(surv_obj ~ DiseaseLevel + arm, data=results)
ggsurvplot(km_full, fun = "cloglog")

# and the Schoenfeld residuals for the continuous residuals (after removing the insignificant variables):
cox_schon = coxph(formula = surv_obj~DiseaseLevel + Age +  arm, data=results)
ggcoxzph(cox.zph(cox_schon), var = "Age")

# and we see that there is insufficient evidence to reject the null hypothesis (of proportional hazards).

#Extra Simulation Method
exp_sim_cox = function(
    nsim,     # Number of simulations
    N,        # Number of participants per group
    rate_C,   # Hazard rate in control group (lambda_C)
    rate_T    # Hazard rate in treatment group (lambda_T)
){
  H0_reject_vec = rep(NA, nsim) # Store 1 if H0 rejected, 0 otherwise
  
  for (i in 1:nsim){
    # Simulate survival times from exponential distributions
    surv_C = rexp(n = N, rate = rate_C)  # Control group
    surv_T = rexp(n = N, rate = rate_T)  # Treatment group
    
    # Create a data frame for Cox regression
    # All events are observed (status = 1), as in the original function
    data = data.frame(
      time = c(surv_C, surv_T),
      status = rep(1, 2 * N),  # All events observed
      arm = rep(c("C", "T"), each = N)  # Control (C) and Treatment (T)
    )
    
    # Surv(data$time, data$status)
    
    # Fit Cox model
    cox_model = coxph(Surv(time, status) ~ arm, data = data)
    
    # Extract p-value for the treatment effect (arm)
    p_val = summary(cox_model)$coefficients["armT", "Pr(>|z|)"]
    
    # Reject H0 if p-value < 0.05
    H0_reject_vec[i] = ifelse(p_val < 0.05, 1, 0)
  }
  
  power.est = mean(H0_reject_vec)
  return(power.est)
}

# Parameters
nsim = 1000
N = 250  # Number of participants per group (example)
rate_C = log(2) / 103  # Hazard rate for control group (median = 103 days)
rate_T = log(2) / 138  # Hazard rate for treatment group (median = 138 days)

# Run the simulation
power = exp_sim_cox(nsim, N, rate_C, rate_T)
cat("Power:", power, "\n")

sim_vec = rep(NA, 100)
for (i in 1:100){
  sim_vec[i] = exp_sim_cox(nsim=100, N=N, rate_C=rate_C, rate_T=rate_T)
}
library(ggplot2)
ggplot(mapping = aes(sim_vec)) + geom_histogram(bins=10)

#Extra Simulation Method with Censoring
exp_sim_cox = function(
    nsim,     # Number of simulations
    N,        # Number of participants per group
    rate_C,   # Hazard rate in control group (lambda_C)
    rate_T    # Hazard rate in treatment group (lambda_T)
){
  H0_reject_vec = rep(NA, nsim) # Store 1 if H0 rejected, 0 otherwise
  
  for (i in 1:nsim){
    # Simulate survival times from exponential distributions
    surv_C = rexp(n = N, rate = rate_C)  # Control group
    surv_T = rexp(n = N, rate = rate_T)  # Treatment group
    
    # Combine survival times
    time = c(surv_C, surv_T)
    arm = rep(c("C", "T"), each = N)
    
    # Create initial status vector (1 for events)
    status = rep(1, 2 * N)
    
    # 5% random censoring for dropouts
    dropout_indices = sample(1:(2 * N), size = round(0.05 * 2 * N))
    status[dropout_indices] = 0
    time[dropout_indices] = runif(length(dropout_indices), 0, time[dropout_indices])
    
    # Censor everything after 274 days
    censored_after_274 = time > 274
    time[censored_after_274] = 274
    status[censored_after_274] = 0
    
    # Create data frame for Cox regression
    data = data.frame(
      time = time,
      status = status,
      arm = arm
    )
    
    # Fit Cox model
    cox_model = coxph(Surv(time, status) ~ arm, data = data)
    
    # Extract p-value for the treatment effect (arm)
    p_val = summary(cox_model)$coefficients["armT", "Pr(>|z|)"]
    
    # Reject H0 if p-value < 0.05
    H0_reject_vec[i] = ifelse(p_val < 0.05, 1, 0)
  }
  
  power.est = mean(H0_reject_vec)
  return(power.est)
}

# Parameters
nsim = 1000
N = 320  # Number of participants per group (example)
rate_C = log(2) / 103  # Hazard rate for control group (median = 103 days)
rate_T = log(2) / 138  # Hazard rate for treatment group (median = 138 days)

# Run the simulation
power = exp_sim_cox(nsim, N, rate_C, rate_T)
cat("Power:", power, "\n")

sim_vec = rep(NA, 100)
for (i in 1:100) {
  sim_vec[i] = exp_sim_cox(nsim = 10000, N = 320, rate_C = rate_C, rate_T = rate_T)
  # Print progress every 10%
  if (i %% 10 == 0) {
    cat("+", i, "% complete\n")
  }
}
library(ggplot2)
ggplot(mapping = aes(sim_vec)) + geom_histogram(bins=10)

table(sim_vec)
mean(sim_vec) # 0.89997

sim_vec_1 = rep(NA, 100)
for (i in 1:100){
  sim_vec_1[i] = exp_sim_cox(nsim=100, N=322, rate_C=rate_C, rate_T=rate_T)
}
library(ggplot2)
ggplot(mapping = aes(sim_vec_1)) + geom_histogram(bins=10)

table(sim_vec_1)
mean(sim_vec_1) # 0.90268

# Run this tomorrow...
sim_vec_2 = rep(NA, 100)
for (i in 1:100){
  sim_vec_2[i] = exp_sim_cox(nsim=10000, N=322, rate_C=rate_C, rate_T=rate_T)
  # Print progress every 10%
  if (i %% 10 == 0) {
    cat("+", i, "% complete\n")
  }
}
library(ggplot2)
ggplot(mapping = aes(sim_vec_2)) + geom_histogram(bins=10)

table(sim_vec_2)
mean(sim_vec_2)
