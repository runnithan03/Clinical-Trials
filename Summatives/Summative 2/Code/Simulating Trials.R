install.packages("ggplot2")

exp_sim = function(
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
    
    # Calculate number of events and total time at risk
    mC = N  # Number of events in control group (all observed)
    mT = N  # Number of events in treatment group (all observed)
    tSUM_C = sum(surv_C)  # Total time at risk in control group
    tSUM_T = sum(surv_T)  # Total time at risk in treatment group
    
    # Total events and total time
    m = mT + mC
    tSUM = tSUM_T + tSUM_C
    
    # Compute the LRT statistic
    # Avoid log(0) or division by 0 by adding a small constant if necessary
    if (mC == 0 || tSUM_C == 0 || mT == 0 || tSUM_T == 0 || m == 0 || tSUM == 0) {
      LRstat = 0  # If any term is 0, set LRstat to 0 (no test possible)
    } else {
      LRstat = 2 * (mC * log(mC / tSUM_C) + mT * log(mT / tSUM_T) - m * log(m / tSUM))
    }
    
    # Compute p-value from chi-square distribution with 1 df
    p_val = pchisq(LRstat, df = 1, lower.tail = FALSE)
    
    # Reject H0 if p-value < 0.05
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
    sim_vec[i] <- exp_sim(nsim = 100, N = N_vals[j],
                          rate_C = rate_C, rate_T = rate_T)
  }
  results$power[j] <- mean(sim_vec)
}

# View results
print(results)

sim_vec = rep(NA, 100)
for (i in 1:100){
  sim_vec[i] = exp_sim(nsim=100, N=250, rate_C=rate_C, rate_T=rate_T)
}
library(ggplot2)
ggplot(mapping = aes(sim_vec)) + geom_histogram(bins=10)

participants_to_allocate <- read.csv("C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\CTs\\Summatives\\Summative 2\\participant_data.csv")

# Adjusting for the notes about sample size calculations for time-to-event data

# Load required library for Cox model
library(survival)

# Simulation function
exp_sim_cox = function(
    nsim,     # Number of simulations
    N,        # Number of participants per group
    rate_C,   # Hazard rate in control group (lambda_C)
    rate_T,   # Hazard rate in treatment group (lambda_T)
    followup_time = 274,  # Follow-up period (9 months = 274 days)
    loss_to_followup = 0.05  # Proportion lost to follow-up
){
  H0_reject_vec = rep(NA, nsim) # Store 1 if H0 rejected, 0 otherwise
  
  for (i in 1:nsim){
    # Simulate survival times from exponential distributions
    surv_C = rexp(n = N, rate = rate_C)  # Control group
    surv_T = rexp(n = N, rate = rate_T)  # Treatment group
    
    # Combine survival times and arm
    eventtime = c(surv_C, surv_T)
    arm = rep(c("C", "T"), each = N)
    
    # Simulate age and bin into two groups: below 70 and 70+
    age = runif(2 * N, min = 50, max = 90)  # Age between 50 and 90
    age_group = ifelse(age < 70, 0, 1)  # 0 for below 70, 1 for 70 and above
    
    # Initial status (assuming your existing logic sets this)
    # If you've already handled censoring for the 274-day follow-up, this will be updated
    status = rep(1, 2 * N)  # Placeholder; assuming your results$status will override this
    
    # Simulate loss to follow-up (5% of participants)
    lost = rbinom(2 * N, size = 1, prob = loss_to_followup)
    lost_time = runif(2 * N, min = 0, max = followup_time)  # Random time for loss
    for (j in 1:(2 * N)) {
      if (lost[j] == 1) {
        # If lost to follow-up, censor at the random time
        if (lost_time[j] < eventtime[j]) {
          eventtime[j] = lost_time[j]
          status[j] = 0
        }
      }
    }
    
    # Create the results data frame
    results = data.frame(
      eventtime = eventtime,
      status = status,
      arm = arm,
      age_group = age_group
    )
    
    # Create the survival object
    surv_obj = Surv(results$eventtime, results$status)
    
    # Fit Cox model using the survival object
    cox_model = coxph(surv_obj ~ arm + age_group, data = results)
    
    # Extract p-value for the treatment effect (arm)
    p_val = summary(cox_model)$coefficients["armT", "Pr(>|z|)"]
    
    # Reject H0 if p-value < 0.05
    H0_reject_vec[i] = ifelse(p_val < 0.05, 1, 0)
  }
  
  power.est = mean(H0_reject_vec, na.rm = TRUE)
  return(power.est)
}

# Set parameters based on the scenario
rate_C = log(2) / 103  # Hazard rate for control (median = 103 days)
rate_T = log(2) / 138  # Hazard rate for treatment (median = 138 days)

# Run the simulation
set.seed(123)  # For reproducibility
power = exp_sim_cox(
  nsim = 1000,  # Number of simulations
  N = 250,      # Number of participants per group
  rate_C = rate_C,
  rate_T = rate_T
)

# Print result
cat("Estimated power:", power, "\n")
