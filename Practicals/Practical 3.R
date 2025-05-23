## C.1 Power simulation for a t-test

# The very simplest way would be the following, where we generate  
# N participants per arm, and their data.

ttest_sim = function(
    nsim,  # the number of simulations to use
    N,     # number of participants per trial arm
    mu,    # mean outcome (for group C / under H0)
    tm,    # minimum detectable effect size
    sd     # sd of outcome
){
  H0_reject_vec = rep(NA, nsim) # to store 1 if H0 rejected, 0 if fail to reject H0
  for (i in 1:nsim){
    out_groupC = rnorm(N, mean=mu, sd=sd)
    out_groupT = rnorm(N, mean=mu + tm, sd=sd)
    ttest = t.test(
      out_groupC, out_groupT, 
      alternative = "two.sided", 
      var.equal = T, conf.level = 0.95)
    H0_reject_vec[i] = ifelse(ttest$p.value<0.05, 1, 0)
  }
  
  power.est = mean(H0_reject_vec)
  power.est
}

# For a sanity check, we can put in the number we got from the equation in Section 2. 
# To generate the data we need to specify a mean

# relies on seeds, so the result changes
ttest_sim(nsim=100, N=113, mu=5, tm=3, sd=8)

# this is close to what we expect from the formula approach. 
# To see how consistent this result is, we can see what happens over many sets of simulations
sim_vec1 = rep(NA, 100)
for (i in 1:100){
  sim_vec1[i] = ttest_sim(nsim=100, N=113, mu=5, tm=3, sd=8)
}
require(ggplot2)
ggplot(mapping = aes(sim_vec1)) + geom_histogram(bins=10)

# To achieve 1−β =  0.9, we need to increase  N. 
# It turns out that around 155 participants in each arm gives a mean power of around 0.9.
sim_vec2 = rep(NA, 100)
for (i in 1:100){
  sim_vec2[i] = ttest_sim(nsim=100, N=155, mu=5, tm=3, sd=8)
}
ggplot(mapping = aes(sim_vec2)) + geom_histogram(bins=10)

# It is sensible to choose a number such that we are sufficiently confident of our study having the power we require. 
# We may for example choose N so that the estimated power is at least 0.9 in some high proportion (say 90%) of the simulations.

# The change we need to make is to include a variance parameter for each group, rather than a common one. 
# Note that we aren’t changing the t.test function to use unequal variances, as it is common practice to keep assuming equal variances, and actually these variances still aren’t that different.
ttest_sim2 = function(
    nsim,  # the number of simulations to use
    N,     # number of participants per trial arm
    mu,    # mean outcome (for group C / under H0)
    tm,    # minimum detectable effect size
    sdC,   # sd of outcome in group C
    sdT    # sd of outcome in group T
){
  H0_reject_vec = rep(NA, nsim) # to store 1 if H0 rejected, 0 if fail to reject H0
  for (i in 1:nsim){
    out_groupC = rnorm(N, mean=mu, sd=sdC)
    out_groupT = rnorm(N, mean=mu + tm, sd=sdT)
    ttest = t.test(
      out_groupC, out_groupT, 
      alternative = "two.sided", 
      var.equal = T, conf.level = 0.95)
    H0_reject_vec[i] = ifelse(ttest$p.value<0.05, 1, 0)
  }
  
  power.est = mean(H0_reject_vec)
  power.est
}

# If the variances are equal, this gives the same as before (approximately)
ttest_sim2(nsim=100, N=155, mu=5, tm=3, sdC=8, sdT=8)

# However, if we increase the outcome variance in group T by 25%, things change:
ttest_sim2(nsim=100, N=155, mu=5, tm=3, sdC=8, sdT=10)

sim_vec3 = rep(NA, 100)
for (i in 1:100){
  sim_vec3[i] = ttest_sim2(nsim=100, N=155, mu=5, tm=3, sdC=8, sdT=10)
}
ggplot(mapping = aes(sim_vec3)) + geom_histogram(bins=10)

# Figure C.3: Power estimates for 100 sets of 100 simulations with larger outcome SD in group T.
# By contrast, if the outcome variance is lower (this is much less likely) then the power will be higher than calculated:
sim_vec4 = rep(NA, 100)
for (i in 1:100){
  sim_vec4[i] = ttest_sim2(nsim=100, N=155, mu=5, tm=3, sdC=8, sdT=6)
}
ggplot(mapping = aes(sim_vec4)) + geom_histogram(bins=10)

# Recreate the function from Exercise C.2, where the two groups have different variances. 
# This time, simulate the data sequentially, using simple random sampling as the allocation method.
# What happens to the power distribution for the scenario we considered in that question?
# What happens for a much smaller N, say N=25?

one_trial_srs = function(
    N,     # number of participants per trial arm
    mu,    # mean outcome (for group C / under H0)
    tm,    # minimum detectable effect size
    sdC,   # sd of outcome in group C
    sdT    # sd of outcome in group T
){
  ## Create empty vectors to contain output for groups C and T
  outC = integer()
  outT = integer()
  
  for (i in 1:(2*N)){
    # allocate using SRS
    arm = sample(c("C", "T"), size=1)
    # generate outcome according to groups' distributions
    if (arm == "C"){
      outi = rnorm(1, mean=mu, sd=sdC)
      outC = c(outC, outi)
    } else if (arm == "T"){
      outi = rnorm(1, mean=mu+tm, sd=sdT)
      outT = c(outT, outi)
    }
  }
  # conduct t-test for this trial
  t.test(x=outC, y=outT, 
         alternative = "two.sided", paired=F, 
         var.equal=T, conf.level=0.95)
}

ttest_sim_srs = function(
    nsim,  # the number of simulations to use
    N,     # number of participants per trial arm
    mu,    # mean outcome (for group C / under H0)
    tm,    # minimum detectable effect size
    sdC,   # sd of outcome in group C
    sdT    # sd of outcome in group T
){
  H0_reject_vec = rep(NA, nsim) # to store 1 if H0 rejected, 0 if fail to reject H0
  
  for (i in 1:nsim){
    trial_i = one_trial_srs(N, mu, tm, sdC, sdT)
    H0_reject_vec[i] = ifelse(trial_i$p.value<0.05, 1, 0)
  }
  
  power.est = mean(H0_reject_vec)
  power.est
}

ttest_sim_srs(nsim=100, N=155, mu=5, tm=3, sdC=8, sdT=10)

sim_vec5 = rep(NA, 100)
for (i in 1:100){
  sim_vec5[i] = ttest_sim_srs(nsim=100, N=155, mu=5, tm=3, sdC=8, sdT=10)
}
ggplot(mapping = aes(sim_vec5)) + geom_histogram(bins=10)

# Figure C.5: Power estimates for 100 sets of 100 simulations with N=155 with larger outcome SD in group T, with SRS allocation.
# and we see that the spread is much higher than when we had two exactly equal groups, as in Figure C.3.
# If N is smaller, the variance in the power estimates will be higher, because the potential for imbalance will be greater.

## C.2 Power simulation for ANCOVA
# only need ANCOVA for power simulation rather than the t-test
# In our analysis section we learned that the ANCOVA model is much more efficient and flexible than a t-test, so now we’ll simulate our data with that analysis method in mind.

## Function to simulate one trial (with ANCOVA analysis)
ancova_trial_srs = function(
    N,       # Number of participants per group
    mu_B,    # baseline mean
    mu,      # outcome mean (control group / H_0)
    rho,     # correlation between baseline and outcome
    tm,      # minimum detectable effect size
    sd_eps,  # SD of error
    sd_B     # SD of baseline measurement
){
  ## Empty data frame for trial data
  trial_mat = matrix(NA, ncol=3, nrow=2*N)
  trial_df = data.frame(trial_mat)
  names(trial_df) = c("baseline", "arm", "outcome")
  
  for (i in 1:(2*N)){
    bas_i = rnorm(1, mean=mu_B, sd=sd_B)
    trial_df$baseline[i] = bas_i
    alloc_i = sample(c("C", "T"), 1) # Using SRS in this function
    trial_df$arm[i] = alloc_i
    eps_i = rnorm(1, mean=0, sd=sd_eps)
    if(alloc_i == "C"){
      out_i = mu + rho*(bas_i - mu_B) + eps_i
    } else if (alloc_i == "T"){
      out_i = mu + tm + rho*(bas_i - mu_B) + eps_i
    }
    trial_df$outcome[i] = out_i
  }
  model.fit = lm(outcome ~ baseline + arm, data=trial_df)
  summary(model.fit)
}

## Function to simulate many trials (with ANCOVA analysis)

ancova_sim_srs = function(
    nsim,  # the number of simulations to use
    N,       # Number of participants per group
    mu_B,    # baseline mean
    mu,      # outcome mean (control group / H_0)
    rho,     # correlation between baseline and outcome
    tm,      # minimum detectable effect size
    sd_eps,  # SD of error
    sd_B     # SD of baseline measurement
){
  H0_reject_vec = rep(NA, nsim) # to store 1 if H0 rejected, 0 if fail to reject H0
  
  for (i in 1:nsim){
    trial_i = ancova_trial_srs(N, mu_B, mu, rho, tm, sd_eps, sd_B)
    H0_reject_vec[i] = ifelse(trial_i$coefficients[3,4]<0.05, 1, 0)
  }
  
  power.est = mean(H0_reject_vec)
  power.est
}

sim_vec6 = rep(NA, 100)
for (i in 1:100){
  sim_vec6[i] = ancova_sim_srs(nsim=100, N=90, mu_B=50, mu=60, 
                               rho=0.65, tm=3, sd_eps=6, sd_B = 8)
}
ggplot(mapping = aes(sim_vec6)) + geom_histogram(bins=10)

# Figure C.6: Power estimates for 100 sets of 100 simulations with larger outcome SD, with SRS allocation, with  
# N=75, assuming an ANCOVA analysis. The power is much higher than for our t-test simulation with N=100. 
# This is because the variance of the treatment effect estimate is lower for ANCOVA for the t-test 
# (by a factor of  (1-p^2), as we saw in Section 4.3.1.2.

# We could extend this in many ways to better understand the uncertainty around the power for a given sample size. For example:
# Trying different allocation methods
# Incorporating uncertainty around quantities like  μ,μB,ρ etc.
# If you still have time you might like to look into these as an extension problem.



