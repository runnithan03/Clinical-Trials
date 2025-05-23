# If you have time, you could now combine the work we’ve done in this practical and the previous one to implement power simulation for a binary outcome trial.
# If you’ve found the above fairly simple, this will be a good exercise to stretch your muscles. If not, then please still give it some thought, but feel free to look at the solution to see how the method applied to another approach.

# Before starting, look back at Practical 3 and the steps to the simulation method and consider what you’ve to do to complete the exercise.
# Exercise D.14 Use simulation to find the appropriate size for a trial with a binary primary outcome variable. You can assume that the proportion of ‘successes’ in the control group is 0.65, 
# and that the least clinically significant proportion we’d like to be able to detect is a treatment group proportion of 0.8. The trial should have a significance level of  α=0.05 and a power of  
# 1−β=0.8. 

# The steps within each simulation are:
# Simulate data using rbinom
# Run some analysis method on that simulated data
# Record whether H0 was rejected
# In the code below I’ve used the simplest type of analysis, a  t-test on the absolute risk difference, but you could use any of the methods we’ve covered.

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
# Let’s now try this for  N=100, as a base point
sim_vec7 = rep(NA, 100)
for (i in 1:100){
  sim_vec7[i] = binom_sim(nsim=100, N=100, piC=0.65, piT=0.8)
}
ggplot(mapping = aes(sim_vec7)) + geom_histogram(bins=10)

# This is clearly too low so we can try a higher number - let’s say  N=150
sim_vec8 = rep(NA, 100)
for (i in 1:100){
  sim_vec8[i] = binom_sim(nsim=100, N=150, piC=0.65, piT=0.8)
}
ggplot(mapping = aes(sim_vec8)) + geom_histogram(bins=10)

#This appears to be much closer - one could now finesse  N to fit some criteria on the proportion of simulations achieving the desired power.
# There are of course lots of ways we could refine and extend this, many of them similar to what we discussed in Section C, for example:
  
# Running more simulations (100 is really not enough, but was chosen so that the file would compile - feel free to increase it)
# Incorporating a better sampling scheme than SRS
# Improving the drop-out / censoring modelling
# Allowing the probability to vary in the data generation (this would be known as a sensitivity analysis).