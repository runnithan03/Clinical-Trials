install.packages("pwr")
library(pwr)

# Define the parameters
alpha <- 0.05
power_original <- 0.80
tau_M <- 1
k_values <- seq(0.1, 0.9, by = 0.1)
sigma_squared <- 1

effect_size <- k * tau_M

# sample size needed for original power
n_original <- pwr.t.test(d = tau_M, sig.level = alpha, power = power_original, type = "two.sample")$n

power_true_effect <- pwr.t.test(d = effect_size, n = n_original, sig.level = alpha, type = "two.sample")$power
print(power_true_effect) # power for k = 0.5

# Initialise a vector to store the power values
power_values <- numeric(length(k_values))

# Loop over the range of k values and calculate the power for each
for (i in seq_along(k_values)) {
  k <- k_values[i]
  effect_size <- k * tau_M
  power_values[i] <- pwr.t.test(d = effect_size, n = n_original, sig.level = alpha, type = "two.sample")$power
}

results <- data.frame(k = k_values, Power = power_values)
print(results)

plot(results$k, results$Power, type = "b", xlab = "k", ylab = "Power", main = "Power vs. k")
