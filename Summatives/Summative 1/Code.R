install.packages(c("naniar", "smdi", "tidyverse", "visdat", "rstanarm", "medicaldata","ggplot2","gridExtra","knitr", "parallel", "dplyr","Minirand","car","gtsummary"))


require(medicaldata)
library(car)
library(naniar)
library(smdi)
library(tidyverse)
library(visdat)
library(rstanarm)
library(ggplot2)
library(gridExtra)
library(knitr)
library(parallel)
library(dplyr)
library(Minirand)
library(gtsummary)


participants_df <- read.csv("C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\CTs\\Summatives\\Summative 1\\participant_data.csv", 
               stringsAsFactors = FALSE)
# sum(is.na(df))
colSums(is.na(participants_df)) # baseline has all the missing variables
colnames(participants_df)

qqplot(participants_df$age, participants_df$baseline)
qqplot(participants_df$BMI, participants_df$baseline)
hist(participants_df$age)
hist(participants_df$BMI)
hist(participants_df$baseline)

baseline <- participants_df$baseline
qqnorm(baseline, main = "Q-Q Plot of Baseline Variable", pch = 16, col = "blue")
qqline(baseline, col = "red", lwd = 2)  # Add reference line

age <- participants_df$age
qqnorm(age, main = "Q-Q Plot of Age Variable", pch = 16, col = "blue")
qqline(age, col = "red", lwd = 2)  # Add reference line

BMI <- participants_df$BMI
qqnorm(BMI, main = "Q-Q Plot of BMI Variable", pch = 16, col = "blue")
qqline(BMI, col = "red", lwd = 2)  # Add reference line

# Sex, Age, Baseline, BMI, Arm
ggplot(data = participants_df, aes(x=age, y=baseline, col = sex)) + 
  geom_point() +
  xlab("Age") + 
  ylab("Baseline") +
  ggtitle("Sex vs. Baseline") +
  theme_bw() # Significant

ggplot(data = participants_df, aes(x=age, y=baseline, col = BMI)) + 
  geom_point() +
  xlab("Age") + 
  ylab("Baseline") +
  ggtitle("BMI vs. Baseline") +
  theme_bw() # Not Significant

ggplot(data = participants_df, aes(x=BMI, y=baseline, col = age)) + 
  geom_point() +
  xlab("BMI") + 
  ylab("Baseline") +
  ggtitle("Age vs. Baseline") +
  theme_bw()

ggplot(df_category, aes(x = sex, y = baseline, fill = sex)) +
  geom_boxplot() +
  labs(title = "Baseline Symptom Distribution by Sex", x = "Sex", y = "Baseline Value") +
  theme_minimal()

## Imputing Missing Data
# vis_dat(participants_df, sort_type=F)
# gg_miss_var(participants_df)
# gg_miss_case(participants_df)
# md.pattern(participants_df)
# miss_case_summary(participants_df)

smdi_hotelling(participants_df) # p < 0.001 -> there is sufficient evidence to suggest MAR or MNAR.

asmd_participants = smdi_asmd(participants_df[ ,-1], includeNA=T) # asmd_min > 0.1
table1 <- kable(asmd_participants$baseline$asmd_table1)
table1

asmd_participants$baseline$asmd_plot
# although the ASMD values are quite large (much bigger than the advised 0.1), because there are a very small number of them they are not statistically significant.

# The only way we could determine whether the mechanisms for participants was MAR or MNAR would be to measure (or otherwise procure) some of the missing data. 
# We simply do not have the necessary information to work out which is the case. If this were a real trial, we would now talk at length with the experts/clinicians, 
# who will have a much better understanding of the probable causes of missingness.

# see Practical 1 for more ... but: mean imputation is not effective because ISETS MCAR
# imputing by logic won't necessarily work as well...

participants_nab <- nabular(participants_df)
ggplot(participants_nab,
       aes(x = BMI,
           fill = baseline_NA)) + 
  geom_histogram()

str(participants_nab)

# We’ve already looked at using logistic regression to understand patterns of missingness, and so it may not come as a surprise that we can use regression models to choose appropriate values for imputation. 
# This won’t be the same model, since we’re now interested in the value, rather than the missingness. Which type of regression model we use depends on the type of the variable we’re imputing values for. 
# If the variable is continuous, linear regression is likely to work well. If the variable is binary, we should try logistic regression. There are plenty of other types of model we could use 
# (as well as a whole host of machine learning type models!) but in this practical we’ll stick to those two.

# The model involves variables with no missingness
baseline_lm1 = stan_glm(
  baseline ~ sex + age + BMI,
  data = participants_df
)
summary(baseline_lm1)

# Sex (Male) is a strong predictor → Males have ~9.3 lower baseline scores than females.
# Age has a small, positive effect, but it's not very strong (~0.2 per year).
# BMI has no meaningful impact on baseline (coefficient ~0).
# Model is well-behaved → Good convergence, efficient sampling, and stable estimates.

baseline_lm2 = stan_glm(
  baseline ~ sex + age + BMI, 
  data = participants_nab
)
summary(baseline_lm2)

participants_baseline_comp = participants_nab[!is.na(participants_nab$baseline),]
participants_baseline_miss = participants_nab[is.na(participants_nab$baseline),]

baseline_imp_lm = predict(baseline_lm2, newdata = participants_baseline_miss) # miss predict
participants_baseline_miss$baseline = baseline_imp_lm

participants_baseline_comp$baseline_imp = predict(baseline_lm2, newdata = participants_baseline_comp)
# The next line just recreates the imputataion so that we can bind the datasets together
participants_baseline_miss$baseline_imp = predict(baseline_lm2, newdata = participants_baseline_miss)

participants_baseline_comp$baseline_imp_rand = participants_baseline_comp$baseline

## This time draw one point at random from the posterior distribution for each participant
## The output is a 1xdraws matrix, which we'll convert to a vector
baseline_rand_draw = posterior_predict(baseline_lm2, newdata = participants_baseline_miss, draws=1)
participants_baseline_miss$baseline_imp_rand = as.numeric(baseline_rand_draw)

## Now we can bind the two dataframes together and plot the randomly imputed / observed
## data against the deterministically fitted data

participants_baseline_imp = rbind(participants_baseline_comp, participants_baseline_miss)

ggplot(data = participants_baseline_imp, aes(x=baseline_imp, y=baseline_imp_rand, col = baseline_NA)) + 
  geom_point() +
  xlab("Regression fit") + 
  ylab("Observed value / randomly imputed baseline") +
  theme_bw()

participants_to_allocate <- participants_baseline_imp %>%
  mutate(baseline = baseline_imp_rand) %>%  # Replace baseline with baseline_imp_rand
  select(-baseline_imp_rand, -ID_NA, -sex_NA, -age_NA, -BMI_NA, -baseline_NA, -baseline_imp)  # Remove unwanted columns
  
# Save the dataset with comma separation
write.table(participants_to_allocate, "C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\CTs\\Summatives\\Summative 1\\participants_to_allocate.csv",
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE, fileEncoding = "UTF-8")

# Check column names
colnames(participants_to_allocate)
# [1] "ID"       "sex"      "age"      "BMI"      "baseline"

## Allocation
# set.seed(1)
# df$arm <- sample(rep(c("A", "B"), each = nrow(df)/2))
# table(df$arm)
df <- participants_to_allocate
table(df$sex) 
summary(df$age)
summary(df$BMI) # notice baseline isn't used - mention why if I use stratified sampling
str(df)

# split the data frame according to levels of factors - binning continuous variables
# strata: "sex"      "age"      "BMI"
# tbl_summary(df, by = "baseline")

# Make a copy of the original dataframe
df_category <- df  

# Define age bins
df_category$age_range[df_category$age < 50] <- "Under 50"  # This shouldn't happen based on your dataset (50-65)
df_category$age_range[df_category$age >= 50 & df_category$age < 55] <- "50-54"
df_category$age_range[df_category$age >= 55 & df_category$age < 60] <- "55-59"
df_category$age_range[df_category$age >= 60] <- "60-65"

# Convert age_range to factor
df_category$age_range <- factor(df_category$age_range, levels = c("50-54", "55-59", "60-65"))

# Define BMI bins
df_category$bmi_category[df_category$BMI < 25] <- "Below Average BMI"
df_category$bmi_category[df_category$BMI >= 25 & df_category$BMI < 29] <- "Average BMI"
df_category$bmi_category[df_category$BMI >= 29] <- "Above Average BMI"

# Convert bmi_category to factor
df_category$bmi_category <- factor(df_category$bmi_category, levels = c("Below Average BMI", "Average BMI", "Above Average BMI"))

# Convert Sex to factor
df_category$sex_category <- factor(df_category$sex, levels = c("M", "F"), labels = c("Male", "Female"))

# Print structure to verify changes
str(df_category)

# View tables to check counts
table(df_category$age_range)
table(df_category$bmi_category)
table(df_category$sex_category)

# View the updated dataframe
View(df_category)

df_category <- df_category %>%
  select(ID, sex_category, age_range, bmi_category, baseline, everything()) %>%  # Reorder columns
  select(-sex, -age, -BMI)  # Drop the original columns
View(df_category)

# library(dplyr)
# # split the data frame according to levels of factors
# strat_gen_age <- df_category %>%
#   group_split(sex_category, age_range) 
# strat_gen_age
# 
# group_sizes = sapply(
#   1:length(strat_gen_age),
#   function(i){
#     nrow(strat_gen_age[[i]])
#   }
# )
# group_sizes
?rep
nsample = nrow(df_category)
covmat=as.matrix(df_category[ ,2:4])
covwt=rep(1,3)/3
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
df_category$treat = res
View(df_category)

## You can now investigate the balance of the design as usual
table(df_category$treat)

# covmat <- covmat %>% 
#   mutate(across(everything(), as.numeric))
# covmat

# Display the number of randomized subjects at covariate factors
balance1 <- randbalance(res, covmat=covmat, ntrt=2, trtseq=c(0,1)) 
balance1

# Calculate the total imbalance of the allocation
totimbal(trt = res, covmat = covmat, covwt = covwt, 
         ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range")

df_category$treat <- ifelse(res == 0, "A", "B")
View(df_category)

table(df_category$treat)
colnames(df_category)
colnames(df)

df$arm <- df_category$treat
table(df$arm)
colnames(df)
View(df)

write.table(df, "C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\CTs\\Summatives\\Summative 1\\allocated_data.csv",
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE, fileEncoding = "UTF-8")

# Sex, Age, Baseline, BMI, Arm
ggplot(data = df, aes(x=age, y=baseline, col = sex)) + 
  geom_point() +
  xlab("age") + 
  ylab("Baseline") +
  ggtitle("Sex vs. Baseline") +
  theme_bw()

ggplot(data = df, aes(x=BMI, y=baseline, col = age)) + 
  geom_point() +
  xlab("BMI") + 
  ylab("Baseline") +
  ggtitle("Age vs. Baseline") +
  theme_bw()

ggplot(data = df, aes(x=age, y=baseline, col = BMI)) + 
  geom_point() +
  xlab("Age") + 
  ylab("Baseline") +
  ggtitle("BMI vs. Baseline") +
  theme_bw()

ggplot(data = df, aes(x=BMI, y=age, col = sex)) + 
  geom_point() +
  xlab("BMI") + 
  ylab("Age") +
  ggtitle("Sex vs. Age") +
  theme_bw()

ggplot(data = df, aes(x=age, y=baseline, col = arm)) + 
  geom_point() +
  xlab("Age") + 
  ylab("Baseline") +
  ggtitle("Arm vs. Baseline") +
  theme_bw()

ggplot(data = df, aes(x=age, y=BMI, col = sex)) + 
  geom_point() +
  xlab("age") + 
  ylab("BMI") +
  ggtitle("Sex vs. BMI") +
  theme_bw()

ggplot(data = df, aes(x=age, y=baseline, col = BMI)) + 
  geom_point() +
  xlab("Age") + 
  ylab("Baseline") +
  ggtitle("Age vs. BMI") +
  theme_bw()

## Trial Results
trial_results <- read.csv("C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\CTs\\Summatives\\Summative 1\\Trial Results.csv", 
                            stringsAsFactors = FALSE)

# Compute summary statistics
summary_stats <- trial_results %>%
  group_by(arm) %>%
  summarise(
    Sample_Size = n(),
    Mean = mean(outcome, na.rm = TRUE),
    SD = sd(outcome, na.rm = TRUE),
    SE_of_Mean = SD / sqrt(Sample_Size)
  )

# Print summary statistics table
print(summary_stats)

# Extract necessary values for t-test
group_A <- summary_stats %>% filter(arm == "A")
group_B <- summary_stats %>% filter(arm == "B")

# Function to compute pooled standard deviation and independent t-test
pooled_t_test <- function(n1, mean1, sd1, n2, mean2, sd2) {
  # Compute pooled standard deviation (sp)
  sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  
  # Compute t-statistic
  t_stat <- (mean1 - mean2) / (sp * sqrt(1/n1 + 1/n2))
  
  # Return results as a list
  return(list(
    pooled_SD = sp,
    t_statistic = t_stat
  ))
}

# Apply function using extracted values
result <- pooled_t_test(
  n1 = group_A$Sample_Size, mean1 = group_A$Mean, sd1 = group_A$SD,
  n2 = group_B$Sample_Size, mean2 = group_B$Mean, sd2 = group_B$SD
)

# Print results
print(result)

# A is the control
# B is the intervention
2 * pt(-abs(result$t_statistic), df = 440 - 5) # not significant

# small different in means as well - only 0.1 difference

# Function to compute confidence interval for two-sample t-test
confidence_interval <- function(nC, meanC, sdC, nT, meanT, sdT, alpha = 0.05) {
  # Compute degrees of freedom
  df <- nC + nT - 5 # CHECK DEGREES OF FREEDOM HERE!
  
  # Compute pooled standard deviation (sp)
  sp <- sqrt(((nC - 1) * sdC^2 + (nT - 1) * sdT^2) / df)
  
  # Compute t-critical value for given confidence level
  t_crit <- qt(1 - alpha/2, df)
  
  # Compute margin of error
  margin_of_error <- t_crit * sp * sqrt((1/nC) + (1/nT))
  
  # Compute confidence interval
  lower_bound <- (meanC - meanT) - margin_of_error
  upper_bound <- (meanC - meanT) + margin_of_error
  
  return(c(lower_bound, upper_bound))
}

confidence_interval(group_A$Sample_Size, group_A$Mean, group_A$SD, group_B$Sample_Size, group_B$Mean, group_B$SD) # [1] 0 is in this, supporting prior assessment

# Using Baseline values
trial_results$difference <- trial_results$outcome - trial_results$baseline
View(trial_results)

difference_summary_stats <- trial_results %>%
  group_by(arm) %>%
  summarise(
    Sample_Size = n(),
    Mean = mean(difference, na.rm = TRUE),
    SD = sd(difference, na.rm = TRUE),
    SE_of_Mean = SD / sqrt(Sample_Size)
  )
difference_summary_stats

group_difference_A <- difference_summary_stats %>% filter(arm == "A")
group_difference_B <- difference_summary_stats %>% filter(arm == "B")

# Apply function using extracted values
difference_result <- pooled_t_test(
  n1 = group_difference_A$Sample_Size, mean1 = group_difference_A$Mean, sd1 = group_difference_A$SD,
  n2 = group_difference_B$Sample_Size, mean2 = group_difference_B$Mean, sd2 = group_difference_B$SD
)

print(difference_result)
2 * pt(-abs(difference_result$t_statistic), df = 440 - 5) # [1] not significant but has reduced a lot

# ANCOVA
trial_results = trial_results[,-8] # remove "difference"
colnames(trial_results)

str(trial_results)

# trial_results$arm <- factor(trial_results$arm, labels = c("Placebo", "Captopril"))
ggplot(trial_results, aes(y = arm, x = baseline, fill = arm)) +
  geom_boxplot(width = 0.5, color = "black", size = 0.5) +  # Black outline around each box
  scale_fill_manual(values = c("cyan", "red")) +  # Fill colors for groups
  scale_color_manual(values = c("black", "black")) +  # Ensure black outline for both groups
  labs(title = "", 
       x = "Baseline Blood Pressure (mmHg)", 
       y = "Group") +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    axis.text.y = element_text(size = 14, color = "black"),  # Increase y-axis label size
    axis.text.x = element_text(size = 12, color = "black")   # Increase x-axis label size
  )

lm_participants = lm(outcome ~ baseline + arm, data = trial_results)
sum <- summary(lm_participants)

vif(lm_participants)

mean_value = sum$coefficients[3,1]
margin_of_error = sum$coefficients[3,2]
t_crit <- qt(1 - 0.05/2, 440-5)

min = mean_value - t_crit*margin_of_error
max = mean_value + t_crit*margin_of_error
min
max

resid_capt = resid(lm_participants)
trial_results$resid= resid_capt
ggplot(data = trial_results, aes(x=baseline, y=resid, col=arm)) +
  geom_point() +
  geom_hline(yintercept=0)+
  xlab("Baseline")+
  ylab("Residual")+theme_bw()

lm_part_int = lm(outcome ~ arm + baseline + baseline:arm, data = trial_results)
tbl_regression(lm_part_int)
# look at the arm:baseline row for p-value significance: 0.97

# Therefore we can be confident that there is no need to fit unequal slopes for this dataset. This fits with our
# earlier conclusion (from inspecting the residuals) that just including first order terms is fine.

trial_results = trial_results[,-8] # remove "resid"
colnames(trial_results)
lm_multiple = lm(outcome ~ baseline + arm + sex + age + BMI, data = trial_results) # BMI significance
summary(lm_multiple)
tbl_regression(lm_multiple)
vif(lm_multiple)


lm_sex_BMI = lm(outcome ~ baseline + arm*sex + arm*BMI, data = trial_results)
summary(lm_sex_BMI)
tbl_regression(lm_sex_BMI)
vif(lm_sex_BMI)

# Fit ANCOVA model with interactions
ancova_model <- lm(outcome ~ arm * baseline + arm * age + arm * BMI + arm * sex, data = trial_results)

# Display results
tbl_regression(ancova_model)
vif(ancova_model)

# Fit ANCOVA model with interactions
final_model <- lm(outcome ~ sex + age + baseline + BMI + arm + baseline:sex, data = trial_results)

# Display results
tbl_regression(final_model)
vif(final_model)

lm_multiple = lm(outcome ~ baseline + arm + sex + age + BMI + sex:baseline, data = trial_results) # BMI significance
summary(lm_multiple)
tbl_regression(lm_multiple)
vif(lm_multiple)

lm_multiple = lm(outcome ~ baseline + arm + age + BMI + sex:baseline, data = trial_results) # BMI significance
vif(lm_multiple)

summary(lm_multiple)
tbl_regression(lm_multiple)

resid_capt = resid(lm_multiple)
trial_results$resid= resid_capt
ggplot(data = trial_results, aes(x=baseline, y=resid, col=arm)) +
  geom_point() +
  geom_hline(yintercept=0)+
  xlab("Baseline")+
  ylab("Residual")+theme_bw()

# Given values
estimate <- 0.52781  # Coefficient estimate
standard_error <- 0.81179  # Standard error
df <- 435  # Degrees of freedom

# Compute the critical t-value for 95% confidence interval (two-tailed)
t_critical <- qt(0.975, df)

# Compute the margin of error
margin_of_error <- t_critical * standard_error

# Compute the confidence interval
ci_lower <- estimate - margin_of_error
ci_upper <- estimate + margin_of_error

# Print the confidence interval
cat("95% Confidence Interval: (", round(ci_lower, 3), ",", round(ci_upper, 3), ")\n")

