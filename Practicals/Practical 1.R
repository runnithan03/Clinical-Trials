install.packages(c("naniar", "smdi", "tidyverse", "visdat", "rstanarm", "medicaldata","ggplot2","gridExtra","knitr", "parallel"))


require(medicaldata)
library(naniar)
library(smdi)
library(tidyverse)
library(visdat)
library(rstanarm)
library(ggplot2)
library(gridExtra)
library(knitr)
library(parallel)

# A.1
data("supraclavicular")

# details of the study
?supraclavicular

# The primary outcome variable is time to 4-nerve sensory block
# This has the column name onset_sensory and the baseline covariates are gender, BMI and age.

# Creating the new dataframe:
sup_df = supraclavicular[ ,c(2:5, 9)]
sup_df = sup_df%>% mutate(across(c(group, gender), factor))

# A.1.1.2 Lung Cancer Dataset
# This is the smdi_data dataset from the package smdi. The intervention group variable is exposure. 
# The outcome data are survival (or ‘time-to-event’) data: in this study, the ‘event’ is death and the follow-up period is five (probably years, but it doesn’t say!).

smdi_data$status = factor(smdi_data$status)

# This is a synthetic (or possibly synthetically altered) dataset designed for use with missing data exploration, 
# so we have the advantage that the help file tells about the missingness mechanisms (more on that soon).

# A.1.1.3 Obstetrics and Periodontal Therapy
# A.2
?opt

opt_df = opt[ ,c(1:4, 10:22, 72)]

# retains most of the baseline covariates, intervention group and outcome variable (birthweight). 
# There is some ‘messiness’ in the data, for example for some variables missing values are recorded as "", rather than NA. 
# We’ll sort that now by running the following command:

opt_df = opt_df %>% 
  mutate(across(c(Use.Tob, Use.Alc, Drug.Add), 
                gsub, pattern = "^ +$", replacement = NA))%>%
  mutate(across(c(Use.Tob, Use.Alc, Drug.Add), factor))

# A.2: Understanding the Patterns of Missingness
# A.2.1 Visualising missingness

# this colours data by type (and missing data in grey)
# sort_type = F means it keeps to the original order
vis_dat(sup_df, sort_type=F)

# The functions gg_miss_var, gg_miss_case and gg_miss_upset (all from naniar) allow us to quickly see how much missing data there is, 
# either by variable, by case or by intersection. Look at the help files to find out what the arguments do. 
# We see that in supraclavicular there are three missing entries for BMI, and all other variables are complete.
gg_miss_var(sup_df)

# Tells us how much missing data there is for each participant
gg_miss_case(sup_df)
md.pattern(sup_df)

# A.3
vis_dat(smdi_data)
# there is missingness for 3 variables

# the naniar functions show us how this is spread across the variables and participants
gg_miss_var(smdi_data, show_pct=T)
gg_miss_case(smdi_data, order_cases = T)

# the function md.pattern shows us how many participants have each combination of missing data:
md.pattern(smdi_data, rotate.names = T)

# Using vis_dat we see that there is a lot of missing data in the opt dataset, 
# and that some variables are particularly badly affected.
vis_dat(opt_df)

gg_miss_var(opt_df, show_pct=T)
gg_miss_case(opt_df, order_cases = T, show_pct=T)

# Not many individual participants have more than about 25% of their data missing, 
# but there are some variables that are missing for nearly all participants (more on this later).

# We may also want to visualise or summarise missingness in a table, and there are various ways to do this, 
# for example miss_case_summary, miss_case_table and miss_var_summary.

# This shows how many variables (both in number and as a percentage of the total number of variables) are missing, 
# for each case (in decreasing order of missingness).
miss_case_summary(opt_df)
miss_case_summary(smdi_data)
miss_case_summary(supraclavicular)

# miss_case_table tabulates the number of cases with n missing variables, for n = 1,..,
miss_case_table(opt_df)
miss_case_table(smdi_data)
miss_case_table(supraclavicular)

# Now we see how many (and what percentage) of cases are missing for each variable.
miss_var_summary(opt_df)
miss_var_summary(smdi_data)
miss_var_summary(supraclavicular)

# A.2.3 Mechanisms of missingness
# A.2.3.1 Missing completely at random (MCAR)

# A.2.3.2 Missing at random (MAR)

# A.2.3.3 Missing not at random (MNAR)

# A.3 Exploring the relationship between missingness and other variables

# In a clinical trial, ultimately what we care about is whether the missingness has changed our understanding of the outcome variable(s).
# In the package naniar you can create a copy of the data frame containing just NA and !NA, indexing exactly where the missing values are.

as_shadow(sup_df)

# To avoid duplication of the original names, the column names are suffixed by ’_NA’. 
# This means it can be appended to the original data frame to create what the authors of naniar call a nabular object.
nabular(sup_df)

# This is useful because we can investigate the values of the actual data while conditioning on the missingness. 
# We can summarize the outcome distribution according to whether a variable is missing or observed:
sup_nab <- nabular(sup_df)
sup_nab %>%
  group_by(bmi_NA) %>%
  summarise_at(.vars = "onset_sensory",
               .funs = c("mean", "sd", "var", "min", "max"),
               na.rm = TRUE)

# and visualise the outcome distribution for missing and non-missing values of a covariate:
ggplot(sup_nab,
         aes(x = onset_sensory,
             fill = bmi_NA)) + 
  geom_histogram()

#We could also explore whether the missingness is related to the values of the other covariates, for example gender (in this data 0=female, 1=male)
ggplot(sup_nab,
       aes(x = gender,
           fill = bmi_NA)) + 
  geom_bar()

# A.5
# This isn’t a great example because there are only three missing BMI values, but you can probably guess what you’ll be doing next…
smdi_nab = nabular(smdi_data)

# This is determining if the missing values are MAR, MCAR or MNAR

# First by ecog_cat. Since there is a lot more missing data here, I’ve chosen to compare the distributions side-by-side using facet_wrap. 
# This also allows me to colour by status (although it is sort of obvious from the large peak at eventtime=5 that those people are still alive). 
# Setting scales = "free_y" makes it easier to compare the overall shape of the distribution. 

smdi_nab %>%
  group_by(ecog_cat) %>%
  summarise_at(.vars = "eventtime",
               .funs = c("mean", "sd", "var", "min", "max"),
               na.rm = TRUE)

ggplot(smdi_nab,
       aes(x = eventtime,
           fill = status)) + 
  geom_histogram() + facet_wrap(~ecog_cat_NA, scales = "free_y")

# Next by egfr_cat:
smdi_nab %>%
  group_by(egfr_cat) %>%
  summarise_at(.vars = "eventtime",
               .funs = c("mean", "sd", "var", "min", "max"),
               na.rm = TRUE)

ggplot(smdi_nab,
       aes(x = eventtime,
           fill = status)) + 
  geom_histogram() + facet_wrap(~egfr_cat_NA, scales = "free_y")

# finally, by pdl1_num:
smdi_nab %>%
  group_by(pdl1_num) %>%
  summarise_at(.vars = "eventtime",
               .funs = c("mean", "sd", "var", "min", "max"),
               na.rm = TRUE)

ggplot(smdi_nab,
       aes(x = eventtime,
           fill = status))+
  geom_histogram() + facet_wrap(~pdl1_num_NA, scales = "free_y")

# It appears that pdl1_num values are likely to be MNAR, since the outcome distribution looks different (proportionally more early deaths) for the missing values.
# We could make many more plots of this type: for example, by plotting the distributions of other variables rather than the outcome.
# For example, we could look at how the missingness of each variable relates to smoking category:

library(gridExtra)
ecog_sm = ggplot(smdi_nab,
                 aes(x = smoking_cat,
                     fill = ecog_cat_NA))+ geom_bar() 


egfr_sm = ggplot(smdi_nab,
                 aes(x = smoking_cat,
                     fill = egfr_cat_NA))+geom_bar() 

pdl1_sm = ggplot(smdi_nab,
                 aes(x = smoking_cat,
                     fill = pdl1_num_NA))+geom_bar()

grid.arrange(ecog_sm, egfr_sm, pdl1_sm, nrow=1)
# So, it looks likely that the missingness of egfr_cat depends on smoking status, 
# but less likely that the missingness of the other two variables does.

# These tests can also be ran on opt but slightly altered... TBD myself
opt_tmp = opt_df[ ,-c(9,12,14,17)]

# A.3.1 Statistical summaries of the effect of missingness

# A.3.1.1 Hotelling’s multivariate t-test - tests for MCAR vs. MNAR and MAR

sum(is.na(sup_df))

# A.6
smdi_hotelling(sup_df)
# For this data, we find that there is insufficient evidence to reject the null hypothesis that the BMI data are MCAR.

smdi_hotelling(smdi_data)
# It appears that ecog_cat may be MCAR, but for egfr_cat and pdl1_num there is sufficient evidence to suggest MAR or MNAR.

smdi_hotelling(opt_tmp)
# These all appear to have a significant departure from MCAR.
# Caution: The power of this test (and others like it) can be strongly influenced by sample size, 
# so it is sensible to combine it with a more detailed approach.

# A.3.1.2 Absolute standardised mean difference

# The function smdi_asmd from smdi creates an asmd object, which has several parts to it. 
# Note that we set includeNA=T so that we can see the effect of missingness, as well as observed values, of other variables. 
# This won’t make a difference for sup_df, but it will in the exercise.
asmd_sup = smdi_asmd(sup_df, includeNA=T)

# There is also a Table 1, so called because a summary table of this nature should always be included when summarising a dataset in terms of the difference between two groups. 
# This is formatted a little strangely, as it is designed for use in printed works. 
# This table includes the result of a statistical test (by default a chi-squared test) showing whether the differences are statistically significant.
kable(asmd_sup$bmi$asmd_table1)

# Finally there is a plot showing each ASMD
asmd_sup$bmi$asmd_plot

# We see that although the ASMD values are quite large (much bigger than the advised 0.1), because there are a very small number of them they are not statistically significant.

# Exercise A.7
# Investigate the ASMD for our datasets smdi_data and opt_tmp. Of the partially observed covariates, which seem to be most strongly related to other covariates? Do any seem to be MCAR?
# Make sure you remove any participant ID variables, since we don’t want to include those in our analysis! Again, it might be a good idea to pair up!

asmd_opt = smdi_asmd(opt_tmp[ ,-1], includeNA=T)
asmd_opt

# These values are all well above 0.1, so it seems unlikely that the data are MCAR.

# We can look in a little more detail, first at BMI:
asmd_opt$BMI$asmd_plot
kable(asmd_opt$BMI$asmd_table1)

# We see that there appears to be a strong relationship between clinic and missingness of BMI. In particular, a disproportionately high number of missing BMI values seem to be from the New York (NY) clinic. 
# There also appears to be less missingness for those with public.asstce=1 (those for whom the government paid for the delivery).

# Next use Use.Alc
kable(asmd_opt$Use.Alc$asmd_table1)
asmd_opt$Use.Alc$asmd_plot

# This time the two most strongly related covariates are Use.Tob and Drug.Add. This is perhaps not surprising because they are quite similar variables. #
# Notice that clinic still has a strong effect (the x-axis goes a lot higher on this plot than the one for BMI). 
# From the table, we see that again a disproportionately high number of missing values come from the NY clinic, and that missingness in Use.Alc is much more likely if Use.Tob and/or Drug.Add are missing.

# A.3.2 Modelling missingness
# One of the most common ways to model missingness is using logistic regression.

# We’ll use our nabular objects for this, since we already have variables denoting missingness. 
# For example, we can see whether missingness of BMI in the supraclavicular dataset relates to any of the other variables.

sup_nab = nabular(sup_df)
sup_glm = glm(bmi_NA ~ group + gender + age + onset_sensory, 
              data = sup_nab, 
              family = binomial(link = "logit"))
summary(sup_glm)

# In this case none of the other variables are significant, so it would be reasonable to proceed as though the missingness of BMI is MCAR (unless there is expert knowledge to suggest that missingness might be linked to some other observed factor). 
# This agrees with what we found in Section A.3.1.
# A word of caution: the default in most R functions, particularly for plotting or fitting models, is to remove all rows with any missingness. 
# In a situation where there are missing values for multiple variables, particularly if the missingness is related, this could in itself introduce bias.

# Exercise A.8

# Lung Cancer Data:
ecog_NA_glm1 = glm(ecog_cat_NA ~ exposure + age_num + female_cat + smoking_cat + physical_cat + alk_cat + histology_cat + ses_cat + copd_cat, 
                   data=smdi_nab, family = binomial(link="logit"))
summary(ecog_NA_glm1)

# Now including the other variables with missing values:
  
ecog_NA_glm2 = glm(ecog_cat_NA ~ exposure + age_num + female_cat + smoking_cat + physical_cat + alk_cat + histology_cat + ses_cat + copd_cat + egfr_cat + pdl1_num, 
                     data=smdi_nab, family = binomial(link="logit"))
summary(ecog_NA_glm2)

# None of the variables are significant, so it is probably reasonable to suppose that ecog_cat is MCAR (indeed the help file tells us it is).

# Now for egfr_cat:
egfr_NA_glm1 = glm(egfr_cat_NA ~ exposure + age_num + female_cat + smoking_cat + physical_cat + alk_cat + histology_cat + ses_cat + copd_cat, 
                     data=smdi_nab, family = binomial(link="logit"))

# This time there is a definite relationship between missingness and the values of the other covariates, suggesting an MAR or MNAR mechanism.
# Finally, let’s model missingness in pdl1_num:
summary(egfr_NA_glm1)

# Similar idea here...
egfr_NA_glm2 = glm(egfr_cat_NA ~ exposure + age_num + female_cat + smoking_cat + physical_cat + alk_cat + histology_cat + ses_cat + copd_cat + ecog_cat + pdl1_num, 
                   data=smdi_nab, family = binomial(link="logit"))
summary(egfr_NA_glm2)

# This time there is a definite relationship between missingness and the values of the other covariates, suggesting an MAR or MNAR mechanism.

# Finally, let’s model missingness in pdl1_num:
pdl1_NA_glm1 = glm(
    pdl1_num_NA ~ exposure + age_num + female_cat + smoking_cat + physical_cat + 
      alk_cat + histology_cat + ses_cat + copd_cat, 
    data=smdi_nab, family = binomial(link="logit"))
summary(pdl1_NA_glm1)

pdl1_NA_glm2 = glm(pdl1_num_NA ~ exposure + age_num + female_cat + smoking_cat + physical_cat + alk_cat + histology_cat + ses_cat + copd_cat + ecog_cat + egfr_cat, 
                   data=smdi_nab, family = binomial(link="logit"))
summary(pdl1_NA_glm2)

# This time there are far fewer significantly related variables, but we can again be confident that the mechanism isn’t MCAR.
# Frustratingly, the only way we could determine whether the mechanisms for egfr_cat and pdl1_num was MAR or MNAR would be to measure (or otherwise procure) some of the missing data. 
# We simply do not have the necessary information to work out which is the case. If this were a real trial, we would now talk at length with the experts/clinicians, who will have a much better understanding of the probable causes of missingness.

# Now we’ll do the same with opt_tmp. We have missingness in several variables: BMI, Use.Tob, Use.Alc, Drug.Add, Birthweight. A model fit to the data in R will remove cases with any of these missing 
# (which for this dataset might mean a lot are removed!). So, as well as building models with all variables, we can build models using only the fully observed data, and use the _NA variables from nabular as covariates:
glm_opt_BMI_allvar = glm(
  BMI_NA ~ Clinic + Group + Age +  
    Education + Public.Asstce + Hypertension + Diabetes + Use.Tob + Use.Alc + Drug.Add + Prev.preg, 
  data = nabular(opt_tmp),
  family = binomial(link = "logit")) 
summary(glm_opt_BMI_allvar)

# Complete cases only:

glm_opt_BMI_comp = glm(
  BMI_NA ~ Clinic + Group + Age +  
    Education + Public.Asstce + Hypertension + Diabetes + Prev.preg, 
  data = nabular(opt_tmp),
  family = binomial(link = "logit")) 
summary(glm_opt_BMI_comp)

# Using missingness of incomplete variables in model:

glm_opt_BMI_NA = glm(
  BMI_NA ~ Clinic + Group + Age +  
    Education + Public.Asstce + Hypertension + Diabetes + Use.Tob_NA + Use.Alc_NA + Drug.Add_NA + Prev.preg, 
  data = nabular(opt_tmp),
  family = binomial(link = "logit")) 
summary(glm_opt_BMI_NA)

# From these three models it appears that BMI is much more likely to be missing for those from the NY clinic (we already knew this from our previous investigations!):
ggplot(data=nabular(opt_tmp), aes(x=Clinic, fill = BMI_NA)) +
  geom_bar()

# Perhaps the most concerning missing data in the opt dataset is in the outcome Birthweight.

glm_opt_BW_allvar = glm(
  Birthweight_NA ~ Clinic + Group + Age +  
    Education + Public.Asstce + Hypertension + Diabetes + BMI + Use.Tob + Use.Alc + Drug.Add + Prev.preg, 
  data = nabular(opt_tmp),
  family = binomial(link = "logit")) 
summary(glm_opt_BW_allvar)

glm_opt_BW_comp = glm(
  Birthweight_NA ~ Clinic + Group + Age +  
    Education + Public.Asstce + Hypertension + Diabetes + Prev.preg, 
  data = nabular(opt_tmp),
  family = binomial(link = "logit")) 
summary(glm_opt_BW_comp)

glm_opt_BW_NA = glm(
  Birthweight_NA ~ Clinic + Group + Age +  
    Education + Public.Asstce + Hypertension + Diabetes + BMI_NA + Use.Tob_NA + Use.Alc_NA + Drug.Add_NA + Prev.preg, 
  data = nabular(opt_tmp),
  family = binomial(link = "logit")) 
summary(glm_opt_BW_NA)

ggplot(data = nabular(opt_tmp), aes(x=Diabetes, fill = Birthweight_NA)) +
  geom_bar()

# A.4 What to do about missing data?!
# There are 2 options: Discard some (non-missing) data and Add in (‘impute’) some synthetic data


# A.4.1 Discarding some (non-missing) data

# A.4.1.1 Complete case analysis
# The very simplest thing we can do is to discard the data for any participant who has some missing data. This is called a complete-case analysis, because we only analyse data for participants whose data are complete.

# There are two main problems with this:
  
# If the missing data are not MCAR, then this can induce bias.
# This approach can drastically reduce the amount of data

# Exercise A.9 
# Number of complete and incomplete cases:
nrow(na.omit(sup_df))
nrow(sup_df)
nrow(na.omit(smdi_data))
nrow(smdi_data)
nrow(na.omit(opt_tmp))
nrow(opt_tmp)

# A.4.2 Imputing data

# A.4.2.1 Mean imputation
# If the data are not MCAR, bias is introduced
# The sample standard deviation is reduced
# Relationships between this and other variables are distorted

# A.4.2.2 Imputing using logic

# Sometimes there is missingness in a dataset that we can fill in using information about how that variable relates to other variables. For example, in the opt data, consider the two columns Use.Tob and BL.Cig.Day:
# Use.Tob: Self-reported participant history of tobacco use, factor; Yes, No; Blank = Missing
# BL.Cig.Day: Self-reported number of cigarettes per day for those with tobacco use history, numeric, range: 1-30; Blank = Missing (variable 16= Yes or blank) or non-smoker (variable 16 = No)`

# We therefore need to condition on the ‘parent’ variable, rather than replace all missing values of the ‘child’ variable. For example:
opt_df_imp = opt_df
opt_df_imp$BL.Cig.Day[opt_df_imp$Use.Tob=="No "] = 0

# By doing this we have ‘fixed’ 704 missing values.
# Similarly we cand fix those for BL.Drks.Day and N.prev.preg:
opt_df_imp$BL.Drks.Day[opt_df_imp$Use.Alc=="No "] = 0
opt_df_imp$N.prev.preg[opt_df_imp$Prev.preg=="No "] = 0

# We can see that this has very much improved our situation! We could do something similar for BL.Diab.Type too since this is linked to Diabetes.
vis_dat(opt_df, sort_type = F)
vis_dat(opt_df_imp, sort_type = F)

# A.4.2.3 Imputation using a regression model
# Using STAN functions from rstanarm

# Let’s suppose we want to use regression to impute values for pdl1_num
colnames(smdi_data)
pdl1_lm = stan_glm(
  pdl1_num ~ exposure + age_num + female_cat + smoking_cat + physical_cat + alk_cat + 
    histology_cat + ses_cat + copd_cat + eventtime + status + ecog_cat + egfr_cat, 
  data = smdi_data
)
summary(pdl1_lm)

# A.11: This approach is going to fail...
# The model involves variables that also have missingness (ecog_cat1 and egfr_cat), and so for any case with either of those missing, we won’t be able to impute a value.

# There are a couple of possibilities:
  
# Use the missingness of those values as an input (by creating a nabular object)
# Remove those variables from the model
# Work iteratively, generating temporary imputed values and cycling round the variables with missingness
# Because all our investigations suggest that ecog_cat1 is more-or-less unrelated to anything, we will remove it from the model. However, because egfr_cat is quite close to being significant in the model, 
# we will keep it in but use missingness in the model.

smdi_nab = nabular(smdi_data)
str(smdi_nab)

pdl1_lm2 = stan_glm(
  pdl1_num ~ exposure + age_num + female_cat + smoking_cat + physical_cat + alk_cat + 
    histology_cat + ses_cat + copd_cat + eventtime + status + egfr_cat_NA, 
  data = smdi_nab
)


summary(pdl1_lm2)

# We can see from the means and standard deviations that most of the coefficients are not even close to ‘significant’, 
# but if we want to investigate some more closely we can visualise them using the plot function (separately in this case because they are quite far apart numerically)
plot(pdl1_lm2, plotfun="areas", prob = 0.95, pars = c("exposure"))
plot(pdl1_lm2, plotfun="areas", prob = 0.95, pars = c("age_num"))

# We will now use this second model to impute values for pdl1_num.

## Split the data according to whether egfr_cat is missing
smdi_pdl1_comp = smdi_nab[!is.na(smdi_nab$pdl1_num),]
smdi_pdl1_miss = smdi_nab[is.na(smdi_nab$pdl1_num),]

## Use the GLM to fit values to egfr_cat
pdl1_imp_lm = predict(pdl1_lm2, newdata = smdi_pdl1_miss)
smdi_pdl1_miss$pdl1_num = pdl1_imp_lm

## Join the data back together again (in a different order, but it doesn't matter)
smdi_imp = rbind(smdi_pdl1_miss, smdi_pdl1_comp)

# This isn’t a very realistic use of a regression model.

# We can see this by comparing the observed values of pdl1_num to their fitted values (which we would have imputed had they been missing).

smdi_pdl1_comp$pdl1_imp = predict(pdl1_lm2, newdata = smdi_pdl1_comp)

# The next line just recreates the imputataion so that we can bind the datasets together
smdi_pdl1_miss$pdl1_imp = predict(pdl1_lm2, newdata = smdi_pdl1_miss)

smdi_pdl1_imp = rbind(smdi_pdl1_comp, smdi_pdl1_miss)

ggplot(data = smdi_pdl1_imp, aes(x=pdl1_imp, y=pdl1_num, col = pdl1_num_NA)) + 
  geom_point() +
  xlab("Regression fit") + 
  ylab("Observed value / imputed value") +
  theme_bw()

# A more realistic approach would be to sample one value from the posterior distribution for each point (this is why we are using rstanarm!)
smdi_pdl1_comp$pdl1_imp_rand = smdi_pdl1_comp$pdl1_num

## This time draw one point at random from the posterior distribution for each participant
## The output is a 1xdraws matrix, which we'll convert to a vector
pdl1_rand_draw = posterior_predict(pdl1_lm2, newdata = smdi_pdl1_miss, draws=1)
smdi_pdl1_miss$pdl1_imp_rand = as.numeric(pdl1_rand_draw)

## Now we can bind the two dataframes together and plot the randomly imputed / observed
## data against the deterministically fitted data

smdi_pdl1_imp = rbind(smdi_pdl1_comp, smdi_pdl1_miss)

ggplot(data = smdi_pdl1_imp, aes(x=pdl1_imp, y=pdl1_imp_rand, col = pdl1_num_NA)) + 
  geom_point() +
  xlab("Regression fit") + 
  ylab("Observed value / randomly imputed value") +
  theme_bw()

# This imputed data looks much more representative of the actual dataset.

# Exercise A.12 Use regression imputation to impute values for BMI in opt_tmp - TBD...

### Do this before doing PRACTICAL 2
opt_nab = nabular(opt_tmp)
opt_tmp_bmi_comp = opt_nab[!is.na(opt_nab$BMI),]
opt_tmp_bmi_miss = opt_nab[is.na(opt_nab$BMI),]

bmi_lm = stan_glm(
  BMI ~. , 
  data = opt_tmp
)
summary(bmi_lm)

bmi_lm2 = stan_glm(
  BMI ~ PID_NA + Clinic_NA + Group_NA + Age_NA + Education_NA + Public.Asstce_NA + Hypertension_NA + Diabetes_NA + BMI_NA + Use.Tob_NA + Use.Alc_NA + Drug.Add_NA + Prev.preg_NA + Birthweight_NA, 
  data = opt_nab
)

str(opt_nab)

opt_tmp_bmi_comp$pdl1_imp = predict(bmi_lm2, newdata = opt_tmp_bmi_comp)
# The next line just recreates the imputataion so that we can bind the datasets together
opt_tmp_bmi_miss$pdl1_imp = predict(pdl1_lm2, newdata = opt_tmp_bmi_miss)

opt_imp = rbind(opt_tmp_bmi_miss, opt_tmp_bmi_comp)

opt_tmp$BMI = opt_tmp_bmi_comp$pdl1_num

## This time draw one point at random from the posterior distribution for each participant
## The output is a 1xdraws matrix, which we'll convert to a vector
bmi_rand_draw = posterior_predict(pdl1_lm2, newdata = smdi_pdl1_miss, draws=1)
smdi_pdl1_miss$pdl1_imp_rand = as.numeric(bmi_rand_draw)

## Now we can bind the two dataframes together and plot the randomly imputed / observed
## data against the deterministically fitted data

smdi_pdl1_imp = rbind(smdi_pdl1_comp, smdi_pdl1_miss)

ggplot(data = smdi_pdl1_imp, aes(x=pdl1_imp, y=pdl1_imp_rand, col = pdl1_num_NA)) + 
  geom_point() +
  xlab("Regression fit") + 
  ylab("Observed value / randomly imputed value") +
  theme_bw()

# A.4.3 Multiple imputation
# Having reached the point of acknowledging the randomness needed in imputation, a natural next step would be to draw multiple values from the posterior distribution, 
# rather than just one. This leads to a method called multiple imputation, which is (arguably) the most widely used approach to imputing missing data




