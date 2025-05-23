## Allocation
install.packages(c("medicaldata", "ggplot2", "gtsummary", "Minirand", 
                   "blockrand", "dplyr", "randomizeR","car"))
library(medicaldata)
library(ggplot2)
library(gtsummary)
library(Minirand)
library(blockrand)
library(dplyr)
library(randomizeR)
library(car)

# Licorice gargle dataset
data("licorice_gargle")

?licorice_gargle

# baseline characteristics used in allocation
str(licorice_gargle)

lic_garg = licorice_gargle[ ,1:8]
# vector of names of columns to be coerced to factor
cols <- c("preOp_gender", "preOp_asa",  
          "preOp_mallampati", "preOp_smoking", "preOp_pain", "treat")
# convert each of those columns to factors
lic_garg[cols] <- lapply(lic_garg[cols], factor) 

# Check the result:
str(lic_garg)

# Demographic tables
# Licorice gargle data contains the allocation used in the trial, 
# so we can create a table summarising the participants in each group.
tbl_summary(lic_garg, by = "treat")

# We can also save this as an object and access the individual sub-tables, for example:
rb_tab = randbalance(
  trt = lic_garg$treat, 
  covmat = lic_garg[,-8], 
  ntrt=2, 
  trtseq = c("0", "1"))
rb_tab$preOp_gender

# There are also packages with functions to output these demographic tables formatted for use in latex documents, for example atable.

# Binning continuous variables
# Binning in the context of clinical trials is a method used to group continuous data into discrete categories (bins). 
lic_garg$preOp_age

lic_garg$age[lic_garg$preOp_age < 50] <- "Under 50"
lic_garg$age[lic_garg$preOp_age >= 50 & lic_garg$preOp_age < 70] <- "50 to 70"
lic_garg$age[lic_garg$preOp_age >= 70] <- "70 plus"
lic_garg$age = factor(lic_garg$age, levels = c("Under 50", "50 to 70", "70 plus"))

table(lic_garg$age)


lic_garg$BMI[lic_garg$preOp_calcBMI < 25] <- "medium_or_low"
lic_garg$BMI[lic_garg$preOp_calcBMI >= 25] <- "high"
lic_garg$BMI = factor(lic_garg$BMI, levels = c("medium_or_low", "high"))

table(lic_garg$BMI)

lg_df = lic_garg[,c(1,2,5,6,7,9,10,8)] # replace preOp BMI and age with categories
lg_df

# A Measure of Imbalance
imbalance = function(
    df,   # participant data frame with allocation column included
    alloc # name of allocation column
){
  alloc_vec = as.factor(df[ ,names(df)==alloc])
  alloc_lev = levels(alloc_vec) # how the treatment groups are coded
  n1 = nrow(df[df[alloc]==alloc_lev[1],])
  n2 = nrow(df[df[alloc]==alloc_lev[2],])
  abs(n1-n2)
}

# imbalance in the allocation recorded in the lic_garg dataset
imbalance(df = lic_garg, alloc = "treat")

# Allocation Methods

# Simple Random Allocation

# In simple random allocation, each participant is allocated to one of the two trial arms with equal probability.
srs = function(
    df, # DF should be the participant data frame. 
    # A column 'treat' will be added
    levels = c("0", "1") # Levels of treat factor
){
  n = nrow(df) # number of rows / participants
  # Create a new column 'treat'
  df$treat = rep(NA, n)
  # work through the rows, randomly allocating patients with probably 1/2
  for (i in 1:n){
    df$treat[i] = sample(levels, size=1, prob = c(0.5, 0.5)) # uses sample so a seed should be set
  }
  df$treat = as.factor(df$treat)
  df
}

lg_srs = srs(
  df = lg_df[,-8], # remove treatment column
  levels = c("T", "C")
)

# Balance Table
tbl_summary(lg_srs, by = "treat")

# Imbalance command
imbalance(lg_srs, alloc = "treat")

# Randomly permuted blocks (RPBs)

# The function for generating RPB designs is also called blockrand. 
# The default for the function blockrand is that it will randomly vary block length within {2,4,6,8}.

blockrand(n=100)

# allocation for the licorice_gargle data (remember our data frame is now lg_df) and produce the balance table
?blockrand

rpb_lg = blockrand(n=235, levels = c("T", "C")) # Remember to define the levels...
# Notice that this doesn’t have 235 rows: the blockrand function will always finish after a whole block.
# block.id labels each block uniquely.
# Each participant is assigned to a specific block before being randomized within that block.
# Why is it needed?
# Tracking assignment → Helps maintain structure and avoid misallocation.
# Ensuring equal distribution → Each block contains a fixed number of participants.

# block.size is the size of each block (e.g., 4, 6, 8 participants).

# Let’s add this to our participant data to create a new data frame lg_rpb:
# create the new data frame, a copy of lg_df
lg_rpb = lg_df  

# Replace the original treat column with the RPB treatment column
# Using only the first 235 allocations
lg_rpb$treat = rpb_lg$treatment[1:235]

# balance table:
tbl_summary(lg_rpb, by = "treat")

# Imbalance command
imbalance(lg_rpb, alloc = "treat")

# The package plotblockrand outputs PDFs of randomization cards, ready to be printed and put into envelopes!

# Biased Coin Design

biased_coin = function(
    data,
    levels = c("T", "C"),
    p=2/3 # makes coin biased...
){
  Dn = 0 # starting value of imbalance
  n = nrow(data)
  alloc = rep(NA, n)
  
  for (i in 1:n){
    if (Dn==0){ # equally balanced
      alloc[i] = sample(levels, size=1, prob=c(0.5, 0.5) )
    } else if(Dn<0){ # More allocations to levels[2] up to this point
      alloc[i] = sample(levels, size=1, prob=c(p, 1-p) )
    } else if(Dn>0){ # More allocations to levels[1] up to this point
      alloc[i] = sample(levels, size=1, prob=c(1-p, p) )
    }
    # Compute imbalance at this stage
    alloc_to_n = alloc[1:i]
    Dn = sum(alloc_to_n==levels[1]) - sum(alloc_to_n == levels[2])
  }
  data$treat = as.factor(alloc)
  data
}

# try p = 0.9
lg_bc1 = biased_coin(
  data=lg_df[,-8],
  p=0.9
)
tbl_summary(lg_bc1, by = "treat")
imbalance(lg_bc1, alloc="treat")

# B 1.2.4 Urn Designs
?udPar
# r = ini
# s = add

# The function udPar creates a model object, storing the parameters for the particular Urn design. 
# To generate sequences of allocations, we use the function genSeq
ud_31 = udPar(235, 3, 1, c("0", "1"))
seq_ud31 = genSeq(ud_31)

# The allocation vector is stored as a row in the matrix M which is a slot in the genSeq object (part of the object). 
# You can access it by
seq_ud31$M[1,]
# The argument r allows you to generate r sequences, in which case the matrix M has r rows.

# B.1.3 Stratifying the dataset
# Add an ID variable so that we can keep track of the order of participants
lg_df$ID = 1:nrow(lg_df)
# split the data frame according to levels of factors
strat_gen_sm <- lg_df %>%
  group_split(preOp_gender, preOp_smoking) 
strat_gen_sm

# will create a list of data frames, one for each combination of preOp_gender and preOp_smoking. 
# In this case there are six data frames, and for example the first (accessed by strat_gen_sm[[1]] contains all participants with preOp_gender=0 and preOp_smoking=1. 
# You can choose different factors to stratify by if you want to!

# We’ll stick with the group as above. To find the numbers of participants in each group we do
group_sizes = sapply(
  1:length(strat_gen_sm),
  function(i){
    nrow(strat_gen_sm[[i]])
  }
)
group_sizes

# B.1.3.1 Another measure of imbalance
# To see how well balanced our allocations are in terms of each covariate, we will define a function summing the marginal imbalance
marg_imbalance = function(
    df,  # participant data frame, including allocation and all factor variables
    alloc, # name of allocation column
    factors # names of prognostic factors to be included
){
  df = as.data.frame(df) # deals with tibbles
  n_fact = length(factors) # the numbers of factors
  imb_sum=0                # a running total of imbalance
  for (i in 1:n_fact){     # loop through the factors 
    ind_i = (1:ncol(df))[names(df)==factors[i]]
    col_i = as.factor(df[ ,ind_i])
    levels_i = levels(col_i)
    nlevels_i = length(levels_i)
    for (j in 1:nlevels_i){ # loop through the levels of factor i
      # df_ij contains just those entries with level j of factor i
      df_ij = df[df[ ,ind_i]==levels_i[j] , ] 
      imb_ij = imbalance(df=df_ij, alloc=alloc) # find the imbalance for the sub-data-frame in which factor i has level j
      imb_sum = imb_sum + imb_ij
    }
  }
  imb_sum
}

# For example, to find the marginal imbalance over the gender and age factors, use
marg_imbalance(df=lg_df, alloc="treat", factors = c("preOp_gender", "age"))

# Note that the larger the total number of factor levels, the larger the marginal imbalance will be, 
# so if you’re comparing between methods, make sure you’re including all the same factors!

# Choose a couple of methods from before, and use them with your stratified dataset. 
# Use the marg_imbalance function to find the marginal imbalance. 
# Try this for just the factors you stratified by, and for other collections of factors. What do you expect to see?
  
# This command creates an empty list, which we will fill with allocation data frames as we go through
alloc_list = list()
# The loop works through the stratified data frames, applies SRS to allocate patients
# and stores them in alloc_list
for (i in 1:length(strat_gen_sm)){
  alloc_list[[i]] = srs(strat_gen_sm[[i]])
}
# bind all the data frames back together again
alloc_full= dplyr::bind_rows(alloc_list)
# re-order according to ID variable
alloc_full[order(alloc_full$ID),]

# It would be silly though to do this with SRS - why?
# Once you’ve performed the allocation, you can find the demographic tables, 
# imbalance and marginal imbalance as before (with alloc_full, or whatever your is called, as the data frame)

# Minimisation

# If we want to try to achieve balance for all prognostic factors, minimisation is a more suitable method. 
# The function Minirand in the package Minirand implements the minimisation algorithm.

## Information about the treatment
ntrt <- 3 # There will three treatment groups
trtseq <- c(1, 2, 3) # the treatment groups are indexed 1, 2, 3
ratio <- c(2, 2, 1)  # the treatment groups will be allocated in a 2:2:1 ratio

## The next few rows generate the participant data frame
nsample <- 120 # we will have 120 participants
c1 <- sample(seq(1, 0), nsample, replace = TRUE, prob = c(0.4, 0.6)) 
c2 <- sample(seq(1, 0), nsample, replace = TRUE, prob = c(0.3, 0.7))
c3 <- sample(c(2, 1, 0), nsample, replace = TRUE, prob = c(0.33, 0.2, 0.5)) 
c4 <- sample(seq(1, 0), nsample, replace = TRUE, prob = c(0.33, 0.67)) 
covmat <- cbind(c1, c2, c3, c4) # generate the matrix of covariate factors for the subjects
# label of the covariates 
colnames(covmat) = c("Gender", "Age", "Hypertension", "Use of Antibiotics") 
covwt <- c(1/4, 1/4, 1/4, 1/4) # equal weights/importance applied to each factor

## Applying the algorithm - start here if you already have participant data!

res <- rep(NA, nsample) # Generate a vector to store the results (the allocations)
# generate treatment assignment for the 1st subject
res[1] = sample(trtseq, 1, replace = TRUE, prob = ratio/sum(ratio)) 
# work through the remaining patients sequentially
for (j in 2:nsample)
{
  # get treatment assignment sequentially for all subjects
  # The vector res is updated and so all previous allocations are accounted for
  # covmat is the data frame of participant data
  res[j] <- Minirand(
    covmat=covmat, j, covwt=covwt, ratio=ratio, ntrt=ntrt, trtseq=trtseq, method="Range", result=res, p = 0.9
  )
}
## Store the allocation vector 'res' as 'trt1'
trt1 <- res

# Display the number of randomized subjects at covariate factors
balance1 <- randbalance(trt1, covmat, ntrt, trtseq) 
balance1

# Calculate the total imbalance of the allocation
totimbal(trt = trt1, covmat = covmat, covwt = covwt, 
         ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range")

# Adapt the code above to apply the minimisation algorithm to the licorice gargle data. Investigate how balanced the data set is.
nsample = nrow(lg_df)
res = rep(NA, nsample)
res[1] = sample(c(0,1), 1, replace = TRUE, prob = c(0.5,0.5)) 
# work through the remaining patients sequentially
for (j in 2:nsample){
  # get treatment assignment sequentially for all subjects
  # The vector res is updated and so all previous allocations are accounted for
  # covmat is the data frame of participant data - including only the covariates
  res[j] <- Minirand(
    covmat=lg_df[ ,1:7], j, covwt=rep(1,7)/7, ratio=c(1,1), ntrt=2, trtseq=c(0,1), method="Range", result=res, p = 0.9
  )
}

lg_df$treat = res

## You can now investigate the balance of the design as usual

# Display the number of randomized subjects at covariate factors
balance1 <- randbalance(lg_df$treat, covmat=lg_df[ ,1:7], ntrt=2, trtseq=c(0,1)) 
balance1

# Calculate the total imbalance of the allocation
totimbal(trt = lg_df$treat, covmat = lg_df[ ,1:7], covwt = rep(1,7)/7, 
         ratio = c(1,1), ntrt = 2, trtseq = c(0,1), method = "Range")

# Analysis

# Polyps Data
?polyps
View(polyps)
data(polyps, package = "medicaldata")

# plot baseline variable
ggplot(data=polyps, aes(x=baseline)) + geom_histogram()

# plot outcome variable
ggplot(data=polyps, aes(x=number12m)) + geom_histogram()

#These are both very right skewed, and exclusively positive (which makes sense given they are counts). 
# A sensible thing to do therefore would be to take the log of these two variables. 
# Since they are counts we will use the log to base 10, so that our results are easier to interpret. 
# If you’d prefer to use natural logs that’s fine, some of your numbers will be different.
polyps$log_baseline = log10(polyps$baseline)
polyps$log_number12m = log10(polyps$number12m)

# Now test the hypothesis that the Sulindac has had some effect.
polyps_df = na.omit(polyps) # two participants (001 and 018) don't have data for 12m, so remove these

# Shortcut for t-test of sulindac vs. placebo
t.test(
  x=polyps_df$log_number12m[polyps_df$treatment == "sulindac"],
  y=polyps_df$log_number12m[polyps_df$treatment == "placebo"],
  alternative = "two.sided",
  var.equal=T, # this makes the method use pooled variances, as we did in lectures
  conf.level = 0.95 # note that this is 1-alpha
)

# t-test, this time using the difference between outcome and baseline. 
# Hint: because we’ve taken logs, you’ll have to think about what variable to use and perhaps experiment with some possibilities

# The first step is to calculate a difference column. This is somewhat complicated by that fact that we have taken logs of the measurements. 
# Taking logs of the difference would not work, since some are negative. Potentially it might work to just work with the (unlogged) differences
polyps_df$diff = polyps_df$number12m - polyps_df$baseline
ggplot(data=polyps_df, aes(x=diff, fill=treatment)) + geom_histogram(position="dodge", bins=10)

# But the outliers look potentially problematic, and the central bulk of each distribution doesn’t look very normal. We could also try the difference of the logged measurements:
polyps_df$diff_log = polyps_df$log_number12m - polyps_df$log_baseline
ggplot(data=polyps_df, aes(x=diff_log, fill=treatment)) + geom_histogram(position="dodge", bins=10)

# This looks a lot better - no more outliers and closer to normal-looking. 
# Obviously with so few observations we won’t have a nice normal curve.
# To do a t-test, we can use R’s in-built function:
t.test(
  x=polyps_df$diff_log[polyps_df$treatment == "sulindac"],
  y=polyps_df$diff_log[polyps_df$treatment == "placebo"],
  alternative = "two.sided",
  var.equal=T, # this makes the method use pooled variances, as we did in lectures
  conf.level = 0.95 # note that this is 1-alpha
)

# In the ANCOVA model we saw in lectures so far, we fit a linear model with the outcome as dependent / target variable, 
# and the trial group and baseline as independent / explanatory variables.
lm_polyp1 = lm(log_number12m ~ treatment + log_baseline, data=polyps_df)
summary(lm_polyp1)

# Based on this, the effect of the drug sulindac is significant (p=0.00595, lower than our previous models). 
# The 95% CI for the treatment effect (which is still log(X_T / X_C now is
c(-0.7046 - qt(0.975, df=17)*0.1675, -0.7046 + qt(0.975, df=17)*0.1675 ) 

# We could create a nicer-looking table of the coefficients using tbl_regression (from gtsummary)
# tbl_regression(polyps_df) - alternative method which doesn't work...

# Unlike in the t-test where we can only compare measurements like-for-like, in ANCOVA we fit a coefficient to the baseline covariate. 
# This means we are no longer limited to comparing the outcome variable to variables on the same scale, but can also include other baseline variables.

# Fit another linear model, this time including the other baseline variables.
lm_polyp2 = lm(log_number12m ~ treatment + log_baseline + sex + age, data = polyps_df)
summary(lm_polyp2)
# In this case it appears that the other two baseline covariates, age and sex, don’t have a significant effect on the outcome.

# Inspect some plots of the residuals, to check whether these models have any systematic problems. Does there appear to be any multicollinearity in the data?
# These solutions will demonstrate some methods with the first model. First of all, we can add columns with the residuals and fitted values from the model.
polyps_df$resid1 = resid(lm_polyp1)
polyps_df$fitted1 = fitted(lm_polyp1)
ggplot(data=polyps_df, aes(x=fitted1, y=resid1, col=treatment)) + geom_point()

ggplot(data=polyps_df, aes(x=log_baseline, y=resid1, col=treatment)) + geom_point()
ggplot(data=polyps_df, aes(x=resid1, fill=treatment)) + geom_histogram(bins=20)

# None of these ring huge alarm bells, but because the dataset is so small it’s quite hard to tell! Arguably the residuals for the sulindac group are more spread out than those for the plaecbo, but there are no obvious systematic trends.
# We can find the variance inflation factor by
vif(lm_polyp1)

# and it indicates no problems. Another sign that the ANCOVA model is well-fitted is that the standard errors of the coefficients are all reasonably small. If there was a problem, for example multicollinearity, some of these would blow up.








