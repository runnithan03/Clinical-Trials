## Extension Problems

# A simulation experiment!

# The general process for this is:

# Choose a numerical summary (probably either the imbalance or the marginal imbalance)
# Choose some large number  n_sim and create a vector of zeroes or NAs that length

# Run through the following process n_sim times:
  # Perform the allocation
  # Find the numerical summary
  # Store the numerical summary in the vector you created in step 2.
# Plot (or otherwise summarise) the vector of summaries.

# Perform the simulation experiment for some of the methods above. 
# Make sure you store the summary vectors by different (and intelligible!) names.
# Based on your results, comment on how plausible it is that Ruetzler et al. (2013) used each method in their trial.

# Hint: think of the simulated summaries as approximations of probability distributions

## TBD this exercise...




# Treatment for maternal periodontal disease
require(medicaldata)
?opt
opt

ggplot(data=opt, aes(x=Birthweight, fill=Group)) + geom_histogram(position="dodge")

# data looks very nice and normal (although there is a bit of a fat tail to the left), so we won’t transform it.

opt_red = opt[ ,c(1:22,72)] 
# Change NAs to "None" for diabetic
# since in opt the non-diabetic people are coded as NA (and therefore excluded from the model)
diab = as.character(opt_red$BL.Diab.Type)
diab[is.na(diab)] = "None"
opt_red$BL.Diab.Type = as.factor(diab)
# similar problem with smokers and how many cigarettes per day

# If people are non-smokers and have missing for number of cigarettes per day
# change their number of cigarettes to zero
sm = opt_red$Use.Tob
cigs = opt_red$BL.Cig.Day
cigs[(is.na(cigs)&(sm=="No "))] = 0
opt_red$BL.Cig.Day = cigs

# Same for alcohol and drinks per day

alc = opt_red$Use.Alc
dr = opt_red$BL.Drks.Day
dr[(is.na(dr)&(alc=="No "))] = 0
opt_red$BL.Drks.Day = dr

# If a participant hasn't had a previous pregnancy, her N.prev.preg should be zero (not NA)

pp = opt_red$Prev.preg
npp = opt_red$N.prev.preg
npp[pp=="No "] = 0
opt_red$N.prev.preg = npp

opt_red

# When we use lm to fit an ANCOVA model, all rows with an NA in will be ignored. 
# It’s therefore important to try to eradicate any NAs where we actually do know the value!

# Perform a  t-test on the outcome Birthweight. 
# To perform the t-test we can use the inbuilt R function

# Check SDs are fairly close before proceeding
sd(opt_red$Birthweight[opt_red$Group == "T"], na.rm=T)
sd(opt_red$Birthweight[opt_red$Group == "C"], na.rm=T)
t.test(
  x=opt_red$Birthweight[opt_red$Group == "T"],
  y=opt_red$Birthweight[opt_red$Group == "C"],
  alternative = "two.sided",
  var.equal=T, # this makes the method use pooled variances, as we did in lectures
  conf.level = 0.95 # note that this is 1-alpha
)

# We find that the difference in Birthweights between the groups is not even close to significant (p=0.456).
# We can’t do a t-test with differences because there is no comparable baseline measurement.

# Now we will move on to ANCOVA.
# To fit a linear model including every variable as a covariate (apart from the target variable), we can do
lm_full = lm(Birthweight ~ ., data = opt_red[ ,-1]) # don't include the ID column!
summary(lm_full)

# We see from the model summary that our model is terrible:  
# R^2 =  0.03241, which means it is explaining about 3% of the variance in Birthweight.
# If you want to use only certain terms, you can include them in the formula, for example
lm_eg = lm(Birthweight ~ Group + Age + Hypertension, data=opt_red)

# We see that, as we expected, the Group variable is not significant (p=0.5644). 
# However, some terms are significant, for example whether or not the participant has diabetes, 
# and how many cigarettes a mother smokes per day - this isn’t surprising given the contextual information we had.

# Perform some diagnostic checks on your model. Do you have any reason to suspect it isn’t adequate?
# The first step is to extract the residuals and the fitted values (which are also useful). 
# We will create a new data frame called opt_diag with these in, so that we can plot things easily but don’t pollute our original dataset (in case we want to fit any more models)
opt_diag = na.omit(opt_red) # lm only fits where all variables are present
opt_diag$resid = resid(lm_full)
opt_diag$fitted = fitted(lm_full)

opt_diag

# plot examples:
ggplot(data=opt_diag, aes(x=resid, fill=Group)) + geom_histogram(position="dodge")
ggplot(data=opt_diag, aes(x=fitted, y=resid, col=Group)) + geom_point()

# Given the results of your ANCOVA model, what do you think the risks would be if the study had been much smaller, or if the allocation had not been well-balanced?
# In this case, it would be possible for the baseline factors that did turn out to be significant to make it appear that the treatment had a significant effect. 
# This wouldn’t be possible with ANCOVA, but it would with a t-test in which the other covariates aren’t considered.
