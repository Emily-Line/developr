# Functions for developr
# Obviously not in package form yet

library(BayesFactor)


# Run script with supporting functions before running this. 

# ------------------------------------------------------------------------------
# FUNCTION: seq.BF()
# ------------------------------------------------------------------------------

# Threshold for Bayes factor is large. Many studies that replicated didn't have BFs this large. 

# calculates Bayes Factor for t-test, anova, or correlation
# Tells user if BF has reached their chosen threshold yet

# TODO: Discuss what should happen with the anova case since there are multiple 
# bayes factors. What should the messages be? 
# TODO: add null hypothesis supported case for anova
# TODO: Add message for all results: BF less than [threshold] do not constitute evidence for the alternative
#       We recommend... 

seq.BF <- function(
  testtype = c("t-test", "anova", "correlation"), 
  x = NULL, 
  y = NULL,
  formula = NULL,
  data = NULL,
  rscale = .5,
  rscaleRandom = "nuisance",
  threshold = 7,
  paired = FALSE,
  mu = NULL,
  whichRandom,
  ...){
  
  if(missing(testtype)){
    stop("Missing argument: Please indicate a test type of ``t-test'', ``anova'', or ``correlation''")
  }
  
  if(testtype == "regression"){
    stop("Sorry, regression is not yet supported with this function.")
  }
  
  if(testtype != "t-test" & testtype != "anova" & testtype != "correlation") {
    stop("Incorrect argument:Please indicate a test type of ``t-test'', ``anova'', or ``correlation''.
         Note that this is case sensitive.")
  }
  
  if(threshold < 1 | !is.numeric(threshold)){
    stop("Please put a number greater than 1 for the threshold parameter")
  }
  
  if(threshold > 20){
    warning("Threshold for Bayes factor is larger than usual. It will take more evidence to exceed this threshold.")
  }
  
  if(threshold < 4){
    warning("Threshold for Bayes factor is smaller than usual. This threshold will not provide strong evidence for the alternative or the null.")
  }
  
  if(rscale > sqrt(2)){
    warning("The rscale parameter is much larger than usual. It is unlikely that you need a prior distribution this wide.")
  }

  
  if(testtype == "t-test"){ # t-test
    BFobj = mod_ttestBF(x = x, y = y, formula = formula, data = data, paired = paired, rscale = rscale)
    
    # if(hasArg(data)) {
    #   BFobj = mod_ttestBF(x = x, y = y, paired = paired, rscale = rscale, data = data)
    # } else {
    #   BFobj = mod_ttestBF(x = x, y = y, formula = formula, paired = paired, rscale = rscale)
    # }
    
    BF = extractBF(BFobj)$bf
    
  } else if(testtype == "anova"){ # anova
    
    BFobj = anovaBF(formula = formula, data = data, rscaleFixed = rscale, rscaleRandom = rscaleRandom, whichRandom = whichRandom)
    BF_df <- extractBF(BFobj)
    BF_vals <- BF_df$bf
    
  } else{ # correlation
    
    BFobj = correlationBF(y=y, x=x, rscale=rscale)
    BF = extractBF(BFobj)$bf
    
  }
  
  bfexceedalt <- paste("Bayes factor has exceeded the threshold of ", threshold, "! Alternative hypothesis is supported.", sep = "")
  bfexceednull <- paste("Bayes factor has exceeded the null threshold of ", round(1/threshold, 2), "! Null hypothesis is supported.", sep = "")
  bfnotexceed <- paste("Bayes factor has not exceeded the threshold of ", threshold, ". " , sep = "")
  end_message <- paste("A Bayes factor less than ", threshold, "or greater than ", round(1/threshold, 2), " does not constitute evidence for the alternative or null hypothesis, respectively. 
  \n We recommend you keep sampling if you want to make inferences about an effect that has a Bayes factor \n
  that does not exceed the threshold for the alternative or the null hypothesis.")
  end_message1 <- "A Bayes factor less than "
  end_message1a <- ". A Bayes factor less than "
  end_message2 <- " does not constitute evidence for the alternative hypothesis. \n 
  We recommend you keep sampling if you want to make inferences about a model that has a Bayes factor \n
  that has not exceeded the threshold for the alternative or the null hypothesis."
  
  # bfexceed_anova <- "All Bayes factors have exceeded the threshold of "
  # bfnotexceed_anova <- "None of the Bayes factors have exceeded the threshold of "
  # bfmixed_anova <- "Some Bayes factors have exceeded the threshold "
  # bfmixed_end <- ". If there are effects you are still interested in that have not exceeded
  # the threshold, keep sampling."
  
  anova_message <- "ANOVA models usually have more than one Bayes Factor. Each Bayes factor corresponds to an \n effect in your model. "
  
  
  if(testtype == "anova"){
    threshold_message <- paste(anova_message, end_message, sep = "\n")
    # BF_exceed <- NA
    # 
    # for(i in 1:length(BF_vals)) {
    #   if(BF_vals[i] > threshold) {
    #     BF_exceed[i] <- 1
    #   } else { 
    #     BF_exceed[i] <- 0
    #     }
    # }
    # 
    # BF_mean <- mean(BF_exceed)
    # 
    # if(BF_mean == 1) {
    #   threshold_message <- paste(bfexceed_anova, threshold, alternative, sep = "")
    # } else if(BF_mean == 0) {
    #   threshold_message <- paste(bfnotexceed_anova, threshold, keep, sep = "")
    # } else if(BF_mean < 1) {
    #   threshold_message <- paste(bfmixed_anova, threshold, bfmixed_end, sep = "")
    # }

  } else {
        if(BF > threshold){
          threshold_message <- paste(bfexceedalt, "", end_message, sep = "\n")
        } else if(BF < 1/threshold){
          threshold_message <- paste(bfexceednull, "", end_message, sep = "\n")
        } else{
          threshold_message <- paste(bfnotexceed, "", end_message, sep = "\n")
        }
  }
  

  
  message(threshold_message)
  return(BFobj)
  
}

# NOTES: 
# 
# BF Threshold: 
# Threshold for minimum BF is set to 10 right now, I have some data from link 
# below of distribution of Bayes Factors, but it's difficult to figure out if a 
# BF other than 10 would be appropriate, majority of the BFs of successfully 
# replicated studies are above BF = 20. 
# https://alexanderetz.com/2015/08/30/the-bayesian-reproducibility-project/
# 
# Cauchy parameter: 
# So apparently with BayesFactor you cannot set a continuous value for the r 
# parameter of the Cauchy distribution, it only has three factor values you can choose,
# with sqrt(2)/2 being the most narrow distribution (which we already decided was not 
# ideal). Setting to that value for now, but it's something to discuss. 


# ------------------------------------------------------------------------------
# EXAMPLES
# ------------------------------------------------------------------------------

# Using data from BayesFactor package until we figure out our own dataset

data(sleep)
sleepx = sleep$extra[sleep$group==1] 
sleepy = sleep$extra[sleep$group==2]

# t-test
# (currently only works with x, y, not with formula)
seq.BF(testtype = "t-test",
       x = sleepx, 
       y = sleepy,
       paired=TRUE, 
       rscale = .5)


seq.BF(testtype = "t-test",
       x = sleepx, 
       y = sleepy,
       data = sleep,
       rscale = .25,
       paired=TRUE)

ttestBF(x = sleepx, y = sleepy, paired=TRUE)

x = sampleA$RT
y = sampleA$shape

seq.BF(testtype = "t-test",
       x = x,
       y = y,
       paired = FALSE)

ttestBF(formula = RT ~ shape, data = puzzles)
ttestBF(x, y)


# anova
# Need to change code to allow for the multiple BFs
# Trying to think of a way to make the formula more user-friendly.. 
data(puzzles)

seq.BF(testtype = "anova",
       formula = RT ~ shape*color,
       data = puzzles,
       whichRandom = "ID",
       rscale = .5)

anovaBF(formula = RT ~ shape*color + ID, data = puzzles, rscaleFixed = 1, rscaleRandom = "nuisance", whichRandom = "ID")

# correlation

seq.BF(testtype = "correlation",
       y = iris$Sepal.Length,
       x = iris$Sepal.Width,
       rscale = .1)

seq.BF(testtype = "regression",
       y = iris$Sepal.Length,
       x = iris$Sepal.Width,
       rscale = .1)

