# SUPPORTING FUNCTIONS FOR seq.bf()

# ------------------------------------------------------------------------------
# BACKGROUND FUNCTIONS
# these are from the Bayes Factor Package
# ------------------------------------------------------------------------------

# from the BayesFactor package - have to include in script for now because 
# it is used in modified functions below
marshallTibble <- function(data) {
  if (inherits(data, 'tbl_df')) {
    data <- as.data.frame(data)
    warning('data coerced from tibble to data frame', call.=FALSE)
  }
  data
}


checkCallback = function(callback, ... ){
  ret = as.integer(callback(...))
  if(ret) stop("Operation cancelled by callback function. ", ret)
  return(ret)
}

makeTtestHypothesisNames = function(rscale, nullInterval=NULL, mu = 0){
  if(is.null(nullInterval)){
    shortName = paste("Alt., r=",rscale,sep="")
    longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, sep="")
  }else{
    if(!is.null(attr(nullInterval,"complement"))){
      shortName = paste("Alt., r=",rscale," !(",nullInterval[1],"<d<",nullInterval[2], ")",sep="")
      longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, " !(",nullInterval[1],"<d<",nullInterval[2],")",sep="")
    }else{
      shortName = paste("Alt., r=",rscale," ",nullInterval[1],"<d<",nullInterval[2],sep="")
      longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, " ",nullInterval[1],"<d<",nullInterval[2],sep="")
    }
  }
  return(list(shortName=shortName,longName=longName))
}

BFoneSample <- function(type, identifier, prior, shortName, longName, analysis = list()){
  new("BFoneSample", type = type,
      identifier = identifier,
      prior = prior,
      shortName = shortName,
      longName = longName,
      analysis = analysis,
      version = BFInfo(FALSE))
}

BFindepSample <- function(type, identifier, prior, shortName, longName, analysis = list()){
  new("BFindepSample", type = type,
      identifier = identifier,
      prior = prior,
      shortName = shortName,
      longName = longName,
      analysis = analysis,
      version = BFInfo(FALSE))
}

# ------------------------------------------------------------------------------
# A little desconstructed code from the Bayes Factor package
# The main reason we have to deconstruct this is so we can set our own priors 
# for the cauchy distribution in the Bayes Factor functions

# -------------------------------------
# changing original ttestBF_oneSample()
# commenting out line that forces rscale to be categorical

mod_ttestBF_oneSample = function(x, mu, nullInterval, rscale, posterior, callback, ... ){
  # rscale = rpriorValues("ttestOne",,rscale)
  hypNames = makeTtestHypothesisNames(rscale, nullInterval, mu = mu)
  
  mod1 = BFoneSample(type = "JZS", 
                     identifier = list(formula = "y ~ 1", nullInterval = nullInterval), 
                     prior=list(rscale=rscale, mu=mu, nullInterval = nullInterval),
                     shortName = hypNames$shortName,
                     longName = hypNames$longName)
  
  if(posterior)
    return(posterior(mod1, data = data.frame(y=x), callback = callback, ...))
  
  bf1 = compare(numerator = mod1, data = data.frame(y=x))
  
  if(!is.null(nullInterval)){
    mod2 = mod1
    attr(mod2@identifier$nullInterval, "complement") = TRUE
    attr(mod2@prior$nullInterval, "complement") = TRUE
    hypNames = makeTtestHypothesisNames(rscale, mod2@identifier$nullInterval, mu = mu)
    mod2@shortName = hypNames$shortName
    mod2@longName = hypNames$longName
    
    bf2 = compare(numerator = mod2, data = data.frame(y=x))
    checkCallback(callback,as.integer(1000))
    return(c(bf1,bf2))
  }else{
    checkCallback(callback,as.integer(1000))
    return(bf1)
  }    
  
} 



# -------------------------------------
# changing original ttestBF_indepSample()
# just commenting out the line that forces rscale to be a categorical input


mod_ttestBF_indepSample = function(formula, data, mu, nullInterval, rscale, posterior, callback, ... ){
  checkFormula(formula, data, analysis = "indept")
  
  #rscale = rpriorValues("ttestTwo",,rscale)
  hypNames = makeTtestHypothesisNames(rscale, nullInterval, mu = mu)
  
  mod1 = BFindepSample(type = "JZS", 
                       identifier = list(formula = stringFromFormula(formula), nullInterval = nullInterval), 
                       prior=list(rscale=rscale, mu=mu, nullInterval = nullInterval),
                       shortName = hypNames$shortName,
                       longName = hypNames$longName
  )
  
  if(posterior)
    return(posterior(mod1, data = data, callback = callback, ...))
  
  bf1 = compare(numerator = mod1, data = data)
  
  if(!is.null(nullInterval)){
    mod2 = mod1
    attr(mod2@identifier$nullInterval, "complement") = TRUE
    attr(mod2@prior$nullInterval, "complement") = TRUE
    hypNames = makeTtestHypothesisNames(rscale, mod2@identifier$nullInterval, mu = mu)
    mod2@shortName = hypNames$shortName
    mod2@longName = hypNames$longName
    
    bf2 = compare(numerator = mod2, data = data)
    checkCallback(callback,as.integer(1000))
    return(c(bf1, bf2))
  }else{
    checkCallback(callback,as.integer(1000))
    return(bf1)
  }  
  
  
} 



# -------------------------------------

# the original ttestBF uses two inner functions that we need to correct to allow for 
# a number rscale parameter (sets the cauchy prior distribution):
# 
# ttestBF_oneSample()
# ttestBF_indepSample()
# 
# This function has been modified to include the modified versions of the above 
# inner functions
mod_ttestBF <- function(x = NULL, y = NULL, formula = NULL, mu = 0, nullInterval = NULL,
                        paired = FALSE, data = NULL, rscale=0.5, posterior=FALSE, callback = function(...) as.integer(0), ...){
  
  data <- marshallTibble(data)
  
  if(!is.null(x) & !is.null(formula)) stop("Only one of x or formula should be defined.")
  
  if(!is.null(x) | !is.null(y))
    if(any(is.na(c(x,y))) | any(is.infinite(c(x,y))))
      stop("x or y must not contain missing or infinite values.")
  
  if(!is.null(nullInterval)){
    nullInterval = range(nullInterval)
    if(identical(nullInterval,c(-Inf,Inf))){
      nullInterval = NULL
    }
  }
  
  checkCallback(callback,as.integer(0))
  
  if( (is.null(formula) & is.null(y)) | (!is.null(y) & paired) ){ # one sample
    if(paired){
      # check that the two vectors have same length
      if(length(x)!=length(y)) stop("Length of x and y must be the same if paired=TRUE.")
      x = x - y
    }
    return( mod_ttestBF_oneSample(x = x, mu = mu,
                                  nullInterval = nullInterval,
                                  rscale = rscale, posterior = posterior,
                                  callback = callback, ... ) )
  }
  if(!is.null(y) & !paired){ # Two-sample; create formula
    if(!is.null(data) | !is.null(formula)) stop("Do not specify formula or data if x and y are specified.")
    data = data.frame(y = c(x,y),
                      group = factor(c(rep("x",length(x)),rep("y",length(y))))
    )
    formula = y ~ group
  }
  if(!is.null(formula)){ # Two-sample
    if(paired) stop("Cannot use 'paired' with formula.")
    if(is.null(data)) stop("'data' needed for formula.")
    if(mu != 0) stop("Use of nonzero null hypothesis not implemented for independent samples test.")
    return(mod_ttestBF_indepSample(formula = formula, data = data, mu = mu,
                                   nullInterval = nullInterval, rscale = rscale,
                                   posterior = posterior, callback = callback, ... ))
  }else{
    stop("Insufficient arguments to perform t test.")
  }
}


checkFormula <- function(formula, data, analysis){
  if(length(formula) < 3) stop("LHS of formula must be given.")
  cnames = colnames(data)
  
  dv = stringFromFormula(formula[[2]])
  
  if(!is.numeric(data[,dv])) stop("Dependent variable must be numeric.")
  if(any(is.na(data[,dv])) | any(is.infinite(data[,dv]))) stop("Dependent variable must not contain missing or infinite values.")
  
  factors = fmlaFactors(formula, data)
  terms = colnames(attr(terms(formula, data = data),"factors"))
  decom = decomposeTerms(terms)
  terms = unlist(decom)
  
  vars = rownames(attr(terms(formula, data = data),"factors"))
  vars = unlist(decomposeTerms(vars))
  
  if(any(is.na(data[,vars]))) stop("Predictors must not contain missing values.")
  
  if(is.null(factors)) return()
  if(factors[1] %in% terms) stop("Dependent variable cannot be a predictor.")
  if(!all(factors %in% cnames)) stop("Some variables missing in data frame.")
  
  if(analysis=="regression"){
    RHS = stringFromFormula(formula[[3]])
    lengths = sapply(decom, length)
    if (any(lengths > 1)) stop("Interactions not allowed in regressionBF (try generalTestBF).")
  }
  
  if(analysis=="lm" | analysis=="anova" | analysis == "regression" | analysis == "indept")
    if(attr(terms(formula, data = data),"intercept") == 0) stop("Formula must include intercept.")
  
  if(analysis=="indept"){
    if( length(decom) > 1 ) stop("Indep. groups t test can only support 1 factor as predictor.")
    if(length(decom[[1]]) > 1) stop("Interaction terms are not allowed in t test.")
    if(nlevels(factor(data[,terms])) > 2) stop("Indep. groups t test requires a factor with exactly 2 levels.")
  }
  invisible()
}


stringFromFormula <- function(formula){
  oneLine = paste(deparse(formula),collapse="")
  sub("\\s\\s+"," ", oneLine, perl=TRUE) # get rid of extra spaces
}


fmlaFactors <- function(formula, data){
  names <- rownames(attr(terms(formula, data = data),"factors"))
  names <- decomposeTerms(names)
  names <- unlist(names)
  names
}



decomposeTerm <- function(term) {
  
  chars <- strsplit(term, '')[[1]]
  components <- character()
  componentChars <- character()
  inQuote <- FALSE
  
  i <- 1
  n <- length(chars)
  
  while (i <= n) {
    char <- chars[i]
    if (char == '`') {
      inQuote <- ! inQuote
    }
    else if (char == '\\') {
      i <- i + 1
      char <- chars[i]
      componentChars <- c(componentChars, char)
    }
    else if (char == ':' && inQuote == FALSE) {
      component <- paste0(componentChars, collapse='')
      components <- c(components, component)
      componentChars <- character()
    }
    else {
      componentChars <- c(componentChars, char)
    }
    i <- i + 1
  }
  
  component <- paste0(componentChars, collapse='')
  components <- c(components, component)
  
  components
}


decomposeTerms <- function(terms) {
  decomposed <- list()
  for (i in seq_along(terms))
    decomposed[[i]] <- decomposeTerm(terms[[i]])
  decomposed
}
