#######################################
## GENERATE CONFIGURATIONS
#######################################

## When called with an unconditional parameter, it
## must return TRUE
conditionsSatisfied <- function (parameters, partialConfiguration, paramName)
{
  condition <- parameters$conditions[[paramName]]
  # If there is no condition, do not waste time evaluating it.
  if (isTRUE(condition)) return(TRUE)

  v <- eval(condition, as.list(partialConfiguration))
  # Return TRUE if TRUE, FALSE if FALSE or NA
  ## FIXME: If we byte-compile the condition, then we should incorporate the
  ## following into the condition directly. See readForbiddenFile.
  v <- !is.na(v) && v
  return(v)
}

new_empty_configuration <- function(parameters)
{
  newConfigurationsColnames <- c(names(parameters$conditions), ".PARENT.")
  return(setNames(as.list(rep(NA, length(newConfigurationsColnames))),
                  newConfigurationsColnames))
}

get.fixed.value <- function(param, parameters)
{
  value <- parameters$domain[[param]][1]
  type <- parameters$types[[param]]
  if (type == "i") {
    return (as.integer(value))
  } else if (type == "c" || type == "o") {
    return (value)
  } else {
    irace.assert (type == "r")
    return (as.double(value))
  }
}

## Calculates the parameter bounds when parameters domain is dependent
getDependentBound <- function(parameters, param, configuration)
{
  values <- parameters$domain[[param]]
  if (is.expression(values)) {
    # Depends contains parameters that enable param and parameters that define
    # its domain. If this is a partial configuration, we need only the latter.
    # Use names() here in case the configuration is simply a list.
    deps <- intersect(names(configuration), parameters$depends[[param]])
    # If it depends on a parameter that is disabled, then this is disabled.
    if (anyNA(configuration[deps])) return(NA)

    values <- sapply(values, eval, configuration)
    irace.assert(all(is.finite(values)))
    # Value gets truncated (defined from robotics initial requirements)
    if (parameters$types[param] == "i") values <- as.integer(values)
    if (values[1] > values[2]) {
      irace.error ("Invalid domain (", paste0(values, collapse=", "),
                   ") generated for parameter '", param,
                   "' that depends on parameters (",
                   paste0(parameters$depends[[param]], collapse=", "),
                   "). This is NOT a bug in irace. Check the definition of these parameters.")
    }
  }

  return(values)
}

### Uniform sampling for the initial generation
sampleUniform <- function (parameters, nbConfigurations, digits,
                           forbidden = NULL, repair = NULL)
{
  if (is.null(repair)) repair <- function(c, p, d) c
  
  namesParameters <- names(parameters$conditions)
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbConfigurations,
                         ncol = length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT."))
                         ), stringsAsFactors=FALSE)
  empty_configuration <- new_empty_configuration(parameters)

  for (idxConfiguration in seq_len(nbConfigurations)) {
    forbidden.retries <- 0
    while (forbidden.retries < 100) {
      configuration <- empty_configuration
      for (p in seq_along(namesParameters)) {
        currentParameter <- namesParameters[p]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          configuration[[p]] <- NA
          next
        }
        # FIXME: We must be careful because parameters$types does not have the
        # same order as namesParameters, because we sample in the order of the
        # conditions.
        currentType <- parameters$types[[currentParameter]]
        if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled          
        } else if (currentType %in% c("i","r")) {
          domain <- getDependentBound(parameters, currentParameter, configuration)
          newVal <- sample_unif(currentType, domain,
                                transf = parameters$transform[[currentParameter]],
                                digits)
        } else {
          irace.assert(currentType %in% c("c","o"))
          possibleValues <- parameters$domain[[currentParameter]]
          newVal <- sample(possibleValues, 1)
        }
        configuration[[p]] <- newVal
      }
      configuration <- as.data.frame(configuration, stringsAsFactors=FALSE)
      configuration <- repair(configuration, parameters, digits)

      if (is.null(forbidden)
          || nrow(checkForbidden(configuration, forbidden)) == 1) {
        newConfigurations[idxConfiguration,] <- configuration
        break
      }
      forbidden.retries <- forbidden.retries + 1
    }
    if (forbidden.retries >= 100) {
      irace.error("irace tried 100 times to sample from the model a configuration not forbidden without success, perhaps your constraints are too strict?")
    }
  }
  return (newConfigurations)
}

# To be called the first time before the second race (with indexIter =
# 2) Nb configurations is the number of configurations at the end
# included the elite ones obtained from the previous iteration
sampleModel <- function (parameters, eliteConfigurations, model,
                         nbNewConfigurations, digits, forbidden = NULL,
                         repair = NULL)
{
  if (is.null(repair)) repair <- function(c, p, d) c
  
  if (nbNewConfigurations <= 0) {
    irace.error ("The number of configurations to generate appears to be negative or zero.")
  }
  namesParameters <- names(parameters$conditions)
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbNewConfigurations,
                         ncol = length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT."))
                         ), stringsAsFactors=FALSE)
  empty_configuration <- new_empty_configuration(parameters)
  
  for (idxConfiguration in seq_len(nbNewConfigurations)) {
    forbidden.retries <- 0
    while (forbidden.retries < 100) {
      # Choose the elite which will be the parent.
      indexEliteParent <- sample.int (n = nrow(eliteConfigurations), size = 1,
                                      prob = eliteConfigurations[[".WEIGHT."]])
      eliteParent <- eliteConfigurations[indexEliteParent, ]
      idEliteParent <- eliteParent[[".ID."]]
      configuration <- empty_configuration
      configuration[[".PARENT."]] <- idEliteParent
      
      # Sample a value for every parameter of the new configuration.
      for (p in seq_along(namesParameters)) {
        # FIXME: We must be careful because parameters$types does not
        # have the same order as parameters$conditions. Ideally, we
        # should fix this or make it impossible to confuse them.
        currentParameter <- namesParameters[p]
        currentType <- parameters$types[[currentParameter]]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          # Some conditions are unsatisfied.
          # Should be useless, NA is (always?) assigned when matrix created
          newVal <- NA
          
        } else if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled
        } else if (currentType %in% c("i", "r")) {
          domain <- getDependentBound(parameters, currentParameter, configuration)
          mean <- as.numeric(eliteParent[currentParameter])
          # If there is not value we obtain it from the model or the mean obtained is
          # not in the current domain, this can happen when using dependent domains
          if (is.na(mean)) mean <- model[[currentParameter]][[as.character(idEliteParent)]][2]
          if (is.na(mean) || !inNumericDomain(mean, domain)) {
            # The elite parent does not have any value for this parameter,
            # let's sample uniformly.
            newVal <- sample_unif(currentType, domain,
                                  transf = parameters$transform[[currentParameter]], digits)
          } else {
            stdDev <- model[[currentParameter]][[as.character(idEliteParent)]][1]
            # If parameters are dependent standard deviation must be computed
            # based on the current domain
            if (parameters$isDependent[currentParameter]) {
              # Conditions should be satisfied for the parameter, thus domain cannot be NA
              stdDev <- (domain[2] - domain[1]) * stdDev
            }
            newVal <- sample_norm(mean, stdDev, currentType, domain,
                                  transf = parameters$transform[[currentParameter]],
                                  digits)
          }
        } else if (currentType == "o") {
          possibleValues <- paramDomain(currentParameter, parameters)  
          value <- eliteParent[currentParameter]
          
          if (is.na(value)) {
            # The elite parent does not have any value for this
            # parameter, let's sample uniformly
            ## FIXME: We should save the last used parameter in the model and use it here.
            newVal <- sample(possibleValues, 1)
          } else {
            # Find the position within the vector of possible
            # values to determine the equivalent integer.
            mean <- match(value, possibleValues) # Return index of value in array
            stdDev <- model[[currentParameter]][[as.character(idEliteParent)]]

            # Sample with truncated normal distribution as an integer.
            # See sample_norm() for an explanation.
            newValAsInt <- floor(rtnorm(1, mean + 0.5, stdDev, lower = 1,
                                        upper = length(possibleValues) + 1L))

            # The probability of this happening is very small, but it can happen.
            if (newValAsInt == length(possibleValues) + 1L)
              newValAsInt <- length(possibleValues)
            
            irace.assert(newValAsInt >= 1L && newValAsInt <= length(possibleValues))
            # Get back to categorical values, find the one corresponding to the
            # newVal
            newVal <- possibleValues[newValAsInt]
          } 
        } else if (currentType == "c") {
          # FIXME: Why is idEliteParent character?
          # FIXME: Why the model is <parameter><Parent>? It makes more sense to be <Parent><parameter>.
          probVector <- model[[currentParameter]][[as.character(idEliteParent)]]
          possibleValues <- paramDomain(currentParameter, parameters)
          newVal <- sample(x = possibleValues, size = 1, prob = probVector)
        } else {
          irace.internal.error("Unexpected condition in sampleModel")
        }
        configuration[[p]] <- newVal
      }
      
      configuration <- as.data.frame(configuration, stringsAsFactors = FALSE)
      configuration <- repair(configuration, parameters, digits)

      if (is.null(forbidden)
          || nrow(checkForbidden(configuration, forbidden)) == 1) {
        newConfigurations[idxConfiguration,] <- configuration
        break
      }
      forbidden.retries <- forbidden.retries + 1
    }
    if (forbidden.retries >= 100) {
      irace.error("irace tried 100 times to sample from the model a configuration not forbidden without success, perhaps your constraints are too strict?")
    }
  }
  return (newConfigurations)
}

transform.from.log <- function(x, transf, lowerBound, upperBound)
{
  trLower <- attr(transf, "lower") 
  trUpper <- attr(transf, "upper")
  x <- exp(trLower + (trUpper - trLower) * x)
  return(x)
}

transform.to.log <- function(x, transf, lowerBound, upperBound)
{
  trLower <- attr(transf, "lower") 
  trUpper <- attr(transf, "upper")
  return((log(x) - trLower)/(trUpper - trLower))
}
## How to sample integer values?
#
# The problem: If we have an integer with domain [1,3] and we sample a real value
# and round, then there are more chances of getting 2 than 1 or 3:
# [1, 1,5) -> 1
# [1.5, 2,5) -> 2
# [2.5, 3) -> 3
#
# The solution: Sample in [lowerbound, upperbound + 1], that is, [1, 4], then floor():
# [1, 2) -> 1
# [2, 3) -> 2
# [3, 4) -> 3
#
# Why floor() and not trunc()?
# Because trunc(-1.5) -> -1, while floor(-1.5) -> -2, so for a domain [-3,-1]:
#
# [-3, -2) -> -3
# [-2, -1) -> -2
# [-1, 0)  -> -1
#
# Issue 1: We can sample 4 (upperbound + 1). In that case, we return 3.
#
# Issue 2: When sampling from a truncated normal distribution, the extremes are
# not symmetric.
#
# nsamples <- 100000
# table(floor(rtnorm(nsamples, mean=1, sd=1, lower=1,upper=4)))/nsamples
# table(floor(rtnorm(nsamples, mean=3, sd=1, lower=1,upper=4)))/nsamples
#
# To make them symmetric, we translate by 0.5, so that the mean is at the
# actual center of the interval that will produce the same value after
# truncation, e.g., given an integer value of 1, then mean=1.5, which is at the
# center of [1,2).
#
# nsamples <- 100000
# table(floor(rtnorm(nsamples, mean=1.5, sd=1, lower=1,upper=4)))/nsamples
# table(floor(rtnorm(nsamples, mean=3.5, sd=1, lower=1,upper=4)))/nsamples
#
# The above reasoning also works for log-transformed domains, because 
# floor() happens in the original domain, not in the log-transformed one,
# except for the case of log-transformed negative domains, where we have to
# translate by -0.5.
# 
numeric_value_round <- function(type, value, lowerBound, upperBound, digits)
{  
  irace.assert(is.finite(value))
  if (type == "i") {
    value <- floor(value)
    # The probability of this happening is very small, but it could happen.
    if (value == upperBound + 1L)
      value <- upperBound
  } else
    value <- round(value, digits)

  irace.assert(value >= lowerBound && value <= upperBound)
  return (value)
}

# Sample value for a numerical parameter.
sample_unif <- function(type, domain, transf, digits)
{
  # Dependent domains could be not available because of inactivity of parameters
  # on which they are depedent. In this case, the dependent parameter becomes 
  # not active and we return NA.
  if (anyNA(domain)) return(NA)

  lowerBound <- domain[1]
  upperBound <- domain[2]

  if (type == "i") {
    # +1 for correct rounding before floor()
    upperBound <- 1L + upperBound
  }
  if (transf == "log") {
    value <- runif(1, min = 0, max = 1)
    value <- transform.from.log(value, transf, lowerBound, upperBound)
  } else {
    value <- runif(1, min = lowerBound, max = upperBound)    
  }
  # We use original upperBound, not the +1L for 'i'.
  value <- numeric_value_round(type, value, lowerBound, upperBound = domain[2], digits)
  return(value)
}

sample_norm <- function(mean, sd, type, domain, transf, digits)
{
  # Dependent domains could be not available because of inactivity of parameters
  # on which they are depedent. In this case, the dependent parameter becomes 
  # not active and we return NA.
  if (anyNA(domain)) return(NA)

  lowerBound <- domain[1]
  upperBound <- domain[2]

  if (type == "i") {
    upperBound <- 1L + upperBound
    # Because negative domains are log-transformed to positive domains.
    mean <- mean + 0.5
  }
  
  if (transf == "log") {
    trMean <- transform.to.log(mean, transf, lowerBound, upperBound)
    value <- rtnorm(1, trMean, sd, lower = 0, upper = 1)
    value <- transform.from.log(value, transf, lowerBound, upperBound)
  } else {
    value <- rtnorm(1, mean, sd, lowerBound, upperBound)
  }
  # We use original upperBound, not the +1L for 'i'.
  value <- numeric_value_round(type, value, lowerBound, upperBound = domain[2], digits)
  return(value)
}
