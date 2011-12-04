###############################################################
# In addition to the Independent Normal class, we also have a
# class for independent log normal, because the true expression
# value follows a log normal distribution in our model.

##=============================================================================
setClass("IndependentLogNormal",
         representation(logmu="numeric",
                        logsigma="numeric"))


##-----------------------------------------------------------------------------
## Generates an IndependentLogNormal object.
IndependentLogNormal <- function(logmu, logsigma) {
  if (length(logsigma) != 1 & length(logsigma) != length(logmu))
    stop("'logmu' and 'logsigma' have different lengths")
  new("IndependentLogNormal",
      logmu=logmu,
      logsigma=logsigma)
}


##-----------------------------------------------------------------------------
## Invoked by validObject() method.
validIndependentLogNormal <- function(object) {
  #cat("validating", class(object), "object", "\n")
  msg <- NULL

  ## Do checks
  if (length(object@logsigma) != 1) {
    if (length(object@logsigma) != length(object@logmu)) {
      msg <- c(msg, "lengths of logsigma and logmu differ")
    }
  }
  if (!all(object@logsigma >= 0)) {
    msg <- c(msg, "all logsigma values must be non-negative")
  }

  ## Pass or fail?
  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
}

setValidity("IndependentLogNormal", validIndependentLogNormal)


##-----------------------------------------------------------------------------
setMethod("summary", signature(object="IndependentLogNormal"),
          function(object, ...) {
  cat("An IndependentLogNormal object, representing a vector\n")
  cat(paste("of length", length(object@logmu),
              "of independent log normal random variables.\n"))
})


##-----------------------------------------------------------------------------
setMethod("rand", signature(object="IndependentLogNormal"),
          function(object, n, ...) {
  p <- length(object@logmu)
  matrix(rlnorm(n*p, object@logmu, object@logsigma), ncol=n)
})


##-----------------------------------------------------------------------------
setMethod("nrow", signature(x="IndependentLogNormal"),
          function(x) {
  length(x@logmu)
})

