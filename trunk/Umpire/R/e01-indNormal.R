###############################################################
# In practice, we use components that are multivariate normal.
# The case of independent normals is treated separately, but
# blocks of independent normals are components in the sense
# described above. Note that a component can be viewed as a
# special case of an engine that only has one component, so
# most methods of engines must also be defined for individual
# components.

##=============================================================================
setClass("IndependentNormal",
         representation(mu="numeric",
                        sigma="numeric"))


##-----------------------------------------------------------------------------
## Generates an IndependentNormal object.
IndependentNormal <- function(mu, sigma) {
  if (length(sigma) != 1 & length(sigma) != length(mu))
    stop("'mu' and 'sigma' have different lengths")
  new("IndependentNormal",
      mu=mu,
      sigma=sigma)
}


##-----------------------------------------------------------------------------
## Invoked by validObject() method.
validIndependentNormal <- function(object) {
  #cat("validating", class(object), "object", "\n")
  msg <- NULL

  ## Do checks
  if (length(object@sigma) != 1) {
    if (length(object@sigma) != length(object@mu)) {
      msg <- c(msg, "lengths of sigma and mu differ")
    }
  }
  if (!(object@sigma >= 0)) {
    msg <- c(msg, "sigma must be non-negative")
  }

  ## Pass or fail?
  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
}

setValidity("IndependentNormal", validIndependentNormal)


##-----------------------------------------------------------------------------
setMethod("summary", signature(object="IndependentNormal"),
          function(object, ...) {
  cat("An IndependentNormal object, representing a vector\n")
  cat(paste("of length", length(object@mu),
            "of independent normal random variables.\n"))
})


##-----------------------------------------------------------------------------
setMethod("rand", signature(object="IndependentNormal"),
          function(object, n, ...) {
  p <- length(object@mu)
  matrix(rnorm(n*p, object@mu, object@sigma), ncol=n)
})


##-----------------------------------------------------------------------------
setMethod("nrow", signature(x="IndependentNormal"),
          function(x) {
  length(x@mu)
})

