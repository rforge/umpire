###############################################################
# In addition to the Independent Normal class, we also have a
# class for independent log normal, because the true expression
# value follows a log normal distribution in our model.

setClass("IndependentLogNormal",
         representation = list(
           logmu="numeric",
           logsigma="numeric"))

IndependentLogNormal <- function(logmu, logsigma) {
  if (length(logsigma) != 1 & length(logsigma) != length(logmu))
    stop("'logmu' and 'logsigma' have different lengths")
  new("IndependentLogNormal", logmu=logmu, logsigma=logsigma)
}

setMethod("summary", "IndependentLogNormal", function(object, ...) {
  cat("An IndependentLogNormal object, representing a vector\n")
  cat(paste("of length", length(object@logmu),
              "of independent log normal random variables.\n"))
})

setMethod("rand", "IndependentLogNormal", function(object, n, ...) {
  p <- length(object@logmu)
  matrix(rlnorm(n*p, object@logmu, object@logsigma), ncol=n)
})

setMethod("nrow", signature(x="IndependentLogNormal"), function(x) {
  length(x@logmu)
})

setValidity("IndependentLogNormal", function(object) {
  sizeIsRight <- length(object@logsigma) == 1 ||
                 length(object@logsigma) == length(object@logmu)
  sizeIsRight && all(object@logsigma >= 0)
})

