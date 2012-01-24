###############################################################
# A NOISE MODEL represents the additional machine noise that is layered
# on top of any biological variabilty when measuring the gene expression
# in a set of samples. We model both additive and multiplicative noise,
# so that the observed expression of gene g in sample i is given by
#    Y_gi = S_gi exp(H_gi) + E_gi
# where
#    Y_gi = observed expression
#    S_gi = true bilogical signal
#    H_gi ~ N(0, phi) defines the multiplicative noise
#    E_gi ~ N(nu, tau) defines the additive noise
# Note that we allow a systematic offset/bias in the additive noise model.

setClass("NoiseModel",
         representation = list(
           additiveOffset = "numeric",
           additiveScale = "numeric",
           multiplicativeScale = "numeric"))

setValidity("NoiseModel", function(object) {
  object@additiveScale >= 0 && object@multiplicativeScale >= 0
})

NoiseModel <- function(nu, tau, phi) {
  new("NoiseModel",
      additiveOffset = nu,
      additiveScale = tau,
      multiplicativeScale = phi)
}

# 'blur' is the main method associated with a noise model.  The main
# operation is given by blur(noiseModel, dataMatrix), which adds
# and multiplies random noise to the dataMatrix containing the true
# signal.

setMethod("blur", "NoiseModel", function(object, x, ...) {
  if(inherits(x, "matrix")) {
    h <- matrix(rnorm(nrow(x)*ncol(x), 0, object@multiplicativeScale),
                nrow=nrow(x))
    e <- matrix(rnorm(nrow(x)*ncol(x), object@additiveOffset,
                      object@additiveScale),
                nrow=nrow(x))
  } else {
    h <- rnorm(length(x), 0, object@multiplicativeScale)
    e <- rnorm(length(x), object@additiveOffset,
               object@additiveScale)
  }
  exp(h)*x + e
})

