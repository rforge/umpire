###############################################################

setClass("CancerEngine",
         representation=list(
           base="character",
           altered="character"))

# The fiction that R has an object-oriented aspect is sometimes amazing.
# I'm going to dispatch on the first argument of the rand function.
# As part of the description of rand, the semantics of 'n' suggest that
# it should be the number of sample (column) vectors to be generated.
# Our model of cancer, however, says it is heterogeneous, in that
# different hits affect different individuals.
# We have two implementation options. We can continue to use 'n' as
# an integer and pass the hit patterns in as another argument, or we can
# just use 'n' for the hit patterns and compute the number of samples
# from that.

# incomplete; does not work
if(0) {
setMethod("rand", "CancerEngine", function(object, n, ...) {
  base <- rand(eval(as.name(object@base)), n)
  altered <- rand(eval(as.name(object@affected)), n)
  isHit <- expand(hitlist)
  base*(1-isHit) + altered*(isHit)
})
}


