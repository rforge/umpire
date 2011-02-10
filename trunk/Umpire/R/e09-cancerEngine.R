###############################################################

setClass("CancerEngine",
         representation=list(
           cm="CancerModel",
#           hitlist="numeric",
           base="character",
           altered="character"))

# sigh.  It would be nice if R were _really_ object-oriented.  I'd like to
# dispatch on the second argument of rand, but that doesn't really work.
# So, we are going to assume that "n" is a list of cancer subtypes instead
# of a simple integer.

setMethod("rand", "CancerEngine", function(object, n, ...) {
  hitlist <- n                            # remember the subtypes
  n <- length(hitlist)                    # convert back to number of vectors to generate
  B <- eval(as.name(object@base))         # base Engine
  U <- unlist(lapply(B@components, nrow)) # size of each component
  ends <- cumsum(U)
  starts <- (1+c(0, ends))[1:length(U)]
  b <- rand(B, n)                             # base simulation
  a <- rand(eval(as.name(object@altered)), n) # altered values
  temp <- object@cm@hitPattern
  # there ought to be a better way to do this
  # idea is to expand the "hit pattern" for each patient ot include
  # all genes in a correlated block
  isHit <- matrix(0, nrow=sum(U), ncol=n)
  for (i in 1:nrow(temp)) {
    for (j in 1:n) {
      if (temp[i, hitlist[j]]==1) {
        isHit[starts[i]:ends[i], j] <- 1
      }
    }
  }
  b*(1-isHit) + a*(isHit)
})



