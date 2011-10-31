###############################################################
# This is tricky.  A CANCER ENGINE is supposed to combine the
# CANCER MODEL (which defines hit patterns and keeps track of
# how each hit is a latent variable linking to survival or
# outcome phenotype data) with two gene expression ENGINE
# objects (that can be used to generate data depending on
# the presence or absence of individual hits).  The tricky part
# arises because we may want to apply this class with either
# an ENGINE or an ENGINE WIT ACTIVITY,and so it is hard to
# define the class represenation.  The first attempt used a
# poor man's refernce schem, by storing character strings with
# the names of objects that could be retrieved by an
# eval(as.name(...)) construction.  This idea fails if you try
# to produce a CANCER ENGINE inside a function, where R fails
# to locate the objects.  So, we replaced that with a local
# environment that stores the actual ENGINEs.


setClass("CancerEngine",
         representation=list(
           cm="CancerModel",
           base="character",
           altered="character",
           localenv="environment"))

CancerEngine <- function(cm, base, altered) {
  localenv <- new.env()
  if (is.character(base)) {
    x <- try(eval(as.name(base), parent.frame()))
    if (inherits(x, "try-error")) {
      stop(paste("Unable to locate base engine via", base))
    }
    if(!inherits(x, "Engine")) {
      stop(paste('base argument (', base, ") must evaluate to an Engine"))
    }
    assign(base, x, envir=localenv)
  } else {
    assign("base", base, envir=localenv)
    base <- "base"
  }
  if (is.character(altered)) {
    x <- try(eval(as.name(altered), parent.frame()))
    if (inherits(x, "try-error")) {
      stop(paste("Unable to locate altered engine via", altered))
    }
    if(!inherits(x, "Engine")) {
      stop(paste('altered argument (', altered, ") must evaluate to an Engine"))
    }
    assign(altered, x, envir=localenv)
  } else {
    assign("altered", altered, envir=localenv)
    altered <- "altered"
  }
  new("CancerEngine", cm=cm, base=base, altered=altered, localenv=localenv)
}

# sigh.  It would be nice if R were _really_ object-oriented.  I'd like to
# dispatch on the second argument of rand, but that doesn't really work.
# So, we are going to assume that "n" is a list of cancer subtypes instead
# of a simple integer.

setMethod("rand", "CancerEngine", function(object, n, ...) {
  hitlist <- n                            # remember the subtypes
  n <- length(hitlist)                    # convert back to number of vectors to generate
  B <- get(object@base, envir=object@localenv)    # base Engine
  A <- get(object@altered, envir=object@localenv) # altered Engine
  b <- rand(B, n)                             # base simulation
  a <- rand(A, n) # altered values
  temp <- object@cm@hitPattern
  # there ought to be a better way to do this
  # idea is to expand the "hit pattern" for each patient ot include
  # all genes in a correlated block
  U <- unlist(lapply(B@components, nrow))      # size of each component
  ends <- cumsum(U)
  starts <- (1+c(0, ends))[1:length(U)]
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



