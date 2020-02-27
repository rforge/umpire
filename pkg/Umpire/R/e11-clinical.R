## Copyright (C) Kevin R. Coombes and Caitlin E. Coombes, 2020

## Default settings for low levels of noise
ClinicalNoiseModel <- function(shape = 1.02, scale = 0.05/shape) {
  new("NoiseModel",
      additiveOffset = 0,
      additiveScale = rgamma(1, shape=shape, scale=scale),
      multiplicativeScale = 0)
}

## percentages or fractions of each data type
mixDataTypes <- function(cont, bina, cat, nFeats){ 
  numCont <- nFeats*cont
  numBin <- nFeats*bina
  numCat <- nFeats*cat
  feats <- sample(1:nFeats, nFeats, replace=FALSE)
  contFeats <- feats[1:numCont]
  binFeats <- feats[(numCont+1):(numCont+numBin)]
  catFeats <- feats[(numCont+numBin+1):(numCont+numBin+numCat)]
  x <- cbind(contFeats, binFeats, catFeats)
  x <- as.data.frame(x)
} # end function


## nExtraBlocks = how many features do you want that are unrelated to clusters?
ClinicalBHP <- function(nFeatures,
                        nExtraBlocks = NULL,
                        meanBlockSize = NULL,
                        sigmaBlockSize =  0,
                        minBlockSize =  NULL,
                        mu0 = 6,
                        sigma0 = 1.5,
                        rate = 28.11,
                        shape = 44.25,
                        p.cor = 0.6,
                        wt.cor = 5) {
  
  BlockHyperParameters(nExtraBlocks,
                       meanBlockSize,
                       sigmaBlockSize,
                       minBlockSize,
                       mu0,
                       sigma0,
                       rate,
                       shape,
                       p.cor,
                       wt.cor
                       )
}

## Original version from CEC master's thesis.
## Arguments:
##    nFeatures = number of features
##    nClusters = number of clusters/subtypes
##    weighted = logical value used to define prevalence of subtypes
##    bHyp = object of class BlockHyperParameters, if not NULL
ClinicalEngine <- function(nFeatures, nClusters, isWeighted, bHyp = NULL){
  hitfn <- ifelse(nFeatures < 15, # small feature size
                  function(n) 2,
                  ifelse(nFeatures < 45,
                         function(n) 3, # medium feature size
                         function(n) 5 # large feature size
                         )
                  )
  ## nPossible = number of possible 'hits' based on multi-hit model
  ## of cancer, where smaller numbers are well separated
  nPos <- function(nF) {
    ifelse(nF < 12,
           nF,   # no correlation; each feature will be a seaprate hit
           ifelse (nF < 50,
                   round(nF/3),
                   ifelse(nF < 100,
                          round(nF/5),
                          20)
                   )
           )
  }

  mod <- CancerModel(name="Clinical Simulation Model (Raw)",
                     nPossible = nPos(nFeatures) # number of possible hits
                     nPattern = nClusters,       # number of subtypes
                     HIT = hitfn,
                     prevalence = Prevalence(isWeighted, nClusters)
                     )
  if (is.null(bHyp)) {
    bHyp <- ClinicalBHP(nFeatures) # use defaults based on features space
  }
  eng <- makeBlockStructure(mod, bHyp)
}


### Below here, code is taken from Simulations_More.Rmd, which was used to
### crete simulation's for CEC's master's thesis.
Prevalence <- function(weighted, k){
  if(weighted=="equal"){
    a = 100
    prev <- rdirichlet(1, rep(a, k))
    as.numeric(prev)
  } else if(weighted=="unequal"){
    if(k<=8){
      a = 10
      prev <- rdirichlet(1, rep(a, k))
      as.numeric(prev)
    } else if (k>=16){
      my.min=-1
      w = 5
      while(my.min<0.01){
      alphas <- c(rep(w*1,k/4), rep(w*2,k/4), rep(w*4,k/4), rep(w*8,k/4))
      prev <- rdirichlet(1, c(alphas))
      my.min <- min(prev)
    } # end while loop
      as.numeric(prev)
    } # end if/else cluster size for unequal weights
  } # end if/else statement
} # end function

EngineBuilder <- function(p, f, k, weighted){
  if(f < 15){
    hitfn = function(n) 2
  } else {
    hitfn = function(n) 5
  }
  mod <- CancerModel(name="Cluster Simulation Model",
                     nPossible=f/3, #number of possible 'hits' based on multi-hit theory of cancer, where smaller numbers are well separated
                     nPattern=k, #number of subtypes
                     HIT = hitfn,
                     prevalence = Prevalence(weighted, k)
                     )
  bHyp <- BlockHyperParameters(nExtraBlocks=0, #how many features do you want that are not related to the clusters? 
                               meanBlockSize=3,
                               sigmaBlockSize= 0,
                               minBlockSize= 3,
                               mu0=6,
                               sigma0=1.5,
                               rate=28.11,
                               shape=44.25,
                               p.cor=0.6,
                               wt.cor=5
                               )
  eng <- makeBlockStructure(mod, bHyp)
}


## KRC: Not needed, since "MakeData <- rand" is equivalent to this definition
MakeData <- function(eng, p){
  sim <- rand(eng, p)
}

## KRC: Probably not needed; only point is that the choice of parameters
## is taken away from the user.
Noisy <- function(sim){
  mod <- ClinicalNoiseModel(30, 40, 0.10)
  noi <- blur(mod, sim)
}

## KRC: Could probably be made more efficient
MakeBinary <- function(bi, bin, ftNum){
  ##Randomly order 1 and 0
  samp <- sample(c(0,1), replace=FALSE)
  lower <- samp[1]
  upper <- samp[2]

  if(bi[ftNum,"BI"] > my.bi){
    bound <- (bi[ftNum,"mu1"]+bi[ftNum,"mu2"])/2
    bin[ftNum,] <- ifelse(bin[ftNum,]<bound, lower, upper)
  } else {
    bound <- sample(5:35, 1)
    bound <- bound*0.01
    low.cutoff <- quantile(bin[ftNum,], bound)
    high.cutoff <- quantile(bin[ftNum,], 1-bound)

    if(bi[ftNum,"pi"]<0.5){
      bin[ftNum,] <- ifelse(bin[ftNum,]<low.cutoff, lower, upper)
    } else {
      bin[ftNum,] <- ifelse(bin[ftNum,]<high.cutoff, lower, upper)
    } # end inner if/else (bi[ftNum, "pi"]<0.5)
  } # end outer if/else (bi[ftNum, "BI] > my.bi)
  bin
} # end function

LabelBinary <- function(bin, ftNum, Type){
  per1 <- sum(bin[ftNum,])/length(bin[ftNum,])
  Type <- ifelse(per1>0.9 | per1<0.1, "asymm", "symm")
  Type
}

# Determine if the feature is nominal or ordinal, and store that information in a vector.
catType <- sample(c(1,2), 1, replace=TRUE) #1 is ordinal, 2 is nominal


MakeCategories <- function(bin, ftNum, nPts, catType){
  
# Sample a number of categories between 3 and 9.
cc <- sample(3:9, 1, replace=TRUE)

# Sample bin sizes from the Dirichlet distribution
a <- 20
r <- rdirichlet(1, rep(a, cc))

# Construct a vector of percentage cutoffs
cuts <- r[1]
for(i in 2:cc){
  cuts <- c(cuts, cuts[i-1]+r[i])
}
cuts <- c(0, round(cuts*nPts))

# Create a numbered list of bin identities.
if(catType == 1){#ordinal
  id <- 1:cc
} else if(catType == 2){#nominal
  id <- sample(1:cc, cc, replace=FALSE)
}

# Link percentage cutoffs to feature number, slice, and assign bin identities.
ranks <- order(bin[ftNum,])

for(kk in 1:cc){
  #print(kk)
for(i in 1:nPts){
above <- ranks[i] > cuts[kk]
if(above == TRUE){
  bin[ftNum,i] <- id[kk]
  
}#close if
}#close for(i)
}#close for(kk)

bin
}#close the function

LabelCategories <- function(catType){
  if(catType == 1){#ordinal
    Type <- "ordinal"
  } else if(catType == 2){#nominal
    Type <- "nominal"
  }
  Type
}

## KRC: Begining version of "chooseDataTypes", though the temptation
## to call it "Alligator" is hard to resist.
Allocator <- function(cont, bina, cat, nFeats){ #percentages or fractions of each data type
  numCont <- nFeats*cont
  numBin <- nFeats*bina
  numCat <- nFeats*cat
  feats <- sample(1:nFeats, nFeats, replace=FALSE)
  contFeats <- feats[1:numCont]
  binFeats <- feats[(numCont+1):(numCont+numBin)]
  catFeats <- feats[(numCont+numBin+1):(numCont+numBin+numCat)]
  x <- cbind(contFeats, binFeats, catFeats)
  x <- as.data.frame(x)
}#close function


