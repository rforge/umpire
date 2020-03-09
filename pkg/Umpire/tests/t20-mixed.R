library(Umpire)

ce <- Umpire:::ClinicalEngine(20, 4, TRUE)
summary(ce)

set.seed(97531)
dset <- rand(ce, 300)
class(dset)
names(dset)
summary(dset$clinical)
dim(dset$data)

cnm <- Umpire:::ClinicalNoiseModel() # default
noisy <- blur(cnm, dset$data)

dt <- Umpire:::setDataTypes(dset$data, 1, 1, 1, 0.3)
cp <- dt$cutpoints

hasBreaks <- all (hb <- sapply(cp, function(X) {
  "breaks" %in% names(X)
}))
summary(hb)
hasLabels <- all (hl <- sapply(cp, function(X) {
  "labels" %in% names(X)
}))
summary(hl)
hasType <- all (ht <- sapply(cp, function(X) {
  "Type" %in% names(X)
}))
summary(ht)

type <- sapply(cp, function(X) { X$Type })
table(type)
sum(is.na(type))
length(type)

mte <- new("MixedTypeEngine",
           ce,
           noise = cnm,
           cutpoints = dt$cutpoints)

R <- rand(mte, 20)
