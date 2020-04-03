library(Umpire)

## Need to generate a gene expression data set to get started.
## So, we must start with a cancer engine
ce <- ClinicalEngine(20, 4, TRUE)
summary(ce)
nrow(ce@cm) # bugged - go bakc and fix
nrow(ce@localenv$eng) # correct
nComponents(ce@localenv$eng) # correct

## Now make a data set
set.seed(97531)
dset <- rand(ce, 300)
class(dset)
names(dset)
summary(dset$clinical)
dim(dset$data) # 50 features, 300 samples

## Must add noise before making a clinical engine
cnm <- ClinicalNoiseModel() # default
noisy <- blur(cnm, dset$data)

## Now we set the data types
dt <- setDataTypes(dset$data, 1/3, 1/3, 1/3, 0.3)
cp <- dt$cutpoints
type <- sapply(cp, function(X) { X$Type })
table(type)
sum(is.na(type))
length(type)
class(dt$binned)
dim(dt$binned)
summary(dt$binned)

## Use the pieces from above to create an MTE.
mte <- new("MixedTypeEngine",
           ce,
           noise = cnm,
           cutpoints = dt$cutpoints)
# and generate some data
R <- rand(mte, 20)
summary(R)

