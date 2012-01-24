NormalVsCancerModel <- function(nBlocks, name="NormalVsCancer") {
  call <- match.call()
  hp <- matrix(rep(0:1, each=nBlocks), ncol=2)
  s <- rnorm(nBlocks, 0, 2)
  o <- rnorm(nBlocks, 0, 2)
  prevalence=c(0.5, 0.5)
  sm <- SurvivalModel()
  new("CancerModel",
      name=name,
      hitPattern=hp,
      survivalBeta=s,
      outcomeBeta=o,
      prevalence=prevalence,
      survivalModel=sm,
      call=call)
}

NormalVsCancerEngine <- function(nBlocks, hyperp=NULL) {
  if (is.null(hyperp)) hyperp <- BlockHyperParameters()
  model <- NormalVsCancerModel(nBlocks)
  makeBlockStructure(model, hyperp)
}
