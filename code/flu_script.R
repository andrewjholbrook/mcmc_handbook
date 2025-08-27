### load libraries
library(ape)
library(caper)
library(MassiveMDS)

### necessary functions
source("fluR_functions.R")

# read tree file
readbeast <- function(priorRootSampleSize = 0.001) {
  # Read BEAST tree
  h1tree <- ape::read.nexus("h1_small_sample.trees")
  N <- length(h1tree$PAUP_1$tip.label)
  
  # Get tree VCV conditional on root
  treeLength <- sum(h1tree$PAUP_1$edge.length)
  treeVcv <- caper::VCV.array(h1tree$PAUP_1) / treeLength
  class(treeVcv) <- "matrix"
  
  # Integrate out fully conjugate root
  treeVcv <- treeVcv + matrix(1 / priorRootSampleSize, ncol = N, nrow = N)
  
  return(list(treeVcv = treeVcv))
}

### application 
# read tree file for h1n1
h1n1_tree_info <- readbeast()

# permuate treePrec
set.seed(666)
#set.seed(12345)
rand <- sample(1:1370, 1370)
h1n1_tree_info$treeVcvRand <- h1n1_tree_info$treeVcv[rand, rand]

# read distance matrix
fluCombi_Deff <- read.delim("fluCombi_Deff.txt", header = TRUE,
                            row.names = 1)
h1n1_dist <- fluCombi_Deff[1:1370, 1:1370]
colnames(h1n1_dist) <- rownames(h1n1_dist)
# h1n1_dist <- read.csv("h1n1_distMat.csv", header = TRUE, row.names = 1)
h1n1_distRand <- h1n1_dist[rand, rand]

# run the sampler
### FULL BMDS
# hmcresult <- hmcsampler(n_iter = 100000, burnIn = 0, data = h1n1_distRand,
#                         beast = h1n1_tree_info, StepSize = 10^-3,
#                         StepSizeSigma = 10^-3, latentDimension = 2,
#                         mdsPrecision = 1.25, sparse = FALSE,
#                         targetAccept = 0.65, targetAcceptSigma = 0.44)

### SURROGATE-TRAJECTORY HMC
sthmcresult <- hmcsampler(n_iter = 200000, burnIn = 0, data = h1n1_distRand,
                          beast = h1n1_tree_info, StepSize = 10^-3,
                          StepSizeSigma = 10^-3, latentDimension = 2,
                          mdsPrecision = 1.25, sparse = TRUE, sparse_bands = 50,
                          targetAccept = 0.65, targetAcceptSigma = 0.44)
