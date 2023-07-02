setwd("~/mcmc_handbook/")

#
####
####### HMC with adaptive stepsize (0.65 acceptance rate target)
####
#

library(coda)
library(parallel)

target <- function(theta,offset=6) {
  # theta is D vector
  # distribution is multimodal in first dimension
  D <- length(theta)
  output <- log(0.5 * exp(- sum(theta[1]^2)/2 ) + 0.5 * exp(-sum((theta[1]-offset)^2)/2) ) +
    - 0.5 * sum(theta[2:D]^2)
  return(output)
}

grad <- function(theta,offset=6) {
  D <- length(theta)
  output1 <- 1/(0.5 * exp(- sum(theta[1]^2)/2 ) + 0.5 * exp(-sum((theta[1]-offset)^2)/2) )* 
    ( (-theta[1])*exp(- sum(theta[1]^2)/2 ) + (-theta[1]+offset) *exp(- sum((theta[1]-offset)^2)/2) )
  output <- c(output1,-theta[2:D])
  return(output)
}

leapfrog <- function(proposalState,momentum,stepSize,L_jitter) {
  momentum <- momentum + 0.5 * stepSize * grad(proposalState)
  for (l in 1:L_jitter) {
    proposalState <- proposalState + stepSize * momentum
    if (l!=L_jitter) momentum <- momentum + stepSize * grad(proposalState)
  }
  momentum <- momentum + 0.5 * stepSize * grad(proposalState)
  return(list(proposalState,momentum))
}


delta <- function(n) {
  return( min(0.01,n^(-0.5)) )
}

adapt_hmc <- function(D, maxIts, targetAccept=0.65, stepSize=1, L=20) {
  
  Us <- runif(maxIts)
  
  chain <- matrix(0,maxIts,D)
  #stepSize <- 1
  chain[1,] <- rnorm(D)
  currentU  <- - target(chain[1,])
  
  totalAccept <- rep(0,maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0
  
  for (i in 2:maxIts) {
    proposalState    <- chain[i-1,]
    momentum         <- rnorm(D)
    currentK   <- sum(momentum^2)/2
    
    # leapfrog steps
    L_jitter <- round(L * rnorm(1) ) 
    
    leapOut <- leapfrog(proposalState,momentum,stepSize,L_jitter)
    proposalState <- leapOut[[1]]
    momentum      <- leapOut[[2]]
    
    # quantities for accept/reject
    proposedU = - target(proposalState)
    proposedK = sum(momentum^2)/2

    if (log(Us[i]) < currentU + currentK - proposedU - proposedK) {
      chain[i,]   <- proposalState
      currentU    <- proposedU
      totalAccept[i] <- 1
      Acceptances = Acceptances + 1
    } else {
      chain[i,] <- chain[i-1,]
    }
    
    SampCount <- SampCount + 1

    # tune
    if (SampCount == SampBound) { 
      AcceptRatio <- Acceptances / SampBound
      if ( AcceptRatio > targetAccept ) {
        stepSize <- stepSize * (1 + delta(i-1))
      } else {
        stepSize <- stepSize * (1 - delta(i-1))
      }
  
      SampCount <- 0
      Acceptances <- 0
    }
    
    
    if (i %% 1000 == 0) cat("Iteration ", i,"\n","stepSize: ", stepSize, "\n") 
  }
  
  cat("Acceptance rate: ", sum(totalAccept)/(maxIts-1))
  return(chain)
}


multiprop_hmc <- function(D, maxIts, targetAccept=0.65, stepSize=1, L=20, P=100) {

  chain <- matrix(0,maxIts,D)
  #stepSize <- 1
  chain[1,] <- rnorm(D)
  currentU  <- - target(chain[1,])

  totalAccept <- rep(0,maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0

  for (i in 2:maxIts) {
    proposalState    <- chain[i-1,]
    momentum         <- rnorm(D)
    currentK   <- sum(momentum^2)/2

    # initial leapfrog steps
    L_jitter <- round(L * rnorm(1) )
    leapOut <- leapfrog(proposalState,momentum,stepSize,L_jitter)
    middle <- leapOut[[1]]
    momentum      <- leapOut[[2]]

    # get magnitude of momentum and P random momenta
    momMag <- sqrt(sum(momentum^2))
    intermediateMoms <- matrix(0,D,P)
    for (j in 1:P) {
      gaus <- rnorm(D)
      gaus <- gaus/sqrt(sum(gaus^2))*momMag
      intermediateMoms[,j] <- gaus
    }

    energies      <- rep(0,P+1)
    energies[1]   <- currentU + currentK
    proposals     <- matrix(0,D,P+1)
    proposals[,1] <- chain[i-1,]

    for(j in 1:P) {
      proposalState   <- middle
      momentum        <- intermediateMoms[,j]
      leapOut         <- leapfrog(proposalState,momentum,stepSize,L_jitter)
      proposalState   <- leapOut[[1]]
      momentum        <- leapOut[[2]]
      proposedU       <- - target(proposalState)
      proposedK       <- sum(momentum^2)/2
      energies[j+1]   <- proposedU + proposedK
      proposals[,j+1] <- proposalState
    }

    gumbs <- - log ( -log(runif(P+1)) )
    selection <- which.max(-energies+gumbs)

    chain[i,]   <- proposals[,selection]
    currentU    <- - target(chain[i,]) # can be made faster
    if(selection != 1) {
      totalAccept[i] <- 1
      Acceptances = Acceptances + 1
    }

    SampCount <- SampCount + 1

    # tune
    # if (SampCount == SampBound) {
    #   AcceptRatio <- Acceptances / SampBound
    #   if ( AcceptRatio > targetAccept ) {
    #     stepSize <- stepSize * (1 + delta(i-1))
    #   } else {
    #     stepSize <- stepSize * (1 - delta(i-1))
    #   }
    #
    #   SampCount <- 0
    #   Acceptances <- 0
    # }


    if (i %% 1000 == 0) cat("Iteration ", i,"\n","stepSize: ", stepSize, "\n")
  }

  cat("Acceptance rate: ", sum(totalAccept)/(maxIts-1))
  return(chain)
}

generate_prop <- function(propIndex,middle,momMag,stepSize,L_jitter) {
  propIndex2 <- propIndex # suppress propIndex not used warning

  gaus <- rnorm(D)
  gaus <- gaus/sqrt(sum(gaus^2))*momMag
  intermediateMom <- gaus
  
  proposalState   <- middle
  momentum        <- intermediateMom
  leapOut         <- leapfrog(proposalState,momentum,stepSize,L_jitter)
  proposalState   <- leapOut[[1]]
  momentum        <- leapOut[[2]]
  proposedU       <- - target(proposalState)
  proposedK       <- sum(momentum^2)/2
  energy          <- proposedU + proposedK
  proposal        <- proposalState

  return(list(energy,proposal))
}

multiprop_hmc_forked <- function(D, maxIts, targetAccept=0.65, stepSize=1, L=20, P=100, nCores=4) {
  
  chain <- matrix(0,maxIts,D)
  chain[1,] <- rnorm(D)
  currentU  <- - target(chain[1,])
  
  totalAccept <- rep(0,maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0
  
  for (i in 2:maxIts) {
    proposalState    <- chain[i-1,]
    momentum         <- rnorm(D)
    currentK   <- sum(momentum^2)/2
    
    # initial leapfrog steps
    L_jitter <- round(L * rnorm(1) )
    leapOut <- leapfrog(proposalState,momentum,stepSize,L_jitter)
    middle <- leapOut[[1]]
    momentum      <- leapOut[[2]]
    
    # get magnitude of momentum and P random momenta
    momMag <- sqrt(sum(momentum^2))
    energies      <- rep(0,P+1)
    energies[1]   <- currentU + currentK
    proposals     <- matrix(0,D,P+1)
    proposals[,1] <- chain[i-1,]
    
    # for(j in 1:P) {
    #   propOut <- generate_prop(j,
    #                            middle=middle,
    #                            momMag=momMag,
    #                            stepSize=stepSize,
    #                            L_jitter=L_jitter)
    #   energies[j+1]   <- propOut[[1]]
    #   proposals[,j+1] <- propOut[[2]]
    # }
    propOut <- mclapply(1:P,FUN=generate_prop, middle=middle,
                        momMag=momMag,stepSize=stepSize,L_jitter=L_jitter,
                        mc.cores = nCores)
    for (j in 1:P) {
      energies[j+1] <- propOut[[j]][[1]]
      proposals[,j+1] <- propOut[[j]][[2]]
    }
    
    gumbs <- - log ( -log(runif(P+1)) )
    selection <- which.max(-energies+gumbs)
    
    chain[i,]   <- proposals[,selection]
    currentU    <- - target(chain[i,]) # can be made faster
    if(selection != 1) {
      totalAccept[i] <- 1
      Acceptances = Acceptances + 1
    }
    
    SampCount <- SampCount + 1
    
    # tune
    # if (SampCount == SampBound) {
    #   AcceptRatio <- Acceptances / SampBound
    #   if ( AcceptRatio > targetAccept ) {
    #     stepSize <- stepSize * (1 + delta(i-1))
    #   } else {
    #     stepSize <- stepSize * (1 - delta(i-1))
    #   }
    #
    #   SampCount <- 0
    #   Acceptances <- 0
    # }
    
    
    if (i %% 1000 == 0) cat("Iteration ", i,"\n","stepSize: ", stepSize, "\n")
  }
  
  cat("Acceptance rate: ", sum(totalAccept)/(maxIts-1))
  return(chain)
}

numberJumps <- function(chain,cutoff=3) {
  N <- dim(chain)[1]
  count <- 0
  
  for (i in 2:N) {
    if ( (chain[i-1,1] < cutoff & chain[i,1] > cutoff) |
         (chain[i-1,1] > cutoff & chain[i,1] < cutoff) ) {
      count <- count + 1
    }
  }
  return(count)
}

################################################################################

D <- 10000
set.seed(1)

for(nProp in 2^(2:10)) {
  
  ptm <- proc.time()[3]
  results <- adapt_hmc(D=D,
                       maxIts=100000,
                       stepSize = 0.22,
                       targetAccept = 0.65,
                       L=20)
  tm1    <- proc.time()[3] - ptm
  nj1    <- numberJumps(results)
  ess1   <- effectiveSize(results[,1])
  ess1_2 <- effectiveSize(results[,D])
  cat(1, tm1, nj1, ess1, ess1_2, "\n",
      file = "output/parallel_results.txt", append = TRUE)
  
  ptm <- proc.time()[3]
  results2 <- multiprop_hmc_forked(D=D,
                                   maxIts=100000,
                                   stepSize = 0.22,
                                   L=20,
                                   P=nProp,
                                   nCores=100)
  tm2 <- proc.time()[3] - ptm
  nj2    <- numberJumps(results2)
  ess2   <- effectiveSize(results2[,1])
  ess2_2 <- effectiveSize(results2[,D])
  cat(nProp, tm2, nj2, ess2, ess2_2, "\n",
      file = "output/parallel_results.txt", append = TRUE)
  
}


