#==============================================================================#
#                                                                              #
#                            SPATIAL CAPTURE-RECAPTURE                         #
#                         Jackdaw-mediated removal of samples                  #
#                                 José Jiménez                                 #
#                               22:37 17/05/2025                               #
#                                                                              #
#==============================================================================#

# Set working directory
setwd('C:/...')

# Load custom SCR functions
source("Funciones_SCR.R")

# Load required libraries
library(spatstat)
library(lattice)
library(coda)

# Define simulation parameters
N <- 10  # Population size
K <- 5   # Number of capture occasions

# Define three jackdaws colony locations
set.seed(123)
colonies <- cbind(runif(3,1,8), runif(3,1,8))

# Function to simulate SCR data under a Poisson model
simSCR0Poisson <- function(discard0 = TRUE, N=N, K=25, p0=0.1, sigma=1,
                           buffer=buffer, array3d = TRUE, rnd = 2013) {
  set.seed(rnd)
  traplocs <- cbind(sort(rep(1:8, 8)), rep(1:8, 8))  # 8x8 trap grid
  Dmat <- e2dist(traplocs, traplocs)
  ntraps <- nrow(traplocs)
  plot(traplocs)

  # Define state space limits
  buffer <- 2.5 * sigma
  Xl <- min(traplocs[, 1] - buffer)
  Xu <- max(traplocs[, 1] + buffer)
  Yl <- min(traplocs[, 2] - buffer)
  Yu <- max(traplocs[, 2] + buffer)

  # Simulate individual activity centers
  sx <- runif(N, Xl, Xu)
  sy <- runif(N, Yl, Yu)
  S <- cbind(sx, sy)

  # Plot traps and activity centers
  plot(traplocs, pch=3, col='blue', xlim=c(Xl, Xu), ylim=c(Yl, Yu))
  points(S, pch=16, col='red')

  # Calculate detection probabilities
  D <- e2dist(S, traplocs)
  alpha1 <- 1 / (2 * sigma^2)
  alpha0 <- log(p0)
  muy <- exp(alpha0) * exp(-alpha1 * D^2)

  # Simulate detections
  Y <- matrix(NA, nrow = N, ncol = ntraps)
  for (i in 1:nrow(Y)) {
    Y[i, ] <- rpois(ntraps, K * muy[i, ])
  }

  # Optionally discard individuals with no detections
  if (discard0) {
    totalcaps <- apply(Y, 1, sum)
    Y <- Y[totalcaps > 0, ]
  }

  dimnames(Y) <- list(1:nrow(Y), paste("trap", 1:ncol(Y), sep = ""))

  # Convert to 3D array if requested
  if (array3d) {
    Y <- array(NA, dim = c(N, K, ntraps))
    for (i in 1:nrow(Y)) {
      for (j in 1:ntraps) {
        Y[i, 1:K, j] <- rpois(K, muy[i, j])
      }
    }
    Y <- aperm(Y, c(1, 3, 2))

    # Remove 75% of samples at specific traps
    Y[, c(14,15,16,22,23,24,30,31,32,39,40,46,47,48,
          54,55,56,63,64, 9,17,18,19,25,26,27,33,34), ] <- 0

    if (discard0) {
      Y2d <- apply(Y, c(1, 3), sum)
      ncaps <- apply(Y2d, 1, sum)
      Y <- Y[ncaps > 0, , ]
    }
  }

  return(list(Y = Y, traplocs = traplocs, xlim = c(Xl, Xu), ylim = c(Yl, Yu),
              N = N, p0 = p0, sigma = sigma, S = S, K = K, buffer = buffer))
}

  # Run 100 simulations
for (sim in 1:100) {
  data <- simSCR0Poisson(N=N, K=K, array3d = TRUE, p0=0.1, sigma=1, discard0=FALSE, rnd=sim)

  S <- data$S
  y3d <- data$Y

  # Identify captured individuals
  captured <- apply(y3d, 1, sum)
  captured[captured > 1] <- 1
  captured <- (1:N) * captured

  # Remove individuals with no detections
  y3d <- y3d[!apply(y3d, 1, sum) == 0, , ]
  y <- apply(y3d, c(1, 2), sum)

  # Summary statistics
  dim(y)[1]  # Number of captured individuals
  sum(y)     # Total number of detections

  detections <- rowSums(y) > 0
  nind <- nrow(y)
  X <- data$traplocs
  K <- data$K
  J <- nrow(X)
  KT <- rep(K, J)

  # Define state space
  xlim <- data$xlim
  ylim <- data$ylim
  area <- diff(xlim) * diff(ylim)

  # Plot detections and activity centers
  plot(X, xlim=xlim, ylim=ylim, pch="+", cex=1, col="blue", asp=TRUE)
  points(S, pch=1, col="red", cex=1.25)
  points(S[captured,], pch=16, col="red", cex=1.25)
  datn <- apply(y3d, c(2,3), sum)
  tot <- apply(datn, 1, sum)
  symbols(X, circles=tot/5, inches=F, bg="#00000022", fg=NULL, add=T)
  points(X, pch="+", cex=1, col="blue")

  # Spiderplot of detections
  spiderplotJJ4(y3d, X, buffer=2, lwd=2)
  symbols(colonies[1,1], colonies[1,2], circles = 2, inches = FALSE, add = TRUE)
  symbols(colonies[2,1], colonies[2,2], circles = 2, inches = FALSE, add = TRUE)
  symbols(colonies[3,1], colonies[3,2], circles = 2, inches = FALSE, add = TRUE)

  # Data augmentation
  M <- 150
  detected <- rowSums(y) > 0
  n0 <- sum(detected)

  # Initial values for activity centers
  sst <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
  for (i in 1:nind) {
    sst[i,1] <- mean(X[y[i,]>0,1])
    sst[i,2] <- mean(X[y[i,]>0,2])
  }

  # Load nimble and define model
  library(nimble)
  #-----------------------------#
  #     NIMBLE MODEL SETUP     #
  #-----------------------------#

  # Define the hierarchical Bayesian model using nimbleCode
  code <- nimbleCode({

    # Priors for detection parameters
    p0 ~ dunif(0,1)              # Baseline detection probability
    alpha1 ~ dnorm(0,0.01)       # Prior for alpha1 (related to sigma)
    sigma <- sqrt(1/(2*alpha1))  # Detection scale parameter
    psi ~ dunif(0,1)             # Inclusion probability for data augmentation
    # Loop over all augmented individuals
    for(i in 1:M){
      z[i] ~ dbern(psi)  # Latent indicator: 1 if individual is part of population

      # Prior for activity center locations
      s[i,1] ~ dunif(xlim[1], xlim[2])
      s[i,2] ~ dunif(ylim[1], ylim[2])

      # Calculate detection rate for each trap
      p[i,1:J] <- GetDetectionRate(s = s[i,1:2], 
                                   X = X[1:J,1:2], 
                                   J = J,
                                   sigma = sigma, 
                                   p0 = p0, 
                                   z = z[i])
    }

    # Observation model for detected individuals
    for(i in 1:n0){
      y[i,1:J] ~ dPoissonVector(lambda = p[i,1:J] * K)
    }

    for(i in (n0+1):M) {
      zeros[i] ~ dpois(sum(p[i,1:J]) * K)
    }

    # Derived quantities
    N <- sum(z[1:M])  # Total population size
    D <- N / area     # Density
  })

  #-----------------------------------------#
  #     CUSTOM FUNCTION AND DISTRIBUTIONS   #
  #-----------------------------------------#

  GetDetectionRate <- nimbleFunction(
    run = function(s = double(1), p0 = double(0), sigma = double(0), 
                   X = double(2), J = double(0), z = double(0)) {
      returnType(double(1))
      if(z == 0) return(rep(0, J))  # No detection if individual is not present
      if(z == 1){
        d2 <- ((s[1] - X[1:J,1])^2 + (s[2] - X[1:J,2])^2)
       ans <- p0 * exp(-d2 / (2 * sigma^2))
        return(ans)
      }
    }
  )

  dPoissonVector <- nimbleFunction(
    run = function(x = double(1), lambda = double(1), log = integer(0, default = 0)) {
      J <- length(x)
      ans <- 0.0
      for(j in 1:J)
        ans <- ans + dpois(x[j], lambda[j], log = TRUE)
      returnType(double())
      if(log) return(ans)
      else return(exp(ans))
    })

  rPoissonVector <- nimbleFunction(
    run = function(n = integer(), lambda = double(1)) {
      J <- length(lambda)
      ans <- numeric(J)
      for(j in 1:J)
        ans[j] <- rpois(1, lambda[j])
      returnType(double(1))
      return(ans)
    })


  registerDistributions(list(
    dPoissonVector = list(
      BUGSdist = "dPoissonVector(lambda)",
      Rdist = "dPoissonVector(lambda)",
      discrete = TRUE,
      range = c(0, Inf),
      types = c('value = double(1)', 'lambda = double(1)')
    )
  ))

  #-----------------------------#
  #     MODEL INITIALIZATION   #
  #-----------------------------#

  # Constants for the model
  constants <- list(M = M, K = K, J = J, n0 = n0, area = area)

  # Data for the model
  data <- list(
    y = y, 
    X = X, 
    xlim = xlim, 
    ylim = ylim,
    zeros = c(rep(NA, n0), rep(0, M - n0)),
    z = c(rep(1, n0), rep(NA, M - n0))
  )

  # Initial values
  inits <- list(
    p0 = 0.5, 
    alpha1 = 1, 
    s = sst,
    z = c(rep(NA, n0), rep(1, M - n0))
  )

  #-----------------------------#
  #     COMPILE AND RUN MCMC   #
  #-----------------------------#

  # Build the model
  Rmodel <- nimbleModel(code = code, constants = constants, data = data, inits = inits,
                        calculate = FALSE, check = FALSE)

  # Compile the model for faster execution
  Cmodel <- compileNimble(Rmodel)

  # Configure MCMC
  conf <- configureMCMC(Rmodel, monitors = c('N', 'D', 'sigma', 'psi', 'p0'))

  # Replace default samplers for activity centers with block samplers
  conf$removeSamplers("s")
  ACnodes <- paste0("s[", 1:constants$M, ", 1:2]")
  for(node in ACnodes) {
    conf$addSampler(target = node,
                    type = "RW_block",
                    control = list(adaptScaleOnly = TRUE),
                    silent = TRUE)
  }

  # Build and compile the MCMC object
  MCMC <- buildMCMC(conf)
  CompMCMC <- compileNimble(MCMC, project = Rmodel)

  # Run the MCMC
  nb <- 1000      # Burn-in iterations
  ni <- 5000 + nb # Total iterations
  nc <- 3         # Number of chains

  start.time2 <- Sys.time()
  outNim <- runMCMC(CompMCMC, niter = ni, nburnin = nb, nchains = nc, inits = inits,
                    setSeed = TRUE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  end.time <- Sys.time()
  end.time - start.time2  # Execution time

  # Save the output
  save(outNim, file = paste("out", sim, "Grajillas.RData", sep = ""))
  print(sim)
}
