# Goal: Estimate 1st stage model
# Dependencies: "datamodel.Rdata"

setwd("~/Replication Package")

library(rjags)
library(coda)
load.module("glm")

library(foreach)
library(doMC)  
registerDoMC(4) 

RNGkind("L'Ecuyer-CMRG") 

set.seed(100)

# Load recoded data

load("datamodel.Rdata")

attach(datamodel)

n <- dim(datamodel)[1]

fixed.unconv <- rep(0, n)
fixed.unconv[unconv.scale == 3 & conv.scale == 0] <- 2
fixed.unconv[unconv.scale == 0 & conv.scale > 6] <- 1
which(fixed.unconv != 0)

fixed.conv <- rep(0, n)
fixed.conv[unconv.scale == 3 & conv.scale == 0] <- 1
fixed.conv[unconv.scale == 0 & conv.scale > 6] <- 2
which(fixed.conv != 0)

Y <- cbind(act.meet.mun, act.cont.mun, act.cont.aut, act.meet.imp, act.meet.pty, act.meet.aso, act.solve.prob, act.work.pty, act.protest, act.strike, act.block)

ntypes <- 2

nact <- dim(Y)[2]

nconv <- nact - 3

wprior <- c(1, 1)

##---------- BUGS/JAGS code ----------##

modelString = "

model{

  for (i in 1:n) { 

    for (j in 1:nact) {

      Y[i,j] ~ dbern(p2[i,j])
      p2[i,j] <- max(0, min(1,p[i,j]))

      logit(p[i,j]) <- alpha[j] + alpha.conv[j] * (conv.type[i] - 1) + alpha.unconv[j] * (unconv.type[i] - 1)
    }

    conv.type[i] <- ifelse(fixed.conv[i] != 0, fixed.conv[i], rconv.type[i])
    rconv.type[i] ~ dcat(P.conv[1:2])

    unconv.type[i] <- ifelse(fixed.unconv[i] != 0, fixed.unconv[i], runconv.type[i])
    runconv.type[i] ~ dcat(P.unconv[1:2])

    }

  for (j in 1:nact) {
    alpha[j] ~ dnorm(0, 0.01)
  }

  for (j in 1:(nact-1)) {
    alpha.conv[j] ~ dlnorm(0, 0.01)
  }

  for (j in 2:nact) {
    alpha.unconv[j] ~ dlnorm(0, 0.01)
  }

  P.conv[1:ntypes] ~ ddirch(wprior[1:ntypes])
  P.unconv[1:ntypes] ~ ddirch(wprior[1:ntypes])

  alpha.conv[nact] <- 0
  alpha.unconv[1] <- 0

}	

"
writeLines(modelString, con = "stage1.bug")

##---------- END OF BUGS/JAGS CODE ----------##

data.jags  <- list("n" = n, "wprior" = wprior, "ntypes" = ntypes, "nact" = nact, "Y"=  as.matrix(Y), "fixed.conv" = fixed.conv, "fixed.unconv" = fixed.unconv)

parameters.jags <- c("alpha", "conv.type", "alpha.conv", "P.conv", "unconv.type", "alpha.unconv", "P.unconv")

samples.stage1 <- foreach(j = 1:3) %dopar% {
  
  inits.jags <- list(list(.RNG.seed = sample(1:10000, 1), .RNG.name = "base::Mersenne-Twister"))
  
  model.jags <- jags.model("stage1.bug", data = data.jags, inits = inits.jags, n.chains = 1, n.adapt = 50000)
  
  mcmc.samples.stage1 <- coda.samples(model.jags, variable.names = parameters.jags, n.iter = 50000, thin = 50)
  
  mcmc.samples.stage1[[1]]
  
}

save(samples.stage1, file = "samples_stage1.Rdata")
