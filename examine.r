##############################################################################
# Metropolis-Hastings MCMC
#
# Runs a Metropolis-Hasting MCMC chain for a given likelihood function.
# Proposal steps are sampled from a Gaussian distribution, either in a single
# step or sequentially over the parameter space.
#
# Input:
#   theta      : starting value of the chain
#   lik.fun    : likelihood function
#   prior.fun  : prior probability function
#   ...        : parameters passed to lik.fun and prior.fun
#   V          : variance-covariance matrix for the proposal
#   maxit      : length of the chain
#   sequential : TRUE or FALSE whether proposals are chosen for all parameters
#                in a single step or sequentially.
#
# Output:
#   list()
#     chain    : matrix with the posterior samples
#     lik      : likelihood values of those posterior samples
##############################################################################

metrop <- function(theta,lik.fun,prior.fun,
                   V=diag(theta),maxit=1000,thin=1,
                   sequential=FALSE,pb=FALSE,...) {
  if (! sequential) {
    if (! require(mvtnorm)) stop("You need to install the mvtnorm library.")
  }
  if (pb) { .pb <- txtProgressBar(0,maxit,style=3) }
  ndim <- length(theta)
  chain <- matrix(NA,nrow=maxit,ncol=ndim)
  lik <- vector(length=maxit)
  last <- theta
  chain[1,] <- theta
  last.lik <- lik.fun(theta,...)
  lik[1] <- last.lik
  last.prior <- prior.fun(theta,...)
  it <- 1
  naccept <- 0
  if (is.finite(last.lik)) {
    for (it in seq(2,maxit)) {
      if (pb) setTxtProgressBar(.pb,it)
      if (sequential) {
        for (j in 1:thin) {
          for (i in 1:length(theta)) {
            if (V[i,i] > 0) {
              accept <- FALSE
              proposal <- last
              proposal[i] <- rnorm(1,mean=last[i],sd=sqrt(V[i,i]))
              proposal.prior <- prior.fun(proposal,...)
              if (is.finite(proposal.prior)) {
                proposal.lik <- lik.fun(proposal,...)
                alpha <- exp(proposal.lik+proposal.prior-last.lik-last.prior)
                if (alpha > runif(1)) accept <- TRUE
              }
              if (accept) {
                last <- proposal
                last.lik <- proposal.lik
                last.prior <- proposal.prior
                naccept <- naccept + 1/sum(diag(V) > 0)
              }
            }
          }
        }
        chain[it,] <- last
        lik[it] <- last.lik
      } else {
        for (j in 1:thin) {
          accept <- FALSE
          proposal <- rmvnorm(1,mean=last,sigma=V)
          proposal.prior <- prior.fun(proposal,...)
          if (is.finite(proposal.prior)) {
            proposal.lik <- lik.fun(proposal,...)
            alpha <- exp(proposal.lik+proposal.prior-last.lik-last.prior)
            if (alpha > runif(1)) accept <- TRUE
          }
          if (accept) {
            last <- proposal
            last.lik <- proposal.lik
            last.prior <- proposal.prior
            naccept <- naccept+1
          }
        }
        chain[it,] <- last
        lik[it] <- last.lik
        if (! pb) {
          message(sprintf("Acceptance ratio = %.4f",naccept/(it*thin)))
        }
      }
    }
  }
  if (pb) close(.pb)
  message(sprintf("Acceptance ratio = %.4f",naccept/(maxit*thin)))
  return(list(chain=chain,lik=lik))
}
 2-dic.R
##############################################################################
# Estimate Deviance Information Criterion (DIC)
#
# References:
#   Bayesian Data Analysis.
#   Gelman, A., Carlin, J., Stern, H., and Rubin D.
#   Second Edition, 2003
#
#   Bayesian predictive information criterion for the evaluation of 
#     hierarchical Bayesian and empirical Bayes models. 
#   Ando, T.
#   Biometrika, 2007 
#
# Input:
#   x       : matrix of posterior samples
#   lik     : vector of the likelihood of the posterior samples
#   lik.fun : function that calculates the likelihood
#   ...     : other parameters that are passed to 'lik.fun'
#
# Output:
#   list()
#     DIC   : Deviance Information Criterion
#     IC    : Bayesian Predictive Information Criterion
#     pD    : Effective number of parameters (pD = Dbar - Dhat)
#     pV    : Effective number of parameters (pV = var(D)/2)
#     Dbar  : Expected value of the deviance over the posterior
#     Dhat  : Deviance at the mean posterior estimate
##############################################################################

calc.dic <- function(x,lik,lik.fun,...) {
  D.bar <- -2*mean(lik)
  theta.bar <- summary(x)$statistics[,"Mean"]
  D.hat <- -2*lik.fun(theta.bar,...)
  pD <- D.bar - D.hat
  pV <- var(-2*lik)/2
  list(DIC=pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
}
 3-marginal.likelihood.R
##############################################################################
# Estimate Marginal Likelihood
#
# Reference:
#   Marginal likelihood from the Metropolis--Hastings output.
#   Chib, S. and Jeliazkov, I.,
#   Journal of the American Statistical Association, 2001
# 
# Input:
#   x             : mcmc object or mcmc.list object, 
#                   the samples from the posterior
#   lik           : vector, likelihood values for the samples in x
#   lik.fun       : function, calculates the likelihood
#   prior.fun     : function, calculates the prior probability
#   num.samples   : integer, number of samples used in the estimate
#   log           : boolean, if log-likelihood or likelihood is used
#
# Output:
#   list()
#     ln.m          : marginal liklihood (log if log-liklihood is used)
#     ln.lik.star   : reference likelihood (log if log-liklihood is used)
#                     The reference sample from the posterior is the
#                     sample that has the largest likelihood.
#     ln.pi.star    : refernce prior probability
#                     (log if log-likelihood is used)
#     ln.pi.hat     : posterior ordinate at reference sample
#                     (log if log-likelihood is used)
##############################################################################

marginal.likelihood <- function(x,lik,V,lik.fun,prior.fun,...,
                               num.samples=1000,log=TRUE) {
  if (class(x) != "mcmc" & class(x) != "mcmc.list") {
    stop("x must be an mcmc or mcmcList object.")
  }
  # get mean parameters
  y <- summary(x)
  theta.star <- y$statistics[,"Mean"]
  lik.star <- lik.fun(theta.star,...)
  # get samples from posterior
  g <- sample(1:nrow(x),num.samples,replace=TRUE)
  theta.g <- x[g,]
  q.g <- dmvnorm(theta.g,mean=theta.star,sigma=V,log=FALSE)
  #lik.g <- apply(theta.g,1,lik.fun,...)
  lik.g <- lik[g]
  alpha.g <- sapply(lik.g,function(l) min(1,exp(lik.star-l)))
  # get samples from proposal
  theta.j <- rmvnorm(num.samples,mean=theta.star,sigma=V)
  lik.j <- apply(theta.j,1,lik.fun,...)
  alpha.j <- sapply(lik.j,function(l) min(1,exp(l-lik.star)))
  pi.hat <- mean(alpha.g*q.g)/mean(alpha.j)
  pi.star <- 0
  if (!is.null(prior.fun)) pi.star <- prior.fun(theta.star,...)
  ln.m <- lik.star + pi.star - log(pi.hat)
  if (! log) ln.m <- exp(ln.m)
  list(ln.m=ln.m,ln.lik.star=lik.star,ln.pi.star=pi.star,ln.pi.hat=log(pi.hat))
}
