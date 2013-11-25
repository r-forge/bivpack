\name{BIVpack-package}
\alias{BIVpack-package}
\alias{BIVpack}
\docType{package}
\title{
Bayesian methods for IV regression
}
\description{
An R package to facilitate Instrumental Variables regression with Bayesian nonparametric methods.
}
\details{
\tabular{ll}{
Package: \tab BIVpack\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-21\cr
License: \tab GPL 2 or later\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Rodney Sparapani

Maintainer: Rodney Sparapani <rsparapa@mcw.edu>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
\code{\link[Rcpp:Rcpp-package]{Rcpp}} 
\code{\link[RcppEigen:RcppEigen-package]{RcppEigen}} 
}
\examples{

require(BIVpack)

N <- 10
p <- 0
q <- 1
p1 <- max(1, p)
r <- p1+q
s <- p1+r

gamma <- 0
delta <- 4 
eta  <- 0
beta <- 0.5
mu   <- 0
rho  <- 0.6

mcmc <- list(M=1, burnin=0, thin=1)

info <- list(theta=list(init=c(rep(gamma, p1), rep(delta, q), rep(eta, p1)),
               prior=list(mean=rep(0., s),
                 prec=diag(0.001, s))),
             beta=list(init=beta, prior=list(mean=0., prec=0.001)),
             Tprec=list(init=solve(matrix(c(1, rho, rho, 1), 2, 2)),
               prior=list(nu=4, Psi=diag(1, 2, 2))),
             mu=list(init=c(mu, mu),
               prior=list(mean=c(0., 0.), prec=diag(0.001, 2))))

sim <-  list(X=matrix(0, nrow=N, ncol=p1),
             Z=matrix(c(-0.24146164, -0.29673723, -0.27538621,
               0.41463628, 0.39023100, -0.22045922, -0.07062149,
               -0.22595298, 0.01247487, -0.14472589), nrow=N, ncol=q),
             t=c(-2.01322819, -2.04167660, -0.56128516,
               0.20783192, 0.31477076, -1.41477107, -0.38701899,
               -0.59955150, 0.01197733, -0.79804809),
             y=c(-1.9924944, -1.9345279, -1.3781082, -0.7646928,
               -0.2881649, 0.1545577, -0.6114224, -0.3703420,
               0.2320320, 0.7451867))

set.seed(42)

nniv(info, sim, mcmc)
##should produce approx...
##          [,1]     [,2]     [,3]      [,4]       [,5]        [,6]     [,7]
## [1,] -17.85732 4.008851 20.01287 0.8217393 -0.3821324 -0.05601638 3.108135
##           [,8]      [,9]    [,10]
## [1,] -1.306804 -1.306804 2.188039

info <- list(beta=list(init=0.,
               prior=list(mean=0., prec=0.001)),
             rho=list(init=0.6),
             mu=list(init=c(0.,0.),
               prior=list(mean=c(0.,0.), prec=diag(0.001,2))),
             tau=list(init=1., prior=list(alpha0=0.1, lambda0=0.1)),
             theta=list(init=c(rep(0., p1), rep(0., q), rep(0., p1)),
               prior=list(mean=rep(0.,s), prec=diag(0.001, s))))

sim <-  list(X=matrix(0, nrow=N, ncol=p1),
             Z=matrix(c(-0.24146164, -0.29673723, -0.27538621,
               0.41463628, 0.39023100, -0.22045922, -0.07062149,
               -0.22595298, 0.01247487, -0.14472589), nrow=N, ncol=q),
             tbin=as.integer(c(0, 0, 0, 1, 1, 0, 0, 0, 1, 0)),
             y=c(-1.9924944, -1.9345279, -1.3781082, -0.7646928,
               -0.2881649, 0.1545577, -0.6114224, -0.3703420,
               0.2320320, 0.7451867))

set.seed(42)

(par.post <- bniv(info, sim, mcmc))
##should produce approx...
##      [,1]     [,2]     [,3]       [,4]      [,5]       [,6]      [,7]     [,8]
## [1,] 47.79852 1.231516 63.82816 -0.1581546 0.1076077 -0.2405642 0.4071552 0.626958

bnivE(par.post[1, ], sim$Z, sim$X)

info <- list(beta=list(init=0, prior=list(mean=0, prec=0.001)),
             rho=list(init=0),
             mu=list(init=c(0, 0)),
             tau=list(init=1),
             theta=list(init=rep(0, s),
               prior=list(mean=rep(0, s),
                   prec=diag(c(rep(0.04, r), rep(0.001, p1)), s, s))),
             dpm=list(m=as.integer(3),
                 alpha=list(fixed=as.integer(0), init=1, prior=list(a=3, b=4)),
                 C=as.integer(0*(1:N)), states=as.integer(N),
                 prior=list(mu0=c(0, 0), T0=diag(0.001, 2), S0=diag(1000, 2),
                    alpha0=0.5, lambda0=0.4)))

set.seed(42)

(non.post <- bnivDPM(info, sim, mcmc))
##should produce approx...
## [[1]]
## [[1]]$beta
## [1] -0.1581546

## [[1]]$theta
## [1] 47.798522  1.231516 63.828162

## [[1]]$C
##  [1] 0 0 0 0 0 0 0 0 0 0

## [[1]]$phi
##            [,1]       [,2]       [,3]      [,4]
## [1,] -0.3704315 -0.9304674 0.08751207 0.4762594

## [[1]]$states
## [1] 10

set.seed(42)

(non.post <- bnivDPM(NULL, sim, mcmc))

bnivE(non.post[[1]], sim$Z, sim$X)

info <- list(theta=list(init=c(rep(gamma, p1), rep(delta, q), rep(eta, p1)),
               prior=list(mean=rep(0., s),
                 prec=diag(0.001, s))),
             beta=list(init=beta, prior=list(mean=0., prec=0.001)),
             rho=list(init=rho),
             mu=list(init=c(mu, mu),
               prior=list(mean=c(0., 0.), prec=diag(0.001, 2))))

sim <-  list(X=matrix(0, nrow=N, ncol=p1),
             Z=matrix(c(-0.24146164, -0.29673723, -0.27538621,
               0.41463628, 0.39023100, -0.22045922, -0.07062149,
               -0.22595298, 0.01247487, -0.14472589), nrow=N, ncol=q),
             tbin=as.integer(c(0, 0, 0, 1, 1, 0, 0, 0, 1, 0)),
             ybin=as.integer(c(0, 0, 0, 0, 0, 1, 0, 0, 1, 1)))

set.seed(42)

(par.post <- bbiv(info, sim, mcmc))
##should produce approx...
##        [,1]     [,2]      [,3]      [,4]       [,5]       [,6]      [,7]
## [1,] 3.41756 3.697175 -15.94475 0.1152464 -0.7233647 -0.8709668 0.3693283

bbivE(par.post[1, ], sim$Z, sim$X)

info <- list(beta=list(init=0, prior=list(mean=0, prec=1)),
             rho=list(init=0),
             mu=list(init=c(0, 0)),
             tau=list(init=1),
             theta=list(init=rep(0, s),
               prior=list(mean=rep(0, s), prec=diag(0.04, s))),
             dpm=list(m=as.integer(3),
                 alpha=list(fixed=as.integer(0), init=1, prior=list(a=3, b=4)),
                 C=as.integer(0*(1:N)), states=as.integer(N),
               prior=list(mu0=c(0, 0), T0=diag(1, 2), S0=diag(1, 2),
                    alpha0=0.5, lambda0=0.4)))

set.seed(42)

(non.post <- bbivDPM(info, sim, mcmc))
##should produce approx...
## [[1]]
## [[1]]$beta
## [1] -0.2440909

## [[1]]$theta
## [1]  0.5403636  1.4322290 -2.5210857

## [[1]]$C
##  [1] 0 0 0 0 1 0 1 0 0 0

## [[1]]$phi
##           [,1]       [,2]  [,3]
## [1,] 0.5215824 -0.9161644 -0.39
## [2,] 0.1711343 -0.9186419 -0.43

## [[1]]$states
## [1] 8 2

## [[1]]$alpha
## [1] 1.029461

set.seed(42)

(non.post <- bbivDPM(NULL, sim, mcmc))

bbivE(non.post[[1]], sim$Z, sim$X)

set.seed(42)

bbivDPM(info, sim, mcmc, stdize=TRUE)
##should produce approx...
## [[1]]
## [[1]]$beta
## [1] 0.09688739

## [[1]]$theta
##           [,1]
## [1,]  1.708780
## [2,]  6.355211
## [3,] -7.972373

## [[1]]$C
##  [1] 0 0 0 0 0 0 0 0 0 0

## [[1]]$phi
##            [,1]       [,2] [,3]
## [1,] -0.1740396 -0.6591617 0.47

## [[1]]$states
## [1] 10
}
