
weak.info <- function(info, data, stdize=FALSE) {
## When binary treatment and/or outcomes are involved, we
## follow the weakly informative advice found in Gelman et al. 2008
## http://www.stat.columbia.edu/~gelman/research/published/priors11.pdf

    N <- nrow(data$X)
    p <- ncol(data$X)
    q <- ncol(data$Z)
    r <- p+q
    s <- p+r

    if(stdize) {
        X <- standard(data$X)
        data$X <- X[[1]]

        Z <- standard(data$Z)
        data$Z <- Z[[1]]

        sqrtT <- diag(0.5*c(X$bt[2, ], Z$bt[2, ], X$bt[2, ]), s, s)
        mu1 <- sqrtT %*% c(X$bt[1, ], Z$bt[1, ], rep(0, p))
        mu2 <- sqrtT %*% c(rep(0, r), X$bt[1, ])

        bt <- list(sqrtT=sqrtT, mu1=mu1, mu2=mu2)
    }
    else bt <- NULL

    if(length(info)==0) {
        ## default prior setup
        if(length(data$tbin)>0 & length(data$ybin)==0) {
            ## binary treatment and numeric outcome
            beta.prec  <- 0.001
            theta.prec <- diag(c(rep(0.04, r), rep(0.001, p)), s, s)
            mu.prec    <- diag(c(1, 0.001), 2, 2)
        }
        else if(length(data$tbin)>0 & length(data$ybin)>0) {
            ## binary treatment and binary outcome
            beta.prec  <- 1
            theta.prec <- diag(0.04, s, s) 
            mu.prec    <- diag(1, 2, 2)
        }
        else {
            ## numeric treatment and numeric outcome
            ## not weakly informative since neither are binary
            beta.prec  <- 0.001
            theta.prec <- diag(0.001, s, s) 
            mu.prec    <- diag(0.001, 2, 2)
        }

        info <- list(theta=list(init=rep(0, s),
                         prior=list(mean=rep(0, s), prec=theta.prec)),
                     beta=list(init=0, prior=list(mean=0, prec=beta.prec)),
                     mu=list(init=c(0, 0),
                         prior=list(mean=c(0, 0), prec=mu.prec)),
                     rho=list(init=0),
                     tau=list(init=1),
                     dpm=list(m=as.integer(3),
                         alpha=list(fixed=as.integer(0), init=1,
                             prior=list(a=3, b=4)),
                         C=as.integer(0*(1:N)), states=as.integer(N),
                         prior=list(mu0=c(0, 0), T0=mu.prec, S0=solve(mu.prec),
                             alpha0=0.5, lambda0=0.4)))
    }
    
    return(list(info=info, data=data, bt=bt))
}
