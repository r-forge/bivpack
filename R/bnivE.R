bnivE <- function(
    one.step,  # one Gibbs step of the posterior
    z,         # a matrix (n by q) of instrumental variables
    x=NULL     # a matrix (n by p) of confounders, can be omitted if p=0
    ) {

    if(length(x)==0) p <- 0
    else p <- ncol(x)
    
    p1 <- max(1, p)
    q <- ncol(z)
    r <- p1+q
    s <- p1+r

    if(class(one.step)=="list") {
        theta <- one.step$theta
        beta  <- one.step$beta
        phi   <- matrix(one.step$phi[ , 1:3], ncol=3)
        states <- one.step$states
    }
    else {
        theta <- one.step[1:s]
        beta  <- one.step[s+1]
        phi   <- matrix(one.step[(s+2):(s+4)], 1, 3)
        states <- nrow(z)
    }

    delta <- matrix(c(theta[(p1+1):(p1+q)]),q,1)
    v <- c(z%*%delta)
    meanv <- mean(v)

    if(p == 0 | length(x)==0){
        vplus <- outer(phi[,1],v,FUN="+") # phi has 3 cols: mu1,mu2,rho
        xeta <- rep(0,length(v))
    }
    else{
        gamma <- theta[1:p1]
        eta <- theta[(p1+q+1):(p1+q+p1)]
        gamma <- matrix(c(gamma),p1,1)
        eta <- matrix(c(eta),p1,1)
        xgamma <- c(x%*%gamma)
        vplus <- outer(phi[,1],v+xgamma,FUN="+") # phi has 3 cols: mu1,mu2,rho
        xeta <- c(x%*%eta)
    }

    eTvx <- apply(pnorm(vplus)*states,2,sum)/sum(states)
    covTV <- mean(v*eTvx)-meanv*mean(eTvx)
    covYV <- beta*(mean(v^2)-meanv^2)

    IVE <- covYV/covTV
    
    return(IVE)
}
