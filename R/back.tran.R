
back.tran <- function(post, bt) {
## When binary treatment and/or outcomes are involved, we
## follow the weakly informative advice found in Gelman et al. 2008
## http://www.stat.columbia.edu/~gelman/research/published/priors11.pdf

    for(i in 1:length(post)) {

        bt.mu1 <- sum(bt$mu1 * post[[i]]$theta)
        bt.mu2 <- sum(bt$mu2 * post[[i]]$theta)
        
        for(k in 1:length(post[[i]]$states)) {
            post[[i]]$phi[k, 1] <- post[[i]]$phi[k, 1]-bt.mu1
            post[[i]]$phi[k, 2] <- post[[i]]$phi[k, 2]-bt.mu2
        }
        
        post[[i]]$theta <- c(bt$sqrtT %*% post[[i]]$theta)
    }

    return(post)
}
