
standard <- function(A) {
## When binary treatment and/or outcomes are involved, we
## follow the weakly informative advice found in Gelman et al. 2008
## http://www.stat.columbia.edu/~gelman/research/published/priors11.pdf

    p <- ncol(A)

    bt <- matrix(c(rep(0, p), rep(1, p)), byrow=TRUE,
                   nrow=2, ncol=p, dimnames=list(c("mean", "std"), NULL))

    for(j in 1:p) {
        bt[1, j] <- mean(A[ , j])
        bt[2, j] <- sd(A[ , j])

        if(bt[2, j] == 0) bt[2, j] <- 1
        else bt[2, j] <- 1/bt[2, j]
        
        A[ , j] <- 0.5*(A[ , j]-bt[1, j])*bt[2, j]
    }

    return(list(A, bt=bt))
}
