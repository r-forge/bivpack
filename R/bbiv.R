
bbiv <- function(info, data, mcmc){
    temp <- .Call( "bbiv", info, data, mcmc, PACKAGE = "BIVpack" )

    p <- ncol(data$X)
    q <- ncol(data$Z)

    dimnames(temp) <- list(NULL, c(paste0("gamma", 1:p),
                                    paste0("delta", 1:q),
                                    paste0("eta", 1:p),
                                    "beta", "mu1", "mu2", "rho"))

    return(temp)
}
