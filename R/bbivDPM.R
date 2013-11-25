
bbivDPM <- function(info, data, mcmc, stdize=FALSE) {
    temp1 <- weak.info(info, data, stdize)

    temp2 <- .Call( "bbivDPM", temp1$info, temp1$data, mcmc, PACKAGE = "BIVpack" )

    if(stdize) temp2 <- back.tran(temp2, temp1$bt)

    return(temp2)
}
