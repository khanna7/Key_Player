#!/usr/bin/env Rscript

# Do not run this test on a login node.

library(sna)
library(network)
#library(ergm) # sampson data set
library(foreach)
library(Rmpi)
library(snow)
library(doParallel)

dyn.load("../lib/bet.so")

asserteq = function(a, b, name) {
    if (sum(abs(a-b)) > 0.00001) {
        print(paste("Not equal! ", name))
        print(a)
        print(b)
        q()
    }
}


btn_wrapper_dopar <- function(dat, p) {
    dat <- as.edgelist.sna(dat)
    n <- attr(dat, "n")
 
    if (n %% p != 0) {
        m <- n + (p - n%%p)
        delta <- m / p
    }
    else
        delta <- n / p

    c <- foreach(i=1:p, .combine='+', .inorder=FALSE) %dopar%
    {
        st <- as.integer((i-1) * delta)
        end <- min(as.integer(st + delta), n)
        .Call("betweenness_partial", dat, n, NROW(dat), 0, st, end)
    }

    return(c)
}
 

data(flo)

a <- betweenness(flo)

np <- mpi.universe.size() - 1
cl <- makeMPIcluster(np)
i <- clusterCall(cl, function() { dyn.load("../lib/bet.so") })

registerDoParallel(cl)

for (p in c(1,2,4,8,16)) {
    b <- btn_wrapper_dopar(flo, p)
    asserteq(a, b, paste("wrapper-dopar ",p))
}

stopCluster(cl)
mpi.exit()

