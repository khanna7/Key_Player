#!/usr/bin/env Rscript

# Do not run this test on a login node.

library(sna)
library(network)
library(foreach)
library(Rmpi)
library(snow)
library(doParallel)

source("../betweenness-parallel.R")

load('fb2.rdata')
dat <- fb.net.deg.gr.2


np <- mpi.universe.size() - 1
cl <- makeMPIcluster(np)
i <- clusterCall(cl, load.custom.lib)

registerDoParallel(cl)

fb.deg2.btwn <- betweenness.parallel(dat, np)
save(fb.deg2.btwn, file="fb-between.rdata")

stopCluster(cl)
mpi.exit()

