#!/usr/bin/env Rscript

# Do not run this test on a login node.

# Make sure snow is loaded after Rmpi
library(sna)
library(foreach)
library(Rmpi)
library(snow)
library(doParallel)

source("../betweenness-parallel.R")

load('fb2.rdata')

# Set up the cluster.
# The path to bet.so is hardcoded, so it'll need to change if you move the file.
np <- mpi.universe.size() - 1 # number of cores = total cores available - 1 (we're already using ONE core)
cl <- makeMPIcluster(np)
# The process running on each core needs to load the function into its environment:
i <- clusterCall(cl, function() { dyn.load("../lib/bet.so") }) 
registerDoParallel(cl)


fb.deg2.btwn <- betweenness.parallel(fb.net.deg.gr.2, np)

save(fb.deg2.btwn, file="fb2.rdata")

# Cleanup
stopCluster(cl)
mpi.exit()

