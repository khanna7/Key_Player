###################################################################
### Test keyplayer function in R
###################################################################

  ## Packages and libraries
  library(statnet)
  library(sna)

  ## Source function
  source("keyplayer.R")

  ## Trial datasets
  data(flo)
  flo.net <- as.network(flo, directed=FALSE)
  flo.net  
  gplot(flo.net,
        #displaylabels=TRUE,
        gmode="graph",
        vertex.cex=0.5,
        label=1:network.size(flo.net),
        edge.col="#00000011",
        )
  
  ## Test function
  list.input <- list(rep(0, 5))
  list.input[[1]] <- flo.net #network
  list.input[[2]] <- 5 #cardinality of KP-set
  list.input[[3]] <- 0.1 #tolerance
  list.input[[4]] <- 560 #number of iterations (or starting sets)
  list.input[[5]] <- 3 #metric

  main(list.input)
  
  list.input[[5]] <- 9 # metric
  main(list.input)

###################################################################
