###################################################################
### Test keyplayer function in R
###################################################################

  ## Packages and libraries
  library(statnet)
  library(sna)

  ## Source function
  source("keyplayer.R")

  ## Trial datasets
  gdata(flo)
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
  list.input[[1]] <- flo.net
  list.input[[2]] <- 3
  list.input[[3]] <- 0.001
  list.input[[4]] <- 10
  list.input[[5]] <- 9

  main(list.input)
  
  list.input[[5]] <- 9
  main(list.input)

###################################################################
