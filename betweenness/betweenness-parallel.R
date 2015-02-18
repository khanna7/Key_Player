# "betweenness" from sna, all default options.
mybetweenness = function (dat, g = 1, nodes = NULL, gmode = "digraph", diag = FALSE, 
    tmaxdev = FALSE, cmode = "directed", geodist.precomp = NULL)
{
    dat <- as.edgelist.sna(dat)
    n <- attr(dat, "n")
    
 
    if (is.null(nodes)) 
        nodes <- 1:n
    if (cmode == "undirected") 
        dat <- symmetrize(dat, rule = "weak", return.as.edgelist = TRUE)
    
    meas <- switch(cmode, undirected = 0, directed = 0, endpoints = 1, 
        proximalsrc = 2, proximaltar = 3, proximalsum = 4, 
        lengthscaled = 5, linearscaled = 6)

    bet <- .Call("betweenness_partial", dat, n, NROW(dat), meas)
    
    if ((cmode == "undirected") || (gmode == "graph")) 
        bet <- bet/2
 
    bet <- bet[nodes]
    return(bet)
}

# brn_wrapper_dopar from the tests
betweenness.parallel <- function(dat, p) {
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
