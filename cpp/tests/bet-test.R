# ignoreeval is TRUE

library(sna)
library(network)
library(ergm) # sampson data set
library(foreach)

dyn.load("../lib/bet.so")


asserteq = function(a, b, name) {
    if (sum(abs(a-b)) > 0.00001) {
        print(paste("Not equal! ", name))
        print(a)
        print(b)
        q()
    }
}

# do p calls to internal betweenness function
btn_wrapper <- function(dat, p) {
    dat <- as.edgelist.sna(dat)
    n <- attr(dat, "n")

    if (n %% p != 0) {
        m <- n + (p - n%%p)
        delta <- m / p
    }
    else
        delta <- n / p
    
    c <- rep(0, n)
    for (i in 1:p) {
        st <- as.integer((i-1) * delta)
        end <- min(as.integer(st + delta), n)
        c <- c + .Call("betweenness_partial", dat, n, NROW(dat), 0, st, end)
    }
    return(c)
}


btn_wrapper_foreach <- function(dat, p) {
    dat <- as.edgelist.sna(dat)
    n <- attr(dat, "n")
 
    if (n %% p != 0) {
        m <- n + (p - n%%p)
        delta <- m / p
    }
    else
        delta <- n / p

    c <- foreach(i=1:p, .combine='+', .inorder=FALSE) %do%
    {
        st <- as.integer((i-1) * delta)
        end <- min(as.integer(st + delta), n)
        .Call("betweenness_partial", dat, n, NROW(dat), 0, st, end)
    }

    return(c)
}
 


data(sampson)

a <- betweenness(samplike)

for (p in c(1,2,4,8,16)) {
    b <- btn_wrapper(samplike, p)
    asserteq(a, b, paste("wrapper ",p))
}

for (p in c(1,2,4,8,16)) {
    b <- btn_wrapper_foreach(samplike, p)
    asserteq(a, b, paste("wrapper-foreach ",p))
}

