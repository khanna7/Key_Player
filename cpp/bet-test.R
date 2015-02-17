# ignoreeval is TRUE

library(sna)
library(network)
library(foreach)
dyn.load("bet.so")


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
    n <- attr(dflo, "n")
    delta <- n / p
    c <- rep(0, n)
    for (i in 1:p) {
        st <- as.integer((i-1) * delta)
        end <- as.integer(st + delta)
        c <- c + .Call("betweenness_partial", dat, n, NROW(dflo), 0, st, end)
    }
    return(c)
}


btn_wrapper_apply <- function(dat, p) {
    dat <- as.edgelist.sna(dat)
    n <- attr(dflo, "n")
    delta <- n / p
    f <- function(i) {
        st <- as.integer((i-1) * delta)
        end <- as.integer(st + delta)
        x <- .Call("betweenness_partial", dat, n, NROW(dflo), 0, st, end)
        return(x)
    }
    c <- foreach(i=1:p, .combine='+') %do% f(i)
    return(c)
}
 


data(flo)

for (p in c(1,2,4,8,16)) {
    b <- btn_wrapper(flo, p)
    asserteq(a, b, paste("wrapper ",p))
}

for (p in c(1,2,4,8,16)) {
    b <- btn_wrapper_apply(flo, p)
    asserteq(a, b, paste("wrapper-sapply ",p))
}
