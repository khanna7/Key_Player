# ignoreeval is TRUE

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
    cl <- lapply(1:p, f)
    rowSums(data.frame(cl))
    #c <- rep(0, n)
    #for (i in 1:p) {
    #    c <- c + cl[[i]]
    #}
    #return(c)
}
 

library(sna)
library(network)
dyn.load("bet.so")

data(flo)
dflo = as.edgelist.sna(flo)
n <- attr(dflo, "n")

a <- betweenness(flo)
b <- .Call("betweenness_partial", dflo, n, NROW(dflo), 0, as.integer(0), as.integer(n))

asserteq(a, b, "btn ccall")

c <- rep(0, n)
for (i in c(0,4,8,12)) {
    c <- c + .Call("betweenness_partial", dflo, n, NROW(dflo), 0, as.integer(i), as.integer(i+4));
}

asserteq(a, c, "btn ccall 4")

for (p in c(1,2,4,8,16)) {
    d <- btn_wrapper(flo, p)
    asserteq(a, d, paste("wrapper ",p))
}

for (p in c(1,2,4,8,16)) {
    d <- btn_wrapper_apply(flo, p)
    asserteq(a, d, paste("wrapper-sapply ",p))
}
