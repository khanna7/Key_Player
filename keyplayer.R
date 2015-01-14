library(statnet)

# pick k numbers less than n. return S, T, where |S|=k and |T|=n-k, S union T = {1,2,..,n} 
# no longer used, switched to using a mask (see get_random_mask)

get_random_set <- function(G, k) {
    # input G as the graph
    n <- network.size(G)
    range <- 1:n
    samp <- sample(n, k)
    s <- range[range %in% samp]
    t <- range[! range %in% samp]
    return (list(s,t));
}


graph_distance <- function(x) { 
    return(geodist(x)$gdist);
}

reachability = function(A) {
    D = graph_distance(A)
    return ((D > 0) & (D < Inf));
}

# trim array based on True, False mask (see get_random_mask above)
trimmed_array = function(A, s) {
    B <- A[s,][,s]
    return(B);
}

# A *metric* is a function that takes an adjancancy matrix G an indicator vector s, where:
#   s = [s1, s2, s3, .. sn] where s_i = TRUE iff node i in S
# and returns a number between 0 and 1.

# KPP-Neg metric:

# number of disconnected pairs (Borgatti (3))
# maximize D =  1 - 2 * sum(i>j) R_ij / (n * (n-1))
# Equation (4) should be more tractable in the optimization context (although lose out on elegance of equation)
metric_3 <- function(G, s) {
    A <- trimmed_array(G, !s)
    n <- length(A[1,])
    R <- reachability(A);
    S <- sum(R[lower.tri(R)]);
    return(1 - 2*S/(n * (n-1)));
}

# Borghatti (9)
# maximize D_F = 1 - 2 * (sum(i>j) (1 / d_ij) / (n * (n-1))
metric_9 <- function(G, s) {
    A <- trimmed_array(G, !s)
    n <- length(A[1,])
    D <- graph_distance(A);
    Dr <- 1/D
    S <- sum(Dr[lower.tri(D)]);
    return(1 - 2*S/(n * (n-1)));
}

# KPP-Pos metric:

# Borghatti (14)
# d_Kj = min { d(i, j) where i in K, j in G\S }
# D_R = sum(j) (1/d_Kj) / n
metric_14 <- function(G, s) {
    n <- length(G[1,])
    D <- graph_distance(G) # future optimization: only do the distance calculation ONE time.

    H <- apply(D[!s, s], 1, min); # 0 < d <= Inf 
    return (sum(1/H)/n)
}


# return something like [TRUE, FALSE, FALSE, TRUE, ...] where there are n total values and k FALSEs
get_random_mask = function(n, k) {
    s <- sample(n, k);
    s <- 1:n %in% s;
    return (s);
}


greedy_optimize = function(A, k, metric, tolerance) {
    n <- length(A[1,])

    s <- get_random_mask(n, k); # random index matrix to start with.
                                # our "target" nodes are {i | s[i] is TRUE}


    fit <- metric(A, s)

    print(paste("Fit: ", fit))
    i <- 0
    while (TRUE) {
        i <- i + 1
        Dfit = 0
        pair = NULL
        for (u in which(s)) {
            for (v in which(!s)) {
                s_ <- s # clone
                s_[v] = TRUE;
                s_[u] = FALSE;
                
                fit_ = metric(A, s_)
                d = fit_ - fit
                if ((d >= 0) && (d > Dfit)) {
                    Dfit = d
                    pair = list(u,v)
                }
            }
        }
        if (Dfit < tolerance)
            break
        print(paste("Iteration",i,":", fit, "=>",fit+Dfit))
        u <- pair[[1]]
        v <- pair[[2]]
        s[v] = TRUE
        s[u] = FALSE
        fit = fit + Dfit
    }

    print(paste("New fit (", i, " iterations): ", fit));
    return(which(s));
}


main = function(argv) {
    ## n = as.numeric(argv[1]) # number of nodes in graph
    ## p = as.numeric(argv[2]) # probability two nodes are connected
    G = as.network(argv[[1]])
    k = as.numeric(argv[[2]]) # number of nodes to find
    tol = as.numeric(argv[[3]]) # tolerance to stop optimize algorithm at
    t = as.numeric(argv[[4]]) # times to repeat the optimize algorithm
    m = as.numeric(argv[[5]]) # metric (3 or 9)
   
    if (m == 3)
        metric <- metric_3
    else if (m == 9)
        metric <- metric_9
    else if (m == 14)
        metric <- metric_14
    else {
        print("Invalid metric! Use 3, 9, or 14 only.")
        q()
    }

    A <- as.matrix(G)
    for (i in 1:t) {
        S <- greedy_optimize(A, k, metric, tol)
        print(sprintf("nodes: %s", paste(S, collapse=",")))
        cat('\n')
    }
}


argv <- commandArgs(trailingOnly = TRUE);
if (length(argv) > 1)
    main(argv)

