# A QUEUE is dequeued from the front, enqueued at the back
# A STACK is pushed from the front, popped at the front
# Both are actually just vectors.

# A is an adjacency matrix, n is the number of nodes in the graph
mybetweenness <- function(A) {
    n <- network.size(A)
    C <- rep(0, n)
    for (s in 1:n) {
        S <- c()
        P <- vector("list", n)
        sig <- rep(0, n)
        sig[s] = 1
        d <- rep(-1, n)
        d[s] = 0
        Q = c(s) # this is a QUEUE
        while (length(Q) > 0) {
            v <- Q[1]
            Q <- Q[-1]
            S <- c(v, S)
            
            neighbors <- get.neighborhood(A, v) #which(A[v,] > 0)
            for (w in neighbors) {
                if (d[w] < 0) {
                    Q <- c(Q, w)
                    d[w] <- d[v] + 1
                }
                if (d[w] == d[v] + 1) {
                    sig[w] <- sig[w] + sig[v]
                    P[[w]] <- c(P[[w]], v)
                }
            }
        }

        delta <- rep(0, n)
        while (length(S) > 0) {
            w <- S[1]
            S <- S[-1]
            for (v in P[[w]]) {
                delta[v] <- delta[v] + ((sig[v]/sig[w])*(1+delta[w]))
            }
            if (w != s) {
                C[w] <- C[w] + delta[w]
            }
        }
    }
    return(C)
}
