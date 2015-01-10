import numpy as np
import bisect
import random
import sys

# pick k numbers less than n
def get_random_set(n, k):
    T = range(n)
    S = []
    for i in range(k):
        j = random.randint(0, len(T)-1)
        bisect.insort(S, T.pop(j)) # insert into sorted list
    return S, T


# in sorted list, take out x and add y
def switch(sl, x, y):
    sl.remove(x)
    bisect.insort(sl, y)


def matsquare(M):
    N = np.asarray(M)
    n, m = M.shape
    for i in range(n):
        for j in range(n):
            for k in range(n):
                N[i,j] = min(N[i,j], M[i, k] + M[k, j])
    return N


def graph_distance(A):
    D = np.asarray(A, dtype=np.double)
    D[D==0] = np.inf
    while True:
        D_ = matsquare(D)
        if (D==D_).all():
            return D
        D = D_

def reachability(A):
    D = graph_distance(A)
    R = np.asarray([1 if d > 0 and d < np.inf else 0 for d in D.flatten()])
    R.shape = D.shape
    return R


# KPP-Neg metric:

# number of disconnected pairs (Borgatti (3))
# minimize sum(i>j) R_ij
# Equation (4) should be more tractable in the optimization context (although lose out on elegance of equation)
def metric_3(A):
    R = reachability(A)
    s = np.sum(np.tril(R, k=-1))
    return s

# Borghatti (9)
# maximize D_F = 1 - 2 * (sum(i>j) (1 / d_ij) / (n * n-1)
# equivalent to minimizing just the sum
def metric_9(A):
    D = graph_distance(A)
    return np.sum(np.tril(1/D, k=-1))


def greedy_optimize(A, k, metric, tolerance):
    S, T = get_random_set(len(A), k) # sorted. S union T = {1 ... |A|}
    
    B = A[T, :][:, T] # A without col i and row i if i in S
    fit = metric(B)

    print "Fit: %s" % (fit)
    while True:
        Dfit = np.inf
        pair = None
        for u in S:
            for v in T:
                T_ = list(T)
                switch(T_, v, u)
                B_ = A[T_, :][:, T_]
                fit_ = metric(B_)
                d = fit - fit_
                if (d >= 0) and (d < Dfit):
                    Dfit = d
                    pair = u,v
        if (Dfit < tolerance):
            break
        u, v = pair
        switch(S, u, v)
        fit = fit - d

    print "New fit: %s" % (fit)
    return S

# random graph adjacency matrix
def rgraph(n, p):
    m = np.zeros(n*n)
    for i in range(n*n):
        if (np.random.rand() > p):
            m[i]=1
    m.shape = (n,n)
    return m

def main(argv):
    n = int(argv[1]) # number of nodes in graph
    p = float(argv[2]) # probability two nodes are connected
    k = int(argv[3]) # number of nodes to find
    tol = float(argv[4]) # tolerance to stop optimize algorithm at
    t = int(argv[5]) # times to repeat the optimize algorithm
    m = int(argv[6]) # metric (3 or 9)
   
    if (m == 3):
        metric = metric_3
    elif (m == 9):
        metric = metric_9
    else:
        print "Invalid metric! Use 3 or 9 only."

    G = rgraph(n, 0.5)

    for i in range(t):
        S = greedy_optimize(G, k, metric, tol)
        print "nodes: %s" % S
        print

if __name__ == "__main__":
    main(sys.argv)



