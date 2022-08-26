using LightGraphs
using GraphIO
using LinearAlgebra
using Random
import Laplacians

function readGraph(gn; dataPath = "../data/") 
    G = SimpleGraph(loadgraph(dataPath*gn*".txt", gn, EdgeListFormat()))
    for i = 1 : nv(G)
        if has_edge(G, i, i); rem_edge!(G, i, i); end
    end
    return G
end

function getDegreeProduct(G)
    n,m = nv(G),ne(G)
    d = degree(G)
    ret = zeros(m)
    for (i,x) in enumerate(collect(edges(G)))
        (u, v) = Tuple(x)
        ret[i] = d[u]*d[v]
    end
    return ret
end

function getDegreeSum(G)
    n,m = nv(G),ne(G)
    d = degree(G)
    ret = zeros(m)
    for (i,x) in enumerate(collect(edges(G)))
        (u, v) = Tuple(x)
        ret[i] = d[u]+d[v]
    end
    return ret
end

function getApproxMyCent(G; eps=0.3)
    n,m = nv(G),ne(G)
    f = Laplacians.approxchol_sddm(laplacian_matrix(G, Float64)+I(n), tol=1e-9)
    B = incidence_matrix(G, Float64, oriented=true)
    k = round(Int, log2(n)/eps^2)
    rho,rho2 = zeros(m),zeros(m)
    for i = 1 : k
        x1,x2 = randn(n),randn(m)
        y1,y2 = f(x1),f(B*x2)
        for (i,x) in enumerate(collect(edges(G)))
            (u, v) = Tuple(x)
            rho[i] += (y1[u]-y1[v])^2 + (y2[u]-y2[v])^2
            rho2[i] += (y1[u]-y1[v])^2
        end
    end
    rho ./= k
    rho2 ./= k
    ret = zeros(m)
    foreach(i -> ret[i] = n*rho2[i]/(1.0-rho[i]), 1:m)
    return ret
end

function getMyCent(G; W=nothing, W2=nothing)
    n,m = nv(G),ne(G)
    if W==nothing; W = inv(laplacian_matrix(G)+I(n)+zeros(n,n)); end
    if W2==nothing; W2 = W * W; end
    ret = zeros(m)
    for (i,x) in enumerate(collect(edges(G)))
        (u, v) = Tuple(x)
        rho = W[u,u]+W[v,v]-2*W[u,v]
        rho2 = W2[u,u]+W2[v,v]-2*W2[u,v]
        ret[i] = n*rho2/(1.0-rho)
    end
    return ret
end

function getEB(G)
    n, m = nv(G), ne(G)
    idx = zeros(Int32, n, n)
    g = Array{Array{Int32, 1}, 1}(undef, n)
    foreach(i -> g[i] = [], 1 : n)
    ID = 0
    for x in collect(edges(G))
        (u, v) = Tuple(x)
        ID += 1
        idx[u, v] = ID
        idx[v, u] = ID
        push!(g[u], v)
        push!(g[v], u)
    end
    Cb = zeros(m)
    p = Array{Array{Int32, 1}, 1}(undef, n)
    d = zeros(Int32, n)
    S = zeros(Int32, n+10)
    sigma = zeros(n)
    d = zeros(Int32, n)
    Q = zeros(Int32, n+10)
    delta = zeros(n)

    for s = 1 : n
        foreach(i -> p[i] = [], 1 : n)
        top = 0
        sigma .= 0
        sigma[s] = 1.0
        d .= -1
        d[s] = 0
        front = 1
        rear = 1
        Q[1] = s

        while front <= rear
            v = Q[front]
            front += 1
            top += 1
            S[top] = v
            for w in g[v]
                if d[w] < 0
                    rear += 1
                    Q[rear] = w
                    d[w] = d[v] + 1
                end
                if d[w] == (d[v] + 1)
                    sigma[w] += sigma[v]
                    push!(p[w], v)
                end
            end
        end

        delta .= 0

        while top > 0
            w = S[top]
            top -= 1
            for v in p[w]
                delta[v] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
                Cb[idx[v, w]] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
            end
        end

    end

    return Cb
end

function getScore(G)
    n = nv(G)
    W = inv(laplacian_matrix(G)+I(n)+zeros(n,n))
    return n*tr(W)-n
end

function getDelta(G0, ans)
    s0 = getScore(G0)
    G = deepcopy(G0)
    for (u,v) in ans
        @assert rem_edge!(G, u, v)
    end
    s1 = getScore(G)
    return s1-s0
end