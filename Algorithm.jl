include("Util.jl")

using Random
using LinearAlgebra

function Greedy(G0, k; method="Betweenness")
    F = Dict("Betweenness"=>getEB, 
             "DegreeProduct"=>getDegreeProduct,
             "DegreeSum"=>getDegreeSum,
             "SPGreedy"=>getMyCent,
             "FastGreedy"=>getApproxMyCent)
    ans = []
    G = deepcopy(G0)
    for rep = 1 : k
        E = collect(edges(G))
        C = F[method](G)
        x = argmax(C)
        (u,v) = Tuple(E[x])
        push!(ans, (u,v))
        rem_edge!(G, u, v)
    end
    return ans
end

function GreedyEPS(G0, k; set_eps=0.1)
    ans = []
    G = deepcopy(G0)
    for rep = 1 : k
        E = collect(edges(G))
        C = getApproxMyCent(G, eps=set_eps)
        x = argmax(C)
        (u,v) = Tuple(E[x])
        push!(ans, (u,v))
        rem_edge!(G, u, v)
    end
    return ans
end

function RandomSelect(G, k)
    E = collect(edges(G))
    m = ne(G)
    ans = []
    A = randperm(m)
    for i = 1 : k
        x = A[i]
        (u,v) = Tuple(E[x])
        push!(ans, (u,v))
    end
    return ans
end

function Opt(G0, k)
    G = deepcopy(G0)
    E = [Tuple(x) for x in edges(G)]
    n,m = nv(G),ne(G)
    del = zeros(Bool, m)
    SW = [inv(laplacian_matrix(G)+I(n)+zeros(n,n))]
    ret = []
    best = 0.0

    get_new(id) = begin
        W = last(SW)
        (u,v) = E[id]
        b_e = W[:,u]-W[:,v]
        alp = 1.0 / (1.0 - (b_e[u]-b_e[v]))
        push!(SW, W + alp * b_e * b_e')
    end
    
    dfs(dep, pre) = begin
        if dep>k
            cnt_score = tr(last(SW))
            if cnt_score>best
                best = cnt_score
                ret = []
                for (i,x) in enumerate(E)
                    if del[i]
                        push!(ret, x)
                    end
                end
            end
            return nothing
        end
        for i = pre : m
            del[i] = true
            get_new(i)
            dfs(dep+1, i+1)
            pop!(SW)
            del[i] = false
        end
    end

    dfs(1, 1)

    return ret
end