include("Util.jl")
include("Algorithm.jl")

using Printf

function sc(x)
    return @sprintf("%.4f", x)
end

function test_small(gname, maxk; lg=stdout)
    G = readGraph(gname)
    println(lg, gname, " ", nv(G), " ", ne(G), " ", maxk)
    println(lg, "OPT : ")
    for k = 1 : maxk
        ans = Opt(G, k)
        result = sc(getDelta(G, ans))
        println(lg, "$k $result")
    end
    ans_rd = RandomSelect(G, maxk)
    ans_sg = Greedy(G, maxk, method="SPGreedy")
    ans_fg = Greedy(G, maxk, method="FastGreedy")
    A = ["Random", "SPGreedy", "FastGreedy"]
    B = [ans_rd, ans_sg, ans_fg]
    for (s,x) in zip(A,B)
        println(lg, s, " : ")
        for k = 1 : maxk
            result = sc(getDelta(G, x[1:k]))
            println(lg, "$k $result")
        end
    end
end

function test_mid(gname, maxk; lg=stdout, query=1:1:maxk)
    G = readGraph(gname)
    println(lg, gname, " ", nv(G), " ", ne(G), " ", maxk)
    ans_rd = RandomSelect(G, maxk)
    ans_ds = Greedy(G, maxk, method="DegreeSum")
    ans_dp = Greedy(G, maxk, method="DegreeProduct")
    ans_eb = Greedy(G, maxk, method="Betweenness")
    ans_sg = Greedy(G, maxk, method="SPGreedy")
    ans_fg = Greedy(G, maxk, method="FastGreedy")
    A = ["Random", "DegreeSum", "DegreeProduct", "Betweenness", "SPGreedy", "FastGreedy"]
    B = [ans_rd, ans_ds, ans_dp, ans_eb, ans_sg, ans_fg]
    for (s,x) in zip(A,B)
        println(lg, s, " : ")
        for k in query
            result = sc(getDelta(G, x[1:k]))
            println(lg, "$k $result")
        end
    end
end

function test_large(gname, k, eps; lg=stdout)
    G = readGraph(gname)
    println(lg, gname, " ", nv(G), " ", ne(G), " ", k, " ", @sprintf("%.1f", eps))
    if nv(G)<50000
        sg_status = @timed ans_sg = Greedy(G, k, method="SPGreedy")
        fg_status = @timed ans_fg = GreedyEPS(G, k, set_eps=eps)
        sg_score,fg_score = getDelta(G, ans_sg),getDelta(G, ans_fg)
        sg_time,fg_time = sg_status.time,fg_status.time
        ratio = fg_score/sg_score
        println(lg, "SPGreedy Time : ", @sprintf("%.2f", sg_time))
        println(lg, "FastGreedy Time : ", @sprintf("%.2f", fg_time))
        println(lg, "SPGreedy Score : ", @sprintf("%.2f", sg_score))
        println(lg, "FastGreedy Score : ", @sprintf("%.2f", fg_score))
        println(lg, "Ratio : ", sc(ratio))
    else
        fg_status = @timed ans_fg = GreedyEPS(G, k, set_eps=eps)
        fg_time = fg_status.time
        println(lg, "FastGreedy Time : ", @sprintf("%.2f", fg_time))
    end
end