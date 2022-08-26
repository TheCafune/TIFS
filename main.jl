include("Exp.jl")

function run_small()
    for gname in ["Tribes", "Karate", "Southernwomen", "Dolphins"]
        lg = open("../results/result_small.txt", "a")
        test_small(gname, 6, lg=lg)
        close(lg)
    end
end

function run_mid()
    for gname in ["EmailUniv", "GridWorm", "GrQc", "WikiElec"]
        lg = open("../results/result_mid.txt", "a")
        test_mid(gname, 50, lg=lg, query=5:5:50)
        close(lg)
    end
end

function run_large()
    gList = ["EmailUniv", "Erdos992", "Bcspwr10", "Reality", "PagesGovernment", "Dmela", "HepPh",
    "Anybeat", "PagesCompany", "AstroPh", "CondMat", "Gplus", "Douban", "Gowalla", "GooglePlus", "Citeseer",
    "MathSciNet", "TwitterFollows", "YoutubeSnap", "Lastfm", "Flixster"]
    for gname in gList
        for eps in [0.1, 0.2, 0.3]
            lg = open("../results/result_large.txt", "a")
            test_large(gname, 50, eps, lg=lg)
            close(lg)
        end
    end
end

run_small()
run_mid()
run_large()