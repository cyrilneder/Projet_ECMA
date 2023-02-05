include("Solvers for comparison/fast_static_solving.jl")
include("Solvers for comparison/fast_dualisation.jl")
include("Solvers for comparison/fast_Branch&Cut.jl")
include("Solvers for comparison/fast_plans_coupants.jl")

# Les 10 premières instances
instances_list_1 = ["10_ulysses_3.tsp", "10_ulysses_6.tsp", "10_ulysses_9.tsp",
"14_burma_3.tsp", "14_burma_6.tsp", "14_burma_9.tsp",
"22_ulysses_3.tsp", "22_ulysses_6.tsp", "22_ulysses_9.tsp",
"26_eil_3.tsp"
]

# Instances 26 <= n <= 30
instances_list_2 = ["26_eil_3.tsp", "26_eil_6.tsp", "26_eil_9.tsp",
"30_eil_3.tsp", "30_eil_6.tsp", "30_eil_9.tsp",
]

# Instances 26 <= n <= 30
instances_list_3 = ["34_pr_3.tsp", "34_pr_6.tsp", "34_pr_9.tsp",
"38_rat_3.tsp", "38_rat_6.tsp", "38_rat_9.tsp"
]

TimeLimit = 600 # secondes
println("Limite de temps considérée : ", TimeLimit, " secondes")

# Warm up
println(fast_dualisation("data 2/10_ulysses_3.tsp", TimeLimit))

for instance in instances_list_1

    println("Traitement de l'instance "*instance)
    Opt, Computing_time = fast_dualisation("data 2/" * instance, TimeLimit)
    println("Meilleure valeur de l'objectif : ", Opt, " atteint en ", Computing_time, " secondes")

end