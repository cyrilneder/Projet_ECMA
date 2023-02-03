include("Solvers for comparison/fast_static_solving.jl")
include("Solvers for comparison/fast_dualisation.jl")
include("Solvers for comparison/fast_Branch&Cut.jl")
include("Solvers for comparison/fast_plans_coupants.jl")

TimeLimit = 5 # secondes

println(readdir("Fichiers Code/data 2"))

instance = "22_ulysses_3.tsp"
#instance = "100_kroA_3.tsp"

#println(fast_static_solving("data 2/" * instance, TimeLimit))
#println(fast_dualisation("data 2/" * instance, TimeLimit))
#println(fast_branch_and_cut("data 2/" * instance, TimeLimit))
println(fast_plans_coupants("data 2/" * instance, TimeLimit))




for instance in readdir("Fichiers Code/data 2", sort=false)
    # Afficher le nom de l'instance
    #println(instance)
    return

end