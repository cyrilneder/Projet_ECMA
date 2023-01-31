#Heuristique de recherche locale permettant la résolution du problème statique

#Opérateurs de voisinage considérés (à condition qu'ils soient permis par les contraintes sur les poids): 
# - Swap de deux sommets
# - Changement d'ensemble d'un sommets


mutable struct Solution
    n::Int64                                    #Amount of vertices
    l::Matrix{Float64}                          #Distance matrix
    vx::Matrix{Int64}                           #Current value of x variable
    vy::Matrix{Int64}                           #Current value of y variable
    K::Int64                                    #Number of sets allowed
    Nsets::Int64                                #Number of sets currently used
    cost::Float64                               #Cost of the solution

    function Solution() = new()
end

function update_Nsets!(sol::Solution)
    K = sol.K
    Nsets = 0
    for k in 1::K
        used_set = (sum(y[k,:]) >= 1)
        if used_set
            Nsets += 1
        end
    end
    sol.Nsets = Nsets
end

function update_cost!(sol::Solution)
    #Returns the cost of the current solution after having modified the vx and vy attributes
    n = sol.n
    sol.cost = sum(sol.l[i,j]*sol.vx[i,j] for i in 1:n, j in i+1:n)
end



function swap!(sol::Solution, s1::Int64, s2::Int64)
    K = sol.K
    k1 = findfirst(lamb -> lamb==1, sol.y[:,s1])
    k2 = findfirst(lamb -> lamb==1, sol.y[:,s2])

    


end