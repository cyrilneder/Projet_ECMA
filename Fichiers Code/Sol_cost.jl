#Heuristique de recherche locale permettant la résolution du problème statique

#Opérateurs de voisinage considérés (à condition qu'ils soient permis par les contraintes sur les poids): 
# - Swap de deux sommets
# - Changement d'ensemble d'un sommets

import Base.copy


mutable struct Solution
    n::Int64                                    #Amount of vertices
    l::Matrix{Float64}                          #Distance matrix
    B::Int64                                    #Maximum allowed weight
    w_v::Vector{Int64}                          #Weight per vertex
    vx::Vector{Vector{Int64}}                   #Current value of x variable
    vy::Vector{Vector{Int64}}                   #Current value of y variable
    K::Int64                                    #Number of sets allowed
    Nsets::Int64                                #Number of sets currently used
    cost::Float64                               #Cost of the solution

    Solution() = new()
end

function parse(inst::String)
    include("data 2/"*inst)
    sol = Solution()
    sol.n = n
    sol.l = [sqrt((coordinates[i, 1]-coordinates[j, 1])^2 + (coordinates[i, 2] - coordinates[j, 2])^2) for i in 1:n, j in 1:n]
    sol.B = B
    sol.w_v = w_v
    sol.vx = Vector{Vector{Int64}}()
    sol.vy = Vector{Vector{Int64}}()
    sol.K = K
    sol.Nsets = 0
    sol.cost = maximum(sol.l)*2*n*(n-1)

    return sol

end

function copy(sol::Solution)
    sol2 = Solution()
    sol2.n = copy(sol.n)
    sol2.l = copy(sol.l)
    sol2.B = copy(sol.B)
    sol2.w_v = copy(sol.w_v)
    sol2.vx = copy(sol.vx)
    sol2.vy = copy(sol.vy)
    sol2.K = copy(sol.K)
    sol2.Nsets = copy(sol.Nsets)
    sol2.cost = copy(sol.cost)
    return sol2
end

function verify_weight(sol::Solution,yk::Vector{Int64})
    #yk is a potential new set of vertices ; it needs to satisfy the weight constraint    
    return (sum(yk[i]*w_v[i] for i in 1:sol.n) <= sol.B)
end

function verify_weight(sol::Solution,y::Vector{Vector{Int64}})
    for yk in y
        if !verify_weight(sol,yk)
            return false
        end
    end
    return true
end

function corresponding_x(newy::Vector{Vector{Int64}})
    n = length(newy[1])
    K = length(newy)
    x = [[0 for j in i+1:n] for i in 1:n-1]
    for k in 1:K
        vk = []
        for i in 1:n
            if newy[k][i] == 1
                push!(vk, i)
            end
        end
        len = length(vk)
        for i in 1:len-1
            for j in i+1:len
                x[vk[i]][vk[j] - vk[i]] = 1
            end
        end
    end
    return x
end

function update_Nsets!(sol::Solution)
    K = sol.K
    Nsets = 0
    for k in 1:K
        used_set = (sum(sol.vy[k]) >= 1)
        if used_set
            Nsets += 1
        end
    end
    sol.Nsets = Nsets
end

function update_cost!(sol::Solution)
    #Returns the cost of the current solution after having modified the vx and vy attributes
    n = sol.n
    res = 0
    for i in 1:n-1
        res += sum(sol.l[i,j]*sol.vx[i][j - i] for j in i+1:n)
    end
    sol.cost = res
end



function swapped(sol::Solution, s1::Int64, s2::Int64)
    println(sol.vy)
    println(s1)
    println(s2)
    k1 = findfirst(lamb -> lamb==1, [sol.vy[k][s1] for k in 1:sol.K])
    k2 = findfirst(lamb -> lamb==1, [sol.vy[k][s2] for k in 1:sol.K])
    yp = copy(sol.vy)

    yp[k1][s1] = 0
    yp[k1][s2] = 1
    yp[k2][s2] = 0
    yp[k2][s1] = 1

    sp = copy(sol)
    if verify_weight(sp,yp)
        sp.vx = corresponding_x(yp)
        sp.vy = yp
        update_cost!(sp)
    end
    return sp    
end

function changed_set(sol::Solution, s1::Int64, k::Int64)
    k1 = findfirst(lamb -> lamb==1, [sol.vy[i][s1] for i in 1:sol.K])
    yp = copy(sol.vy)

    yp[k1][s1] = 0
    yp[k][s1] = 1
    println("yp :",yp)
    println("vy :",sol.vy)

    sp = copy(sol)
    if verify_weight(sp,yp)
        sp.vx = corresponding_x(yp)
        sp.vy = yp
        update_cost!(sp)
    end
    return sp
end


function real_sol(n::Int64, l::Matrix{Float64}, B::Int64, w_v::Vector{Int64}, K::Int64)
    sol = Solution()
    sol.n = n
    sol.l = l
    sol.B = B
    sol.w_v = w_v
    sol.K = K
    sol.cost = maximum(sol.l)*2*n*(n-1)

    #Heuristique ne garantissant pas une solution réalisable (plutôt à faire avec la PPC)

    to_place = Vector{Int64}(1:n)
    to_place_w = copy(w_v)
    k = 1
    y = [[0 for _ in 1:n] for _ in 1:K]

    while !isempty(to_place) && k <= K
        yk = copy(y[k])
        current_weight = 0
        addable = findall(lamb -> lamb <= B, to_place_w)
        if isempty(addable)
            fillable = false
        else
            s = findfirst(lamb -> lamb == maximum(to_place_w[findall(lamb -> lamb <= B, to_place_w)]), to_place_w)
            fillable = true
        end
        while fillable
            yk[to_place[s]] = 1
            current_weight += w_v[s]
            deleteat!(to_place_w, s)
            deleteat!(to_place, s)

            addable = findall(lamb -> lamb + current_weight <= B, to_place_w)
            if isempty(addable)
                fillable = false
            else
                s = findfirst(lamb -> lamb == maximum(to_place_w[findall(lamb -> lamb + current_weight <= B, to_place_w)]), to_place_w)
            end
        end
        y[k] = copy(yk)
        k += 1

    end
    if isempty(to_place)
        println("Placement fini")
        sol.vy = y
        sol.vx = corresponding_x(y)
        update_Nsets!(sol)
        update_cost!(sol)
        return sol
    else
        println("Solution réalisable non trouvée.")
        return -1
    end

end

function real_sol!(sol::Solution)
    #Ici sol fait office d'instance de départ, sans valeur attribuée à vx et vy
    
    #Heuristique ne garantissant pas une solution réalisable (plutôt à faire avec la PPC)
    n = sol.n
    B = sol.B
    w_v = sol.w_v
    K = sol.K

    to_place = Vector{Int64}(1:n)
    to_place_w = copy(w_v)
    k = 1
    y = [[0 for _ in 1:n] for _ in 1:K]

    while !isempty(to_place) && k <= K
        yk = copy(y[k])
        current_weight = 0
        addable = findall(lamb -> lamb <= B, to_place_w)
        if isempty(addable)
            fillable = false
        else
            s = findfirst(lamb -> lamb == maximum(to_place_w[findall(lamb -> lamb <= B, to_place_w)]), to_place_w)
            fillable = true
        end
        while fillable
            yk[to_place[s]] = 1
            current_weight += w_v[to_place[s]]
            #println("current weight:",current_weight)
            #println("added vertex:",to_place[s])
            deleteat!(to_place_w, s)
            deleteat!(to_place, s)

            addable = findall(lamb -> lamb + current_weight <= B, to_place_w)
            if isempty(addable)
                fillable = false
            else
                s = findfirst(lamb -> lamb == maximum(to_place_w[findall(lamb -> lamb + current_weight <= B, to_place_w)]), to_place_w)
            end
        end
        y[k] = copy(yk)
        k += 1

    end
    if isempty(to_place)
        println("Placement fini")
        sol.vy = y
        sol.vx = corresponding_x(y)
        update_Nsets!(sol)
        update_cost!(sol)
        return sol
    else
        println("Solution réalisable non trouvée.")
        return -1
    end


end

function print_sets(sol::Solution)
    for k in 1:sol.K
        vk = Vector{Int64}([])
        for i in 1:sol.n
            if sol.vy[k][i] == 1
                push!(vk,i)
            end
        end
        println(vk)
        if isempty(vk)
            println(0)
        else
            println(sum(sol.w_v[vk]))
        end
    end
end

function vns(sol::Solution ; dur::Int64 = 30)
    #A ce stade sol est une solution réalisable de valeur très peu optimisée
    #On réalise une descente à voisinages variables déterministe pour naviguer entre optima locaux

    startime = time_ns()/1000000000
    nit = 0

    cursol = copy(sol)
    bestsol = copy(sol)
    testsol = copy(sol)
    n = sol.n
    K = sol.K

    currentime = time_ns()/1000000000
    finished = (currentime - startime >= dur)

    while !finished
        nit += 1
        #Essai de shift 
        for s in 1:n
            for k in 1:K
                testsol = changed_set(cursol, s, k)
                update_cost!(testsol)
                update_Nsets!(testsol)
                if testsol.cost < bestsol.cost
                    bestsol = copy(testsol)
                end
                testsol = copy(cursol)
            end
        end

        #Essai de swap
        for s1 in 1:n
            for s2 in 1:n
                testsol = swapped(cursol, s1, s2)
                update_cost!(testsol)
                update_Nsets!(testsol)
                if testsol.cost < bestsol.cost
                    bestsol = copy(testsol)
                end
                testsol = copy(cursol)
            end
        end

        cursol = copy(bestsol)

        currentime = time_ns()/1000000000
        finished = (currentime - startime >= dur)
    end
    println("Nombre d'itérations :",nit)

    return cursol
end
