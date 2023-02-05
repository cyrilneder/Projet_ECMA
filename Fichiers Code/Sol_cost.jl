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
    W_v::Vector{Float64}                        #Indecision factor for weight
    lh::Vector{Int64}                           #Indecision factor for length
    W::Int64                                    #Maximum indecision factor for weight
    L::Int64                                    #Maximum indecision factor for length
    vx::Vector{Vector{Int64}}                   #Current value of x variable
    vy::Vector{Vector{Int64}}                   #Current value of y variable
    delta1::Vector{Vector{Float64}}             #Current value of delta_one variable
    K::Int64                                    #Number of sets allowed
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
    sol.W_v = W_v
    sol.lh = lh
    sol.W = W
    sol.L = L
    sol.vx = Vector{Vector{Int64}}()
    sol.vy = Vector{Vector{Int64}}()
    sol.delta1 = Vector{Vector{Float64}}()
    sol.K = K
    sol.cost = maximum(sol.l)*2*n*(n-1)

    return sol

end

function copy(sol::Solution)
    sol2 = Solution()
    sol2.n = copy(sol.n)
    sol2.l = copy(sol.l)
    sol2.B = copy(sol.B)
    sol2.w_v = copy(sol.w_v)
    sol2.W_v = copy(sol.W_v)
    sol2.lh = copy(sol.lh)
    sol2.W = copy(sol.W)
    sol2.L = copy(sol.L)
    sol2.vx = deepcopy(sol.vx)
    sol2.vy = deepcopy(sol.vy)
    sol2.delta1 = deepcopy(sol.delta1)
    sol2.K = copy(sol.K)
    sol2.cost = copy(sol.cost)
    return sol2
end

function verify_weight(sol::Solution,yk::Vector{Int64})
    if sum(yk) == 0
        return true
    end
    #yk is a potential new set of vertices ; it needs to satisfy the weight constraint
    #On force la robustesse ici : yk doit vérifier la contrainte de poids dans le cas de la pire instanciation possibles
    
    #Pour construire la pire instanciation de poids, on agit récursivement : on trie les poids par ordre décroissant,
    #et on augmente au maximum le coefficient delta2 du sommet de poids maximum (tant qu'il est inférieur à W_v et que
    #la somme des coefficients est inférieure à W) et ainsi de suite jusqu'à ne plus pouvoir augmenter

    delta2 = [0.0 for _ in 1:sol.n]
    sorted_w = sort(sol.w_v, rev = true)
    addable_w = sol.W

    vk = Vector{Int64}([])
    for i in 1:sol.n
        if yk[i] == 1
            push!(vk,i)
        end
    end
    sorted_w = sort(sol.w_v[vk], rev = true)

    worst_case = false
    while !worst_case
        max_w = popfirst!(sorted_w)
        ind = findfirst(lamb -> sol.w_v[lamb] == max_w, vk)
        delta2[vk[ind]] = min(sol.W_v[vk[ind]], addable_w)
        addable_w -= delta2[vk[ind]]

        worst_case = isempty(sorted_w) || (addable_w == 0)
    end

    worst_weight = sum(yk[i]*w_v[i]*(1 + delta2[i]) for i in 1:sol.n)
    return (worst_weight <= sol.B)
end

function verify_weight(sol::Solution,y::Vector{Vector{Int64}})
    for yk in y
        if !verify_weight(sol,yk)
            return false
        end
    end
    return true
end

function verify_weights(sol::Solution)
    return verify_weight(sol, sol.vy)
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

function update_cost!(sol::Solution)
    #Returns the cost of the current solution after having modified the vx and vy attributes
    n = sol.n
    res = 0
    for i in 1:n-1
        res += sum((sol.l[i,j] + sol.delta1[i][j - i]*(sol.lh[i] + sol.lh[j]))*sol.vx[i][j - i] for j in i+1:n)
    end
    sol.cost = res
end



function swapped(sol::Solution, s1::Int64, s2::Int64)
    sp = copy(sol)

    k1 = findfirst(lamb -> lamb==1, [sol.vy[k][s1] for k in 1:sol.K])
    k2 = findfirst(lamb -> lamb==1, [sol.vy[k][s2] for k in 1:sol.K])

    if k1 != k2
        yp = deepcopy(sol.vy)

        yp[k1][s1] = 0
        yp[k2][s2] = 0
        yp[k1][s2] = 1
        yp[k2][s1] = 1

        if verify_weight(sp,yp)
            sp.vx = corresponding_x(yp)
            sp.vy = yp
            update_cost!(sp)
        end
    end
    return sp    
end

function changed_set(sol::Solution, s1::Int64, k::Int64)

    k1 = findfirst(lamb -> lamb==1, [sol.vy[i][s1] for i in 1:sol.K])
    sp = copy(sol)

    if k1 != k
        yp = deepcopy(sol.vy)

        yp[k1][s1] = 0
        yp[k][s1] = 1

        if verify_weight(sp,yp)
            sp.vx = corresponding_x(yp)
            sp.vy = yp
            update_cost!(sp)
        end
    end
    return sp
end


function real_sol!(sol::Solution)
    #Ici sol fait office d'instance de départ, sans valeur attribuée à vx et vy
    #Le problème considéré de prime abord est le problème statique
    
    #Heuristique ne garantissant pas une solution réalisablea priori (plutôt à faire avec la PPC), mais marche bien sur les instances
    n = sol.n
    B = sol.B
    w_v = sol.w_v
    W_v = sol.W_v
    W = sol.W
    K = sol.K

    to_place = Vector{Int64}(1:n)
    to_place_w = copy(w_v)
    k = 1
    y = [[0 for _ in 1:n] for _ in 1:K]

    while !isempty(to_place) && k <= K
        yk = copy(y[k])
        current_worse_w = 0
        robust_w = 0
        addable = findall(lamb -> to_place_w[lamb]*(1 + min(W_v[to_place[lamb]], W - robust_w)) <= B, 1:length(to_place))
        if isempty(addable)
            fillable = false
        else
            robust_weights = [to_place_w[i]*(1 + min(W_v[to_place[i]], W - robust_w)) for i in 1:length(to_place)]
            s = findfirst(lamb -> lamb == maximum(robust_weights[addable]), robust_weights)
            fillable = true
        end
        while fillable
            yk[to_place[s]] = 1
            delta2s = min(W_v[to_place[s]], W - robust_w)
            robust_w += delta2s
            current_worse_w += w_v[to_place[s]]*(1 + delta2s)
            deleteat!(to_place_w, s)
            deleteat!(to_place, s)

            addable = findall(lamb -> to_place_w[lamb]*(1 + min(W_v[to_place[lamb]], W - robust_w)) + current_worse_w <= B, 1:length(to_place))
            if isempty(addable)
                fillable = false
            else
                robust_weights = [to_place_w[i]*(1 + min(W_v[to_place[i]], W - robust_w)) for i in 1:length(to_place)]
                s = findfirst(lamb -> lamb == maximum(robust_weights[addable]), robust_weights)
            end
        end
        y[k] = copy(yk)
        k += 1

    end
    if isempty(to_place)
        println("Placement fini")
        sol.vy = y
        sol.vx = corresponding_x(y)
        sol.delta1 = [[0 for j in i+1:n] for i in 1:n-1]
        update_cost!(sol)
        return sol
    else
        println("Solution réalisable non trouvée.")
        return -1
    end


end

function print_sets(sol::Solution)
    if !verify_weights(sol)
        println("\n \n \n \n ALERTE : ENSEMBLES NON ROBUSTES \n \n \n \n")
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
                println("Poids statique de l'ensemble:",sum(sol.w_v[vk]))

                delta2 = [0.0 for _ in 1:sol.n]
                sorted_w = sort(sol.w_v, rev = true)
                addable_w = sol.W
                sorted_w = sort(sol.w_v[vk], rev = true)
            
                worst_case = false
                while !worst_case
                    max_w = popfirst!(sorted_w)
                    ind = findfirst(lamb -> sol.w_v[lamb] == max_w, vk)
                    delta2[vk[ind]] = min(sol.W_v[vk[ind]], addable_w)
                    addable_w -= delta2[vk[ind]]
            
                    worst_case = isempty(sorted_w) || (addable_w == 0)
                end
            
                worst_weight = sum(sol.vy[k][i]*w_v[i]*(1 + delta2[i]) for i in 1:sol.n)
                println("Poids dans le pire cas :",worst_weight)
                if worst_weight > sol.B
                    print("\n ERREUR \n POIDS DE L'ENSEMBLE A VERIFIER \n \n")
                end
            end

        end
    end
end

function vns_descent(sol::Solution ; dur::Int64 = 30)
    #A ce stade sol est une solution réalisable de valeur très peu optimisée
    #On réalise une descente à voisinages variables déterministe pour naviguer entre optima locaux
    #Cet algorithme présente la particularité d'aboutir assez rapidement à un plateau

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
                if testsol.cost < bestsol.cost
                    bestsol = copy(testsol)
                end
                testsol = copy(cursol)
            end
        end

        #Essai de swap
        for s1 in 1:n-1
            for s2 in s1+1:n
                testsol = swapped(cursol, s1, s2)
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


function vns(sol::Solution ; dur::Int64 = 10)
    #A ce stade sol est une solution réalisable de valeur très peu optimisée
    #On réalise une recherche à voisinage variable avec un aspect aléatoire

    startime = time_ns()/1000000000
    nit = 0

    cursol = copy(sol)
    randsol = copy(sol)
    bestsol = copy(sol)
    testsol = copy(sol)
    n = sol.n
    K = sol.K

    currentime = time_ns()/1000000000
    finished = (currentime - startime >= dur)

    while !finished
        nit += 1

        neigh = rand(1:2)

        if neigh == 1
            #Voisin obtenu à l'aide d'un shift
            s = rand(1:n)
            k = rand(1:K)
            randsol = changed_set(cursol, s, k)
        else
            #Voisin obtenu à l'aide d'un swap
            s1 = rand(1:n-1)
            s2 = rand(s1+1:n)
            randsol = swapped(cursol, s1, s2)
        end
        
        #Essai de shift 
        for s in 1:n
            for k in 1:K
                testsol = changed_set(randsol, s, k)
                update_cost!(testsol)
                if testsol.cost < bestsol.cost
                    bestsol = copy(testsol)
                end
                testsol = copy(randsol)
            end
        end

        #Essai de swap
        for s1 in 1:n-1
            for s2 in s1+1:n
                testsol = swapped(randsol, s1, s2)
                if testsol.cost < bestsol.cost
                    bestsol = copy(testsol)
                end
                testsol = copy(randsol)
            end
        end

        cursol = copy(bestsol)

        currentime = time_ns()/1000000000
        finished = (currentime - startime >= dur)
    end
    println("Nombre d'itérations :",nit)

    return cursol
end


function robust_model!(sol::Solution)
    #Le but de cette fonction est d'instancier les valeurs de delta1 à leur "pire" valeur possible, c'est-à-dire celle
    #maximisant la valeur de l'objectif

    n = sol.n
    robust_length = sol.L

    sol.delta1 = [[0 for j in i+1:n] for i in 1:n-1]

    robust_fac = Dict((i,j) => sol.lh[i] + sol.lh[j] for i in 1:n-1 for j in i+1:n)
    sorted_rf = sort(collect(robust_fac), by = lamb -> lamb[2], rev = true)
    while length(sorted_rf) > 0 && robust_length > 0
        arc = popfirst!(sorted_rf)
        (i,j) = arc[1]
        if sol.vx[i][j - i] == 1
            addable_length = min(3, robust_length)
            sol.delta1[i][j - i] = addable_length
            robust_length -= addable_length
        end
    end

    update_cost!(sol)
    println("Coût du pire scenario :",sol.cost)
end


function robust_opt(inst::String ; dur::Int64 = 60)
    sol = parse(inst)
    real_sol!(sol)
    #A ce stade sol est une solution respectant la robustesse des contraintes de poids

    bestsol = vns(sol ; dur = 5)

    robust_model!(bestsol)
    #A ce stade le coût de bestsol est bien en accord avec la robustesse des longueurs
    println("Coût de la solution robuse :", bestsol.cost)

    cursol = copy(bestsol)
    startime = time_ns()/1000000000

    finished = false

    while !finished

        cursol = vns(bestsol ; dur = 5)
        #A ce stade cursol est de coût inférieur à celui de bestsol

        robust_model!(cursol)
        #La robustesse des longueurs peut avoir fait de cursol une moins bonne solution que bestsol

        if cursol.cost < bestsol.cost
            bestsol = copy(cursol)
        end

        currentime = time_ns()/1000000000
        finished = (currentime - startime >= dur)
        println("Coût de la meilleure solution robuste :",bestsol.cost)
    end

    return sol
end


