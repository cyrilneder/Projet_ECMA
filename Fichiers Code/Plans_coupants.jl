using JuMP
using CPLEX
include("sp_solving.jl")

path = "data 2/10_ulysses_3.tsp"

function plans_coupants(inputFile::String)

    include(inputFile)

    #definition de l
    l = [sqrt((coordinates[i, 1]-coordinates[j, 1])^2 + (coordinates[i, 2] - coordinates[j, 2])^2) for i in 1:n, j in 1:n]

    U1 = [l]
    U2 = [convert(Vector{Float64}, w_v)]

    # Première itération

    m = Model(CPLEX.Optimizer)

    #Variables
    @variable(m, z >= 0)
    @variable(m, x[i in 1:n, j in i+1:n], Bin)
    @variable(m, y[k in 1:K, i in 1:n], Bin)

    #Contraintes Xcomb
    @constraint(m, [l1 in U1], z >= sum(l1[i,j]*x[i,j] for i in 1:n,j in i+1:n))

    @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] + y[k,j] <= 1+x[i,j])
    @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] - y[k,j] <= 1-x[i,j])
    @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], -y[k,i] + y[k,j] <= 1-x[i,j])

    @constraint(m, [i in 1:n], sum(y[k,i] for k in 1:K) == 1)

    @constraint(m, [k in 1:K, w2_v in U2], sum(w2_v[i]*y[k,i] for i in 1:n) <= B)

    #Objectif
    @objective(m, Min, z)

    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)

    optimize!(m)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
            # Récupération des valeurs d’une variable
            x_val = JuMP.value.(x)
            y_val = JuMP.value.(y)
            z_val = JuMP.objective_value(m)
    end

    # Résolution de SP1 et récupération des données
    z1, delta1_val = SP1_solve(inputFile, x_val, l)

    # Actualisation de U1
    if z1 > z_val
        l1 = l .+ delta1_val .*(lh .+ transpose(lh))
        println("this is diff : ", l1-l)
        push!(U1, l1)
        println("Add constraint type SP1 because of objective breach ", z1)
    end

    # Résolution de SP2_k et récupération des données
    z2, delta2_val = SP2k_solve(inputFile, y_val)

    # Actualisation de U2
    if maximum(z2) > B
        for k in 1:K
            w2 = w_v .* (delta2_val[k,:] .+ 1)
            push!(U2, w2)
        end
        println("Add constraint type SP2k because of cluster weight breach : ", maximum(z2))
    end

    # Nouvelles itérations de résolution du problème maître

    while z1 > z_val || maximum(z2) > B
        
        m = Model(CPLEX.Optimizer)

        #Variables
        @variable(m, z >= 0)
        @variable(m, x[i in 1:n, j in i+1:n], Bin)
        @variable(m, y[k in 1:K, i in 1:n], Bin)

        #Contraintes Xcomb
        for l1 in U1
            @constraint(m, z >= sum(l1[i,j]*x[i,j] for i in 1:n,j in i+1:n))
        end

        @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] + y[k,j] <= 1+x[i,j])
        @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] - y[k,j] <= 1-x[i,j])
        @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], -y[k,i] + y[k,j] <= 1-x[i,j])

        @constraint(m, [i in 1:n], sum(y[k,i] for k in 1:K) == 1)

        for w2 in U2
            @constraint(m, [k in 1:K], sum(w2[i]*y[k,i] for i in 1:n) <= B)
        end

        #Objectif
        @objective(m, Min, z)

        # Désactive les sorties de CPLEX (optionnel)
        set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)

        optimize!(m)

        # Récupération du status de la résolution
        feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
        isOptimal = termination_status(m) == MOI.OPTIMAL
        if feasibleSolutionFound
                # Récupération des valeurs d’une variable
                x_val = JuMP.value.(x)
                y_val = JuMP.value.(y)
                z_val = JuMP.objective_value(m)
        end

        # Résolution de SP1 et récupération des données
        z1, delta1_val = SP1_solve(inputFile, x_val, l)

        # Actualisation de U1
        if z1 > z_val
            l1 = l .+ delta1_val .*(lh .+ transpose(lh))
            push!(U1, l1)
            println("Add constraint type SP1 because of objective breach ", z1)
            #println("diff is :", l1)
        end

        # Résolution de SP2_k et récupération des données
        z2, delta2_val = SP2k_solve(inputFile, y_val)

        # Actualisation de U2
        if maximum(z2) > B
            for k in 1:K
                w2 = w_v .* (delta2_val[k,:] .+ 1)
                push!(U2, w2)
            end
            println("Add constraint type SP2k because of cluster weight breach : ", maximum(z2))
        end        
    end

    #Affichage de solution
    println("Valeur optimale :", z_val," obtenue avec les clusters : ")
    for k in 1:K
        println(findall(!iszero, y_val[k,:]))
    end
end

plans_coupants(path)