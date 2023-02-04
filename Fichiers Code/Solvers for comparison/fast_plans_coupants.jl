using JuMP
using CPLEX
include("fast_sp_solving.jl")

function fast_plans_coupants(inputFile::String, TimeLimit::Int64)

    include(inputFile)

    start = time()
    remaining_time = int(TimeLimit - (time() - start))

    #definition de l
    l = [sqrt((coordinates[i, 1]-coordinates[j, 1])^2 + (coordinates[i, 2] - coordinates[j, 2])^2) for i in 1:n, j in 1:n]

    U1 = [l]
    U2 = [convert(Vector{Float64}, w_v)]

    # Première itération

    m = Model(CPLEX.Optimizer)

    #Variables
    @variable(m, z >= 0)
    @variable(m, x[i in 1:n, j in i+1:n], Bin)
    @variable(m, y[k in 1:K, i in k:n], Bin)

    #Contraintes Xcomb
    @constraint(m, [l1 in U1], z >= sum(l1[i,j]*x[i,j] for i in 1:n,j in i+1:n))

    @constraint(m, [k in 1:K, i in k:n, j in i+1:n], y[k,i] + y[k,j] <= 1+x[i,j])

    @constraint(m, [i in 1:K], sum(y[k,i] for k in 1:i) == 1)
    @constraint(m, [i in K+1:n], sum(y[k,i] for k in 1:K) == 1)

    @constraint(m, [k in 1:K, w2_v in U2], sum(w2_v[i]*y[k,i] for i in k:n) <= B)

    #Objectif
    @objective(m, Min, z)

    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)

    # Définit une limite de temps
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", remaining_time)

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

    remaining_time = int(TimeLimit - (time() - start))

    # Résolution de SP1 et récupération des données
    z1, delta1_val = fast_SP1_solve(inputFile, x_val, l, remaining_time)

    # Actualisation de U1
    if z1 > z_val
        l1 = l .+ delta1_val .*(lh .+ transpose(lh))
        push!(U1, l1)
    end

    remaining_time = int(TimeLimit - (time() - start))

    # Résolution de SP2_k et récupération des données
    z2, delta2_val = fast_SP2k_solve(inputFile, y_val, remaining_time)

    # Actualisation de U2
    if maximum(z2) > B
        for k in 1:K
            w2 = w_v .* (delta2_val[k,:] .+ 1)
            push!(U2, w2)
        end
    end

    # Nouvelles itérations de résolution du problème maître

    while (z1 > z_val || maximum(z2) > B) && time() - start < TimeLimit
        
        m = Model(CPLEX.Optimizer)

        #Variables
        @variable(m, z >= 0)
        @variable(m, x[i in 1:n, j in i+1:n], Bin)
        @variable(m, y[k in 1:K, i in k:n], Bin)

        #Contraintes Xcomb
        for l1 in U1
            @constraint(m, z >= sum(l1[i,j]*x[i,j] for i in 1:n,j in i+1:n))
        end

        @constraint(m, [k in 1:K, i in k:n, j in i+1:n], y[k,i] + y[k,j] <= 1+x[i,j])

        @constraint(m, [i in 1:K], sum(y[k,i] for k in 1:i) == 1)
        @constraint(m, [i in K+1:n], sum(y[k,i] for k in 1:K) == 1)

        for w2 in U2
            @constraint(m, [k in 1:K], sum(w2[i]*y[k,i] for i in k:n) <= B)
        end

        #Objectif
        @objective(m, Min, z)

        # Désactive les sorties de CPLEX (optionnel)
        set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)

        # Définit une limite de temps
        remaining_time = int(TimeLimit - (time() - start))
        set_optimizer_attribute(m, "CPX_PARAM_TILIM", remaining_time)

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
        remaining_time = int(TimeLimit - (time() - start))

        z1, delta1_val = fast_SP1_solve(inputFile, x_val, l, remaining_time)

        # Actualisation de U1
        if z1 > z_val
            l1 = l .+ delta1_val .*(lh .+ transpose(lh))
            push!(U1, l1)
        end

        # Résolution de SP2_k et récupération des données
        remaining_time = int(TimeLimit - (time() - start))
        z2, delta2_val = fast_SP2k_solve(inputFile, y_val, remaining_time)

        # Actualisation de U2
        for k in 1:K
            if z2[k] > B
                w2 = w_v .* (delta2_val[k,:] .+ 1)
                push!(U2, w2)
            end
        end
    end

    computation_time = time() - start

    return z_1, computation_time
end