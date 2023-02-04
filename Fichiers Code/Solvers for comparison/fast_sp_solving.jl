using JuMP
using CPLEX

function fast_SP1_solve(inputFile::String, x_val, l, TimeLimit::Int64)

    include(inputFile)

    #Sous-problème 1
    SP1 = Model(CPLEX.Optimizer)

    @variable(SP1, delta1[i in 1:n, j in i+1:n] >= 0)

    @constraint(SP1, sum(delta1[i,j] for i in 1:n,j in i+1:n) <= L)
    @constraint(SP1, [i in 1:n, j in i+1:n], delta1[i,j] <= 3)

    @objective(SP1, Max, sum((l[i,j] + delta1[i,j]*(lh[i]+lh[j]))*x_val[i,j]  for i in 1:n,j in i+1:n))

    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(SP1, "CPX_PARAM_SCRIND", 0)

    # Définit une limite de temps
    set_optimizer_attribute(SP1, "CPX_PARAM_TILIM", TimeLimit)

    #Résolution
    optimize!(SP1)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(SP1) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(SP1) == MOI.OPTIMAL
    if feasibleSolutionFound
            # Récupération des valeurs d’une variable
            delta1_val = JuMP.value.(delta1)
            z1 = JuMP.objective_value(SP1)
    end

    # Conversion de delta1_val car type pas approprié pour plans_coupants
    delta1_val_converted = zeros(n,n)

    for i in 1:n
        for j in i+1:n
            delta1_val_converted[i,j] = delta1_val[i,j]
        end
    end
    return z1, delta1_val_converted
end

function fast_SP2k_solve(inputFile::String, y_val, TimeLimit::Int64)

    include(inputFile)

    #Sous-problème 2,k
    delta2_val = Matrix{Float64}(undef,K,n)
    z2 =[]

    for k in 1:K
        SP2k = Model(CPLEX.Optimizer)

        @variable(SP2k, delta2k[i in k:n] >= 0)

        @constraint(SP2k, sum(delta2k[i] for i in k:n) <= W)
        @constraint(SP2k, [i in k:n], delta2k[i] <= W_v[i])

        @objective(SP2k, Max, sum(w_v[i]*(1 + delta2k[i])*y_val[k,i]  for i in k:n))

        # Désactive les sorties de CPLEX (optionnel)
        set_optimizer_attribute(SP2k, "CPX_PARAM_SCRIND", 0)

        # Définit une limite de temps
        set_optimizer_attribute(SP2k, "CPX_PARAM_TILIM", TimeLimit)

        #Résolution
        optimize!(SP2k)

        # Récupération du status de la résolution
        feasibleSolutionFound = primal_status(SP2k) == MOI.FEASIBLE_POINT
        isOptimal = termination_status(SP2k) == MOI.OPTIMAL
        if feasibleSolutionFound
                # Récupération des valeurs d’une variable
                delta2_val[k,k:end] = JuMP.value.(delta2k)
                push!(z2, JuMP.objective_value(SP2k))
        end
    end

    return z2, delta2_val
end
