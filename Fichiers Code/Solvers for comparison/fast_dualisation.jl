using JuMP
using CPLEX

function fast_dualisation(inputFile::String, TimeLimit::Int64)

    include(inputFile)

    start = time()

    #definition de l
    l = [sqrt((coordinates[i, 1]-coordinates[j, 1])^2 + (coordinates[i, 2] - coordinates[j, 2])^2) for i in 1:n, j in 1:n]

    m = Model(CPLEX.Optimizer)

    #Variables primales
    @variable(m, x[i in 1:n, j in i+1:n], Bin)
    @variable(m, y[k in 1:K, i in 1:n], Bin)

    #Variables duales
    @variable(m, beta[i in 1:n, j in i+1:n] >= 0)
    @variable(m, alpha >= 0)
    @variable(m, gamma[k in 1:K] >= 0)
    @variable(m, zeta[k in 1:K, i in 1:n] >= 0)

    #Contraintes Xcomb
    @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] + y[k,j] <= 1+x[i,j])
    #@constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] - y[k,j] <= 1-x[i,j])
    #@constraint(m, [i in 1:n, j in i+1:n, k in 1:K], -y[k,i] + y[k,j] <= 1-x[i,j])

    @constraint(m, [i in 1:n], sum(y[k,i] for k in 1:K) == 1)

    #Contraintes Xnum et autres dualisées
    @constraint(m, [i in 1:n, j in i+1:n], alpha + beta[i,j] >= (lh[i] + lh[j])*x[i,j])
    @constraint(m, [k in 1:K], W*gamma[k] + sum(w_v[i]*y[k,i] + W_v[i]*zeta[k,i] for i in 1:n) <= B)

    @constraint(m, [i in 1:n, k in 1:K], gamma[k] + zeta[k,i] >= w_v[i]*y[k,i])


    #Objectif
    @objective(m, Min, L*alpha + sum(l[i,j]*x[i,j] + 3*beta[i,j] for i in 1:n,j in i+1:n))

    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", TimeLimit)

    #Résolution
    optimize!(m)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
        computation_time = time() - start

        # Récupération de la valeur de l'objectif
        vOpt = JuMP.objective_value(m)
    end

    return vOpt, computation_time
end