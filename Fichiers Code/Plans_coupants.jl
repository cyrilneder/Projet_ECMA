using JuMP
using CPLEX

path = "data 2/10_ulysses_3.tsp"

function callbackMain(inputFile::String)

    include(inputFile)

    #definition de l
    l = [sqrt((coordinates[i, 1]-coordinates[j, 1])^2 + (coordinates[i, 2] - coordinates[j, 2])^2) for i in 1:n, j in 1:n]

    m = Model(CPLEX.Optimizer)

    # Il est imposé d'utiliser 1 seul thread en Julia avec CPLEX pour
    # utiliser les callbacks
    MOI.set(m, MOI.NumberOfThreads(), 1)

    #Variables
    @variable(m, z >= 0)
    @variable(m, x[i in 1:n, j in i+1:n], Bin)
    @variable(m, y[k in 1:K, i in 1:n], Bin)

    #Contraintes Xcomb
    @constraint(m, z >= sum(l[i,j]*x[i,j] for i in 1:n,j in i+1:n))

    @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] + y[k,j] <= 1+x[i,j])
    @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] - y[k,j] <= 1-x[i,j])
    @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], -y[k,i] + y[k,j] <= 1-x[i,j])

    @constraint(m, [i in 1:n], sum(y[k,i] for k in 1:K) == 1)

    @constraint(m, [k in 1:K], sum(w_v[i]*y[k,i] for i in 1:n) <= B)

    #Objectif
    @objective(m, Min, z)

    #Définition du Callback
    function mon_super_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if isIntegerPoint(cb_data,context_id)

            CPLEX.load_callback_variable_primal(cb_data, context_id)

            # On récupère la valeur de x, y et de z
            x_val = callback_value.(cb_data, x)
            y_val = callback_value.(cb_data, y)
            z_val = callback_value(cb_data, z)

            #Sous-problème 1
            SP1 = Model(CPLEX.Optimizer)

            @variable(SP1, delta1[i in 1:n, j in i+1:n] >= 0)

            @constraint(SP1, sum(delta1[i,j] for i in 1:n,j in i+1:n) <= L)
            @constraint(SP1, [i in 1:n, j in i+1:n], delta1[i,j] <= 3)

            @objective(SP1, Max, sum((l[i,j] + delta1[i,j]*(lh[i]+lh[j]))*x_val[i,j]  for i in 1:n,j in i+1:n))

            # Désactive les sorties de CPLEX (optionnel)
            set_optimizer_attribute(SP1, "CPX_PARAM_SCRIND", 0)

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

            if z1 > z_val
                cstr = @build_constraint(z >= sum((l[i,j] + delta1_val[i,j]*(lh[i]+lh[j]))*x[i,j]  for i in 1:n,j in i+1:n))
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
                println("Add constraint type SP1 because of objectif breach ", z1)
            end

            #Sous-problème 2,k
            delta2 = Matrix{Float64}(undef,K,n)
            z2 =[]

            for k in 1:K
                SP2k = Model(CPLEX.Optimizer)

                @variable(SP2k, delta2k[i in 1:n] >= 0)

                @constraint(SP2k, sum(delta2k[i] for i in 1:n) <= W)
                @constraint(SP2k, [i in 1:n], delta2k[i] <= W_v[i])

                @objective(SP2k, Max, sum(w_v[i]*(1 + delta2k[i])*y_val[k,i]  for i in 1:n))

                # Désactive les sorties de CPLEX (optionnel)
                set_optimizer_attribute(SP2k, "CPX_PARAM_SCRIND", 0)

                #Résolution
                optimize!(SP2k)

                # Récupération du status de la résolution
                feasibleSolutionFound = primal_status(SP2k) == MOI.FEASIBLE_POINT
                isOptimal = termination_status(SP2k) == MOI.OPTIMAL
                if feasibleSolutionFound
                        # Récupération des valeurs d’une variable
                        delta2[k,:] = JuMP.value.(delta2k)
                        push!(z2, JuMP.objective_value(SP2k))
                end
            end

            if maximum(z2) > B
                for k in 1:K
                    cstr = @build_constraint(sum(w_v[i]*(1 + delta2[k,i])*y[k,i]  for i in 1:n) <= B)
                    MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)

                end
                println("Add constraint type SP2k because of cluster weight breach : ", maximum(z2))
            end
        end
    end

    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)

    # On précise que le modèle doit utiliser notre fonction de callback
    MOI.set(m, CPLEX.CallbackFunction(), mon_super_callback)
    optimize!(m)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
            # Récupération des valeurs d’une variable
            vx = JuMP.value.(x)
            vy = JuMP.value.(y)
            vOpt = JuMP.objective_value(m)
    end

    #Affichage de solution
    println("Valeur optimale :", vOpt," obtenue avec les clusters : ")
    for k in 1:K
        println(findall(!iszero, vy[k,:]))
    end

end

# Fonction  permettant  de  déterminer  si  c'est  l'obtention  d'une
# solution entière qui a entraîné l'appel d'un callback
# (il n'est pas nécessaire d'en comprendre le fonctionnement)
function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)

    # context_id  == CPX_CALLBACKCONTEXT_CANDIDATE si le  callback est
    # appelé dans un des deux cas suivants :
    # cas 1 - une solution entière a été obtenue; ou
    # cas 2 - une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end

    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)

    # S'il n'y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end

callbackMain(path)
