using JuMP
using CPLEX

path = "data 2/10_ulysses_3.tsp"

function static_solving(inputFile::String)

      include(inputFile)

      #definition de l
      l = [sqrt((coordinates[i, 1]-coordinates[j, 1])^2 + (coordinates[i, 2] - coordinates[j, 2])^2) for i in 1:n, j in 1:n]

      m = Model(CPLEX.Optimizer)

      #Variables
      @variable(m, x[i in 1:n, j in i+1:n], Bin)
      @variable(m, y[k in 1:K, i in 1:n], Bin)

      #Contraintes Xcomb
      @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] + y[k,j] <= 1+x[i,j])
      @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], y[k,i] - y[k,j] <= 1-x[i,j])
      @constraint(m, [i in 1:n, j in i+1:n, k in 1:K], -y[k,i] + y[k,j] <= 1-x[i,j])

      @constraint(m, [i in 1:n], sum(y[k,i] for k in 1:K) == 1)

      @constraint(m, [k in 1:K], sum(w_v[i]*y[k,i] for i in 1:n) <= B)

      #Objectif
      @objective(m, Min, sum(l[i,j]*x[i,j] for i in 1:n,j in i+1:n))

      #Résolution
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
      println(m)
end

static_solving(path)