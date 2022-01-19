module MinimumJerk
using LinearAlgebra, Convex, Gurobi, Plots
export test_analytical, test_primal

pyplot()
Plots.PyPlotBackend()


function eval_traj_point(t, coefficients, derivative_order, po=5)
  if po > 5
    println("polynomial order > 5")
    return
  end
  if derivative_order == 0
    equation = [1 t t^2 t^3 t^4 t^5]
  elseif derivative_order == 1
    equation = [0 1 2t 3t^2 4t^3 5t^4]
  elseif derivative_order == 2
    equation = [0 0 2 6t 12t^2 20t^3]
  elseif derivative_order == 3
    equation = [0 0 0 6 24t 60t^2]
  elseif derivative_order == 4
    equation = [0 0 0 0 24 120t]
  elseif derivative_order == 5
    equation = [0 0 0 0 0 120]
  elseif derivative_order == 6
    equation = [0 0 0 0 0 0]
  else
    println("derivative order is too high")
    return 0
  end
  dot(equation[1:po+1], coefficients)
end


function get_matrices(xi, xf, ts, po=5)
  if po > 5
    println("polynomial order > 5")
    return
  end
  # po = polynomial order
  ti = first(ts)
  tf = last(ts)
  P = [
    0 0 0              0               0              0;
    0 0 0              0               0              0;
    0 0 0              0               0              0;
    0 0 0      36(tf-ti)   72(tf^2-ti^2) 120(tf^3-ti^3);
    0 0 0  72(tf^2-ti^2)  192(tf^3-ti^3) 360(tf^4-ti^4);
    0 0 0 120(tf^3-ti^3)  360(tf^4-ti^4) 720(tf^5-ti^5);
  ]
  A = [
    1 ti ti^2  ti^3   ti^4   ti^5;
    0  1  2ti 3ti^2  4ti^3  5ti^4;
    0  0    2   6ti 12ti^2 20ti^3;
    1 tf tf^2  tf^3   tf^4   tf^5;
    0  1  2tf 3tf^2  4tf^3  5tf^4;
    0  0    2   6tf 12tf^2 20tf^3;
  ]
  b = collect(skipmissing([xi xf]))
  A_row_idx = findall(!ismissing, vec([xi xf]))
  P[1:po+1, 1:po+1], A[A_row_idx, 1:po+1], b
end


function solve_primal(P, A, b)
  x = Variable(size(P)[1])
  objective = quadform(x, P) # 1/2x'Px
  constraints = [A * x == b]
  problem = minimize(objective, constraints)
  solve!(problem, Gurobi.Optimizer)
  evaluate(x)
end


function test_analytical()
  xi = [0 0 0]
  xf = [1 0 0]
  ts = Vector(LinRange(0, 1, 10))
  derivative_order = 0
  polynomial_order = 5
  _P, A, b = get_matrices(xi, xf, ts, polynomial_order)
  coefficients = inv(A) * b
  position = map(t -> eval_traj_point(t, coefficients, derivative_order, polynomial_order), ts)
  plot(position)
end


function test_primal()
  xi = [0 0 0]
  xf = [1 missing 0]
  ts = Vector(LinRange(0, 1, 10))
  derivative_order = 0
  polynomial_order = 5
  P, A, b = get_matrices(xi, xf, ts, polynomial_order)
  coefficients = solve_primal(P, A, b)
  position = map(t -> eval_traj_point(t, coefficients, derivative_order, polynomial_order), ts)
  plot(position)
end


end

