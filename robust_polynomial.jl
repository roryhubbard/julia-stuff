module RobustPolynomial
using Convex, Gurobi, LinearAlgebra, Plots
export get_trajectory, test

pyplot()
Plots.PyPlotBackend()


function get_trajectory(xi, xf, ts, traj_type="position")
  coefficients = solve_for_coefficients(xi, xf, ts)
  map(t -> eval_traj_point(t, coefficients, traj_type), ts)
end


function eval_traj_point(t, coefficients, traj_type)
  if traj_type == "position"
    equation = [1 t t^2 t^3 t^4 t^5]
  elseif traj_type == "velocity"
    equation = [0 1 2*t 3*t^2 4*t^3 5*t^4]
  else
    equation = [0 0 2 6*t 12*t^2 20*t^3]
  end
  first(equation * coefficients)
end


function solve_for_coefficients(xi, xf, ts)
  ti = first(ts)
  tf = last(ts)
  c = [0 0 0 0 0 0 6(tf-ti) -6(tf-ti) 12(tf^2-ti^2) -12(tf^2-ti^2) 20(tf^3-ti^3) -20(tf^3-ti^3)]'
  A = [
    1 -1 ti -ti ti^2 -ti^2  ti^3  -ti^3   ti^4   -ti^4   ti^5   -ti^5;
    0  0  1  -1  2ti  -2ti 3ti^2 -3ti^2  4ti^3  -4ti^3  5ti^4  -5ti^4;
    0  0  0   0    2    -2   6ti   -6ti 12ti^2 -12ti^2 20ti^3 -20ti^3;
    1 -1 tf -tf tf^2 -tf^2  tf^3  -tf^3   tf^4   -tf^4   tf^5   -tf^5;
    0  0  1  -1  2tf  -2tf 3tf^2 -3tf^2  4tf^3  -4tf^3  5tf^4  -5tf^4;
    0  0  0   0    2    -2   6tf   -6tf 12tf^2 -12tf^2 20tf^3 -20tf^3;
  ]
  b = [xi xf]'
  x = Variable(12)
  objective = c' * x
  constraints = [A * x == b,
                     x >= 0]
  problem = minimize(objective, constraints)
  solve!(problem, Gurobi.Optimizer)
  x = evaluate(x)
  [ x[1] -  x[2],
    x[3] -  x[4],
    x[5] -  x[6],
    x[7] -  x[8],
    x[9] - x[10],
   x[11] - x[12],
  ]
end


function robust_b_analysis(xi, xf, ts, uncertain_idx)
  ti = first(ts)
  tf = last(ts)
  A = [
    1 -1 ti -ti ti^2 -ti^2  ti^3  -ti^3   ti^4   -ti^4   ti^5   -ti^5;
    0  0  1  -1  2ti  -2ti 3ti^2 -3ti^2  4ti^3  -4ti^3  5ti^4  -5ti^4;
    0  0  0   0    2    -2   6ti   -6ti 12ti^2 -12ti^2 20ti^3 -20ti^3;
    1 -1 tf -tf tf^2 -tf^2  tf^3  -tf^3   tf^4   -tf^4   tf^5   -tf^5;
    0  0  1  -1  2tf  -2tf 3tf^2 -3tf^2  4tf^3  -4tf^3  5tf^4  -5tf^4;
    0  0  0   0    2    -2   6tf   -6tf 12tf^2 -12tf^2 20tf^3 -20tf^3;
  ]
  b = [xi xf]'
  basic_idx = [2 3 6 7 10 11]
  B = reduce(hcat, [A[:, i] for i in basic_idx])
  basic_variables = inv(B) * b
  delta_bounds = get_delta_bounds(B[:, uncertain_idx], basic_variables[uncertain_idx])
  x = zeros(12)
  x[basic_idx] = basic_variables
  [ x[1] -  x[2],
    x[3] -  x[4],
    x[5] -  x[6],
    x[7] -  x[8],
    x[9] - x[10],
   x[11] - x[12],
  ]
end


function get_delta_bounds(g, x)
end


function test()
  initial_state = [0 0 0]
  final_state = [1 0 0]
  ts = Vector(LinRange(0, 1, 10))
  traj_type = "position"
  position = get_trajectory(initial_state, final_state, ts, traj_type)
  plot(position)
end


function test_robust()
  initial_state = [0 0 0]
  final_state = [1 0 0]
  ts = Vector(LinRange(0, 1, 10))
  uncertain_idx = 4
  robust_b_analysis(initial_state, final_state, ts, uncertain_idx)
end


end

