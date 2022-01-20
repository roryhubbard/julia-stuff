module RobustPolynomial
include("utils.jl")
using .Utils
using LinearAlgebra, Convex, SCS, Plots

pyplot()
Plots.PyPlotBackend()


function plot_trajectories(coefficients, ts)
  p1 = plot(map(t -> eval_traj_point(t, coefficients, 0), ts))
  p2 = plot(map(t -> eval_traj_point(t, coefficients, 1), ts))
  p3 = plot(map(t -> eval_traj_point(t, coefficients, 2), ts))
  p4 = plot(map(t -> eval_traj_point(t, coefficients, 3), ts))
  plot(p1, p2, p3, p4, layout = (2, 2), legend = false)
end


function get_matrices(xi, xf, ts)
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
  c, A, b
end


function solve_for_coefficients(c, A, b)
  x = Variable(12)
  objective = c' * x
  constraints = [A * x == b,
                     x >= 0]
  problem = minimize(objective, constraints)
  solve!(problem, SCS.Optimizer)
  x = evaluate(x)
  basic_indices = findall(!iszero, x)
  println(basic_indices)
  coefficients = [
     x[1] -  x[2],
     x[3] -  x[4],
     x[5] -  x[6],
     x[7] -  x[8],
     x[9] - x[10],
    x[11] - x[12],
  ]
  println("optimal cost: ", c' * x)
  coefficients, problem.optval
end


function robust_b_analysis(xi, xf, ts, uncertain_idx)
  c, A, b = get_matrices(xi, xf, ts)
  basic_idx = [2 3 6 7 10 11]
  B = reduce(hcat, [A[:, i] for i in basic_idx])
  basic_variables = inv(B) * b
  delta_bounds = get_delta_bounds(B[:, uncertain_idx], basic_variables)
  x = zeros(12)
  x[basic_idx] = basic_variables
  coefficients = [
     x[1] -  x[2],
     x[3] -  x[4],
     x[5] -  x[6],
     x[7] -  x[8],
     x[9] - x[10],
    x[11] - x[12],
  ]
  println("optimal cost: ", c' * x)
  coefficients, delta_bounds
end


function get_delta_bounds(g, x)
  ratio = -x ./ g
  lower_bound = (any(ratio .< 0) ? maximum(ratio[ratio .< 0]) : -Inf)
  upper_bound = (any(ratio .> 0) ? minimum(ratio[ratio .> 0]) : Inf)
  [lower_bound upper_bound]
end


function test()
  xi = [0 0 0]
  xf = [1 0 0]
  ts = Vector(LinRange(1, 2, 20))
  c, A, b = get_matrices(xi, xf, ts)
  coefficients, optimal_cost = solve_for_coefficients(c, A, b)
  plot_trajectories(coefficients, ts)
end


function test_robust()
  xi = [0 0 0]
  xf = [1 0 0]
  ts = Vector(LinRange(1, 2, 20))
  uncertain_idx = 1
  coefficients, delta_bounds = robust_b_analysis(xi, xf, ts, uncertain_idx)
  println(delta_bounds)
  plot_trajectories(coefficients, ts)
end


end
