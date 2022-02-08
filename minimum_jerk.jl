module MinimumJerk
include("utils.jl")
using .Utils
using LinearAlgebra, Convex, SCS, Plots
export test_analytical, test_primal, plot_coefficients


function solve_primal(P, A, b)
  x = Variable(size(P)[1])
  objective = quadform(x, P) # x'Px
  constraints = [A * x == b]
  problem = minimize(objective, constraints)
  solve!(problem, SCS.Optimizer)
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
  position = map(t -> evaluate_polynomial(t, coefficients, derivative_order), ts)
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
  position = map(t -> evaluate_polynomial(t, coefficients, derivative_order), ts)
  plot(position)
end


function plot_coefficients()
  ts = Vector(LinRange(0, 1, 10))
  derivative_order = 0
  polynomial_order = 5
  n = 10
  coefficients = Vector{Float64}[]
  for v in LinRange(-1, 1, n)
    xi = [0 0 0]
    xf = [1 v 0]
    _P, A, b = get_matrices(xi, xf, ts, polynomial_order)
    push!(coefficients, inv(A) * b)
  end
  coefficients = reduce(hcat, coefficients)
  p = plot()
  for i in 1 : polynomial_order + 1
    plot!(p, coefficients[i, :])
  end
  display(p)
end


function plot_all(coefficients, ts)
  position = map(t -> evaluate_polynomial(t, coefficients, 0), ts)
  velocity = map(t -> evaluate_polynomial(t, coefficients, 1), ts)
  acceleration = map(t -> evaluate_polynomial(t, coefficients, 2), ts)
  jerk = map(t -> evaluate_polynomial(t, coefficients, 3), ts)
  p1 = plot(position)
  p2 = plot(velocity)
  p3 = plot(acceleration)
  p4 = plot(jerk)
  display(plot(p1, p2, p3, p4, layout=(2, 2)))
end


end

