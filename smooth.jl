module Smooth
include("utils.jl")
using .Utils
using LinearAlgebra, Convex, Gurobi, Plots
export main

gr()
Plots.GRBackend()


function main()
  xi = [0 0 0]
  xf = [1 0 0]
  ts = Vector(LinRange(0, 1, 10))
  polynomial_order = 5

  P, A, b = get_matrices(xi, xf, ts, polynomial_order)
  x = Variable(size(P)[1])
  objective = quadform(x, P) # x'Px
  constraints = [A * x == b]

  halfway = evaluate_polynomial(ts[end ÷ 2], x, 0)
  constraints += [halfway >= .8]

  problem = minimize(objective, constraints)
  solve!(problem, Gurobi.Optimizer)
  coefficients = evaluate(x)
  plot_all(coefficients, ts)
end


function plot_all(coefficients, ts, po=5)
  position = map(t -> evaluate_polynomial(t, coefficients, 0, po), ts)
  velocity = map(t -> evaluate_polynomial(t, coefficients, 1, po), ts)
  acceleration = map(t -> evaluate_polynomial(t, coefficients, 2, po), ts)
  jerk = map(t -> evaluate_polynomial(t, coefficients, 3, po), ts)
  p1 = plot(position)
  p2 = plot(velocity)
  p3 = plot(acceleration)
  p4 = plot(jerk)
  display(plot(p1, p2, p3, p4, layout=(2, 2)))
end


end

