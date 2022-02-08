module TrajectoryTracking
using Convex, SCS, Plots
include("wheeled_models.jl")
include("utils.jl")
using .Utils
export differential_flatness_trajectory


function add_boundary_constraints(problem, x, p₀, p₁, v₀, v₁, a₀, a₁, t)
  problem.constraints += [
    evaluate_polynomial(0, x, 0) == p₀,
    evaluate_polynomial(0, x, 1) == v₀,
    #evaluate_polynomial(0, x, 2) == a₀,
    evaluate_polynomial(t, x, 0) == p₁,
    evaluate_polynomial(t, x, 1) == v₁,
    evaluate_polynomial(t, x, 2) == a₁,
  ]
end

function polynomial_trajectory(polynomial_order, p₀, p₁, v₀, v₁, a₀, a₁, θ₀, θ₁, t)
  if polynomial_order > 5
    error("polynomial order > 5")
  end
  P = [
    0 0 0      0       0      0;
    0 0 0      0       0      0;
    0 0 0      0       0      0;
    0 0 0    36t   72t^2 120t^3;
    0 0 0  72t^2  192t^3 360t^4;
    0 0 0 120t^3  360t^4 720t^5;
  ]
  po = polynomial_order
  x_poly = Variable(po+1)
  y_poly = Variable(po+1)
  objective = quadform(x_poly, P[1:po+1, 1:po+1]) + quadform(y_poly, P[1:po+1, 1:po+1])
  problem = minimize(objective)

  add_boundary_constraints(problem, x_poly, p₀[1], p₁[1], v₀[1], v₁[1], a₀[1], a₁[1], t)
  add_boundary_constraints(problem, y_poly, p₀[2], p₁[2], v₀[2], v₁[2], a₀[2], a₁[2], t)
  problem.constraints += [
    evaluate_polynomial(t/2, x_poly, 0) == 1,
  ]

  solve!(problem, SCS.Optimizer)
  evaluate(x_poly), evaluate(y_poly)
end


function recover_θ(t, x_poly, y_poly)
  x_dot = evaluate_polynomial(t, x_poly, 1)
  y_dot = evaluate_polynomial(t, y_poly, 1)
  atan(y_dot, x_dot)
end


function recover_v(t, x_poly, y_poly)
  x_dot = evaluate_polynomial(t, x_poly, 1)
  y_dot = evaluate_polynomial(t, y_poly, 1)
  √(x_dot^2 + y_dot^2)
end


function recover_w(t, x_poly, y_poly)
  x_dot = evaluate_polynomial(t, x_poly, 1)
  y_dot = evaluate_polynomial(t, y_poly, 1)
  x_ddot = evaluate_polynomial(t, x_poly, 2)
  y_ddot = evaluate_polynomial(t, y_poly, 2)
  (y_ddot * x_dot - x_ddot * y_dot) / (x_dot^2 + y_dot^2)
end


function differential_flatness_trajectory()
  poly_order = 5
  p₀ = [0. 0.]
  p₁ = [0. 1.]
  v₀ = [0. 0.]
  v₁ = [0. 0.]
  a₀ = [0. 0.]
  a₁ = [0. 0.]
  θ₀ = 0.
  θ₁ = 0.
  t = 1.
  x_poly, y_poly = polynomial_trajectory(poly_order, p₀, p₁, v₀, v₁, a₀, a₁, θ₀, θ₁, t)
  x = Vector{Float64}()
  y = Vector{Float64}()
  θ = Vector{Float64}()
  for t in LinRange(0, 1, 100)
    x_evaluated = evaluate_polynomial(t, x_poly, 0)
    y_evaluated = evaluate_polynomial(t, y_poly, 0)
    θ_evaluated = recover_θ(t, x_poly, y_poly)
    push!(x, x_evaluated)
    push!(y, y_evaluated)
    push!(θ, θ_evaluated)
  end
  p1 = plot(x, y, aspect_ratio=:equal)
  p2 = plot(θ)
  display(plot(p1, p2, layout=(1, 2)))
end


end
