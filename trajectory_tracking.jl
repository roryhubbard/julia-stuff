module TrajectoryTracking
using Convex, SCS, Plots
include("wheeled_models.jl")
include("utils.jl")
using .Utils
export differential_flatness_trajectory


function polynomial_trajectory(polynomial_order, p₀, p₁, v₀, v₁, θ₀, θ₁, t, is_x=false)
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
  x = Variable(po+1)
  objective = quadform(x, P[1:po+1, 1:po+1])
  constraints = [
    evaluate_polynomial(0, x, 0, po) == p₀,
    evaluate_polynomial(0, x, 1, po) == v₀,
    evaluate_polynomial(t, x, 0, po) == p₁,
    evaluate_polynomial(t, x, 1, po) == v₁,
  ]
  if is_x
    constraints += [
      evaluate_polynomial(0, x, 1, po) == v₀*cos(θ₀),
      evaluate_polynomial(t, x, 1, po) == v₁*cos(θ₁),
    ]
  else
    constraints += [
      evaluate_polynomial(0, x, 1, po) == v₀*sin(θ₀),
      evaluate_polynomial(t, x, 1, po) == v₁*sin(θ₁),
    ]
  end
  problem = minimize(objective, constraints)
  solve!(problem, SCS.Optimizer)
  evaluate(x)
end


function differential_flatness_trajectory()
  poly_order = 3
  x_poly = polynomial_trajectory(poly_order, 0, 0, 0, 0, 0, 0, 1, true)
  y_poly = polynomial_trajectory(poly_order, 0, 1, 0, 0, 0, 0, 1, false)
  x = Vector{Float64}()
  y = Vector{Float64}()
  θ = Vector{Float64}()
  for t in LinRange(0, 1, 100)
    x_evaluated = evaluate_polynomial(t, x_poly, 0, poly_order)
    y_evaluated = evaluate_polynomial(t, y_poly, 0, poly_order)
    push!(x, x_evaluated)
    push!(y, y_evaluated)
  end
  plot(x, y, aspect_ratio=:equal)
  #p1 = plot(x, y, aspect_ratio=:equal)
  #p2 = plot(θ)
  #display(plot(p1, p2, layout=(1, 2)))
end


end
