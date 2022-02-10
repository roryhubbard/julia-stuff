module TrajectoryPlanning
using Convex, SCS, Plots
include("utils.jl")
using .Utils
export differential_flatness_trajectory


function add_continuity_constraints!(problem, poly1, poly2, t)
  problem.constraints += [
    evaluate_polynomial(t, poly1, 0) == evaluate_polynomial(0, poly2, 0),
    evaluate_polynomial(t, poly1, 1) == evaluate_polynomial(0, poly2, 1),
    evaluate_polynomial(t, poly1, 2) == evaluate_polynomial(0, poly2, 2),
    evaluate_polynomial(t, poly1, 3) == evaluate_polynomial(0, poly2, 3),
    evaluate_polynomial(t, poly1, 4) == evaluate_polynomial(0, poly2, 4),
  ]
end


function polynomial_trajectory(po, t)
  if po > 5
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
  xpoly1 = Variable(po+1)
  xpoly2 = Variable(po+1)
  ypoly1 = Variable(po+1)
  ypoly2 = Variable(po+1)
  objective = (quadform(xpoly1, P[1:po+1, 1:po+1])
               + quadform(ypoly1, P[1:po+1, 1:po+1])
               + quadform(xpoly2, P[1:po+1, 1:po+1])
               + quadform(ypoly2, P[1:po+1, 1:po+1]))
  problem = minimize(objective)

  problem.constraints += [
    evaluate_polynomial(0., xpoly1, 0) == 0.,
    evaluate_polynomial(t, xpoly1, 0) == 1.,
    evaluate_polynomial(0., xpoly1, 1) == 0.,
    evaluate_polynomial(0., xpoly1, 2) == 0.,
    evaluate_polynomial(t, xpoly2, 0) == 0.,
    evaluate_polynomial(t, xpoly2, 1) == 0.,
    evaluate_polynomial(t, xpoly2, 2) == 0.,
    evaluate_polynomial(0., ypoly1, 0) == 0.,
    evaluate_polynomial(0., ypoly1, 1) == 0.,
    evaluate_polynomial(0., ypoly1, 2) == 0.,
    evaluate_polynomial(t, ypoly2, 0) == 1.,
    evaluate_polynomial(t, ypoly2, 1) == 0.,
    evaluate_polynomial(t, ypoly2, 2) == 0.,
  ]

  add_continuity_constraints!(problem, xpoly1, xpoly2, t)
  add_continuity_constraints!(problem, ypoly1, ypoly2, t)

  solve!(problem, SCS.Optimizer)
  evaluate(xpoly1), evaluate(ypoly1), evaluate(xpoly2), evaluate(ypoly2)
end


function recover_θ(t, xpoly, ypoly)
  x_dot = evaluate_polynomial(t, xpoly, 1)
  y_dot = evaluate_polynomial(t, ypoly, 1)
  atan(y_dot, x_dot)
end


function recover_v(t, xpoly, ypoly)
  x_dot = evaluate_polynomial(t, xpoly, 1)
  y_dot = evaluate_polynomial(t, ypoly, 1)
  √(x_dot^2 + y_dot^2)
end


function recover_w(t, xpoly, ypoly)
  x_dot = evaluate_polynomial(t, xpoly, 1)
  y_dot = evaluate_polynomial(t, ypoly, 1)
  x_ddot = evaluate_polynomial(t, xpoly, 2)
  y_ddot = evaluate_polynomial(t, ypoly, 2)
  (y_ddot * x_dot - x_ddot * y_dot) / (x_dot^2 + y_dot^2)
end


function differential_flatness_trajectory()
  poly_order = 5
  tp = 1.
  xpoly1, ypoly1, xpoly2, ypoly2 = polynomial_trajectory(poly_order, tp)
  x = Vector{Float64}()
  y = Vector{Float64}()
  θ = Vector{Float64}()
  v = Vector{Float64}()
  w = Vector{Float64}()
  for t in LinRange(0, 2tp, 21)
    (xpoly, ypoly) = t < tp ? (xpoly1, ypoly1) : (xpoly2, ypoly2)
    t = mod(t, tp)
    x_evaluated = evaluate_polynomial(t, xpoly, 0)
    y_evaluated = evaluate_polynomial(t, ypoly, 0)
    θ_evaluated = recover_θ(t, xpoly, ypoly)
    v_evaluated = recover_v(t, xpoly, ypoly)
    w_evaluated = recover_w(t, xpoly, ypoly)
    push!(x, x_evaluated)
    push!(y, y_evaluated)
    push!(θ, θ_evaluated)
    push!(v, v_evaluated)
    push!(w, w_evaluated)
  end
  dt = 2tp / 20
  x, y, θ, v, w, dt
end


function test_differential_flatness_trajectory()
  x, y, θ, v, w, _ = differential_flatness_trajectory()
  p1 = plot(x, y, aspect_ratio=:equal)
  p2 = plot(θ)
  display(plot(p1, p2, layout=(1, 2)))

  #p3 = plot(v)
  #p4 = plot(w)
  #display(plot(p3, p4, layout=(1, 2)))
end


end
