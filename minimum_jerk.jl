module MinimumJerk
export get_trajectory

function get_trajectory(xi, xf, ts, traj_type="position", solve_method="analytical")
  if solve_method == "analytical"
    coefficients = MJAnalytical.solve_for_coefficients(xi, xf, ts)
  elseif solve_method == "lp"
    coefficients = MJLinearProgram.solve_for_coefficients(xi, xf, ts)
  elseif solve_method == "dual"
    coefficients = MJDual.solve_for_coefficients(xi, xf, ts)
  else
    print("unrecognized solve method")
    return 0
  end
  # [eval_traj_point(t, coefficients, traj_type) for t in ts]
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


module MJAnalytical
  using LinearAlgebra

  function solve_for_coefficients(xi, xf, t)
    ti = first(t) # t[begin]
    tf = last(t) # t[end]
    b = [xi; xf]
    A = [
      1 ti ti^2 ti^3   ti^4    ti^5;
      0 1  2*ti 3*ti^2 4*ti^3  5*ti^4;
      0 0  2    6*ti   12*ti^2 20*ti^3;
      1 tf tf^2 tf^3   tf^4    tf^5;
      0 1  2*tf 3*tf^2 4*tf^3  5*tf^4;
      0 0  2    6*tf   12*tf^2 20*tf^3;
    ]
    inv(A) * b
  end

end


module MJLinearProgram
  using Convex, Gurobi

  function solve_for_coefficients(xi, xf, t)
    x = Variable(6)
    ti = first(t) # t[begin]
    tf = last(t) # t[end]
    c = [0 0 0 6*(tf-ti) 12*(tf^2-ti^2) 20*(tf^3-ti^3)]'
    A = [
      1 ti ti^2 ti^3   ti^4    ti^5;
      0 1  2*ti 3*ti^2 4*ti^3  5*ti^4;
      0 0  2    6*ti   12*ti^2 20*ti^3;
      1 tf tf^2 tf^3   tf^4    tf^5;
      0 1  2*tf 3*tf^2 4*tf^3  5*tf^4;
      0 0  2    6*tf   12*tf^2 20*tf^3;
    ]
    b = [xi; xf]
    objective = c' * x
    constraints = [A * x == b]
    problem = minimize(objective, constraints)
    solve!(problem, Gurobi.Optimizer)
    evaluate(x)
  end

end


module MJDual
  using Convex, Gurobi, LinearAlgebra

  function solve_for_coefficients(xi, xf, t)
    # TODO: doesn't work yet
    v = Variable(6)
    ti = first(t) # t[begin]
    tf = last(t) # t[end]
    c = [0 0 0 6*(tf-ti) 12*(tf^2-ti^2) 20*(tf^3-ti^3)]'
    A = [
      1 ti ti^2 ti^3   ti^4    ti^5;
      0 1  2*ti 3*ti^2 4*ti^3  5*ti^4;
      0 0  2    6*ti   12*ti^2 20*ti^3;
      1 tf tf^2 tf^3   tf^4    tf^5;
      0 1  2*tf 3*tf^2 4*tf^3  5*tf^4;
      0 0  2    6*tf   12*tf^2 20*tf^3;
    ]
    b = [xi; xf]
    objective = -v' * b
    constraints = [A' * v == -c]
    problem = maximize(objective, constraints)
    solve!(problem, Gurobi.Optimizer)
    evaluate(v)
  end

end


end

