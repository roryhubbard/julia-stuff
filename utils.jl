module Utils
using LinearAlgebra
export evaluate_polynomial, get_matrices


function evaluate_polynomial(t, coefficients, derivative_order, po=5)
  if po > 5
    error("polynomial order > 5")
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
    return
  end
  dot(equation[1:po+1], coefficients)
end


function get_matrices(xi, xf, ts, po=5)
  if po > 5
    error("polynomial order > 5")
  end
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


end

