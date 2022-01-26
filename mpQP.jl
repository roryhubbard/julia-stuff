module mpQP
using Convex, Gurobi, LinearAlgebra
export main


struct ExplicitSolution
  problem::Problem
end


function explicit_solution(problem::Problem)
  0
end

function main()
  t = 1
  Q = [
    0 0 0      0       0      0;
    0 0 0      0       0      0;
    0 0 0      0       0      0;
    0 0 0    36t   72t^2 120t^3;
    0 0 0  72t^2  192t^3 360t^4;
    0 0 0 120t^3  360t^4 720t^5;
  ]
  c = zeros(size(Q)[1])
  A = [
    1 0   0    0     0     0;
    0 1   0    0     0     0;
    1 t t^2  t^3   t^4   t^5;
    0 1  2t 3t^2  4t^3  5t^4;
  ]
  xi = [0 0]
  xf = [1 0]
  b = [xi xf]'
  F = Matrix{Float64}(I, length(b), length(b))
  x = Variable(size(Q)[1])
  θ = Variable(length(b))
  objective = (1/2)*quadform(x, Q) + c'x
  constraints = [A*x <= b + F*θ]
  problem = minimize(objective, constraints)
end


end

