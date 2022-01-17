module Tst
include("minimum_jerk.jl")
using Plots, .MinimumJerk

pyplot()
Plots.PyPlotBackend()

export mj_analytical, mj_lp

function mj_analytical()
  initial_state = [0 0 0]'
  final_state = [1 0 0]'
  ts = Vector(LinRange(0, 1, 10))
  traj_type = "position"
  solve_method = "analytical"
  position = get_trajectory(initial_state, final_state, ts, traj_type, solve_method)
  plot(position)
end

function mj_lp()
  initial_state = [0 0 0]'
  final_state = [1 0 0]'
  ts = Vector(LinRange(0, 1, 10))
  traj_type = "position"
  solve_method = "lp"
  position = get_trajectory(initial_state, final_state, ts, traj_type, solve_method)
  plot(position)
end

function mj_dual()
  initial_state = [0 0 0]'
  final_state = [1 0 0]'
  ts = Vector(LinRange(0, 1, 10))
  traj_type = "position"
  solve_method = "dual"
  position = get_trajectory(initial_state, final_state, ts, traj_type, solve_method)
  plot(position)
end

end

