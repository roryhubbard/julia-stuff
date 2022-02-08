module DifferentialDriveControl
using Plots
export test_unicycle_model, test_differential_drive_model


mutable struct UnicycleKinematicModel
  x
  y
  yaw
  r
end


mutable struct DifferentialDriveKinematicModel
  x
  y
  yaw
  r
  L
end


function update!(model::UnicycleKinematicModel, v, w, dₜ=0.01)
  yaw = model.yaw
  G = [cos(yaw) 0;
       sin(yaw) 0;
              0 1]
  qdot = G * [v; w]
  qd = qdot * dₜ
  model.x += qd[1]
  model.y += qd[2]
  model.yaw += qd[3]
end


function update!(model::DifferentialDriveKinematicModel, v, w, dₜ=0.01)
  # http://msl.cs.uiuc.edu/planning/node659.html
  yaw = model.yaw
  r = model.r
  L = model.L
  uᵣ, uₗ = calculate_wheel_velocities(model, v, w)
  uₜ = (uₗ+ uᵣ) / 2
  uᵩ = uᵣ - uₗ
  G = [r*cos(yaw)   0;
       r*sin(yaw)   0;
                0 r/L]
  qdot = G * [uₜ; uᵩ]
  qd = qdot * dₜ
  model.x += qd[1]
  model.y += qd[2]
  model.yaw += qd[3]
end


function calculate_wheel_velocities(model::DifferentialDriveKinematicModel, v, w)
  # http://faculty.salina.k-state.edu/tim/robotics_sg/Control/kinematics/unicycle.html
  r = model.r
  L = model.L
  J⁻¹ = [1/r  L/(2r);
         1/r -L/(2r)]
  J⁻¹ * [v; w]
end


function test_unicycle_model()
  dd = UnicycleKinematicModel(0., 0., 0., .1)
  v = 1.
  w = π / 10
  x = Vector{Float64}()
  y = Vector{Float64}()
  push!(x, dd.x)
  push!(y, dd.y)
  for _ in 1:100
    update!(dd, v, w)
    push!(x, dd.x)
    push!(y, dd.y)
  end
  plot(x, y, aspect_ratio=:equal)
end


function test_differential_drive_model()
  dd = DifferentialDriveKinematicModel(0., 0., 0., .1, .5)
  v = 1.
  w = π / 10
  x = Vector{Float64}()
  y = Vector{Float64}()
  push!(x, dd.x)
  push!(y, dd.y)
  for t in 1:100
    update!(dd, v, w)
    push!(x, dd.x)
    push!(y, dd.y)
  end
  plot(x, y, aspect_ratio=:equal)
end

end

