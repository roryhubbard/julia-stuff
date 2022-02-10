module WheeledModels
using Plots
export UnicycleKinematicModel, DifferentialDriveKinematicModel, update!


mutable struct UnicycleKinematicModel
  x
  y
  θ
end


mutable struct DifferentialDriveKinematicModel
  x
  y
  θ
  r
  L
end


function update!(model::UnicycleKinematicModel, v, w, dₜ=0.1)
  θ = model.θ
  G = [cos(θ) 0;
       sin(θ) 0;
            0 1]
  qdot = G * [v; w]
  qd = qdot * dₜ
  model.x += qd[1]
  model.y += qd[2]
  model.θ += qd[3]
  model.θ = atan(sin(model.θ), cos(model.θ))
end


function update!(model::DifferentialDriveKinematicModel, v, w, dₜ=0.1)
  # http://msl.cs.uiuc.edu/planning/node659.html
  θ = model.θ
  r = model.r
  L = model.L
  uᵣ, uₗ = calculate_wheel_velocities(model, v, w)
  uₜ = (uₗ+ uᵣ) / 2
  uᵩ = uᵣ - uₗ
  G = [r*cos(θ)   0;
       r*sin(θ)   0;
              0 r/L]
  qdot = G * [uₜ; uᵩ]
  qd = qdot * dₜ
  model.x += qd[1]
  model.y += qd[2]
  model.θ += qd[3]
  model.θ = atan(sin(model.θ), cos(model.θ))
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
  um = UnicycleKinematicModel(0., 0., 0.)
  v = 1.
  w = π
  x = Vector{Float64}()
  y = Vector{Float64}()
  θ = Vector{Float64}()
  push!(x, um.x)
  push!(y, um.y)
  push!(θ, um.θ)
  for _ in 1:100
    update!(um, v, w)
    push!(x, um.x)
    push!(y, um.y)
    push!(θ, um.θ)
  end
  p1 = plot(x, y, aspect_ratio=:equal)
  p2 = plot(θ)
  display(plot(p1, p2, layout=(1, 2)))
end


function test_differential_drive_model()
  ddm = DifferentialDriveKinematicModel(0., 0., 0., .1, .5)
  v = 1.
  w = π
  x = Vector{Float64}()
  y = Vector{Float64}()
  θ = Vector{Float64}()
  push!(x, ddm.x)
  push!(y, ddm.y)
  push!(θ, ddm.θ)
  for _ in 1:100
    update!(ddm, v, w)
    push!(x, ddm.x)
    push!(y, ddm.y)
    push!(θ, ddm.θ)
  end
  p1 = plot(x, y, aspect_ratio=:equal)
  p2 = plot(θ)
  display(plot(p1, p2, layout=(1, 2)))
end


end

