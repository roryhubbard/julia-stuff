module TrajectoryTracking
using LinearAlgebra, ControlSystems, Plots
include("wheeled_models.jl")
using .WheeledModels
include("trajectory_planning.jl")
using .TrajectoryPlanning
export track_trajectory


function express_in_rotated_frame(x, θ)
  R = [ cos(θ) sin(θ) 0;
       -sin(θ) cos(θ) 0;
             0      0 1]
  R * x
end


function get_control_input(xᵣ, yᵣ, θᵣ, vᵣ, wᵣ, q, t)
  A = [    1 wᵣ*t    0;
       -wᵣ*t    1 vᵣ*t;
           0    0    1]
  B = [-t  0;
        0  0;
        0 -t]
  Q = I
  R = I
  K = dlqr(A, B, Q, R)
  qᵣ = [xᵣ yᵣ θᵣ]'
  e = express_in_rotated_frame(qᵣ - q, q[3])
  u₁, u₂ = -K * e
  v = vᵣ * cos(e[3]) + u₁
  w = wᵣ + u₂
  v, w
end


function track_trajectory()
  xᵣ, yᵣ, θᵣ, vᵣ, wᵣ, t = differential_flatness_trajectory()
  ddm = DifferentialDriveKinematicModel(xᵣ[2], yᵣ[2], θᵣ[2], 0.1, 0.5)
  x = Vector{Float64}()
  y = Vector{Float64}()
  θ = Vector{Float64}()
  push!(x, ddm.x)
  push!(y, ddm.y)
  push!(θ, ddm.θ)
  for i in 2:length(xᵣ)
    q = [ddm.x ddm.y ddm.θ]'
    v, w = get_control_input(xᵣ[i], yᵣ[i], θᵣ[i], vᵣ[i], wᵣ[i], q, t)
    update!(ddm, v, w, t)
    #update!(ddm, vᵣ[i], wᵣ[i], t) # open loop
    push!(x, ddm.x)
    push!(y, ddm.y)
    push!(θ, ddm.θ)
  end
  p1 = plot(xᵣ, yᵣ, aspect_ratio=:equal)
  plot!(p1, x, y)
  p2 = plot(θᵣ)
  plot!(p2, θ)
  display(plot(p1, p2, layout=(1, 2)))
end

end
