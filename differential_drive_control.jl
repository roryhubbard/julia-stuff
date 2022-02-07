module DifferentialDriveControl


struct UnicycleKinematicModel
  x
  y
  yaw
  r
end


struct DifferentialDriveKinematicModel
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
  ̇q = G * [v; w]
  q += ̇q* dₜ
  model.x = q[0]
  model.y = q[1]
  model.yaw = q[2]
end


function update!(model::DifferentialDriveKinematicModel, v, w, dₜ=0.01)
  # http://msl.cs.uiuc.edu/planning/node659.html
  r = model.r
  L = model.L
  uₗ, uᵣ = calculate_wheel_velocities(v, w)
  uₜ = (uₗ+ uᵣ) / 2
  uᵩ = uᵣ - uₗ
  G = [r*cos(yaw)   0;
       r*sin(yaw)   0;
                0 r/L]
  ̇q = G * [uₜ; uᵩ]
  q += ̇q * dₜ
  model.x = q[0]
  model.y = q[1]
  model.yaw = q[2]
end


function calculate_wheel_velocities(model::DifferentialDriveKinematicModel, v, w)
  # http://faculty.salina.k-state.edu/tim/robotics_sg/Control/kinematics/unicycle.html
  r = model.r
  L = model.L
  J⁻¹ = [1/r  L/(2r);
         1/r -L/(2r)]
  J⁻¹ * [v; w]
end

