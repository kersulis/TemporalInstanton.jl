time_steps = 10
fixed_wind = []
idx = 1

traj = temperature_trajectory(i, o, idx, time_steps)
# using Plots
# gr(
#     label="",
#     lw=2,
#     color=:black
# )
# tt = range(0; stop=i.time_values[end], length=length(traj)) ./ 60
# plot(tt, traj; xlabel="Time (min)", ylabel="Temperature (C)")

# estimate_temperature_limits(nd)

@testset "Temperature" begin
    atol = 1e-4
    @test traj[5] ≈ 59.5963 atol=atol
    @test traj[end] ≈ 100.0011 atol=atol
end
