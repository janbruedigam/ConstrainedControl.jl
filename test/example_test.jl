using BenchmarkTools

files = [
    "lqr_acrobot"
    "lqr_cartpole_n_pendulum"
    "lqr_cartpole"
    "lqr_deltabot"
    "lqr_pendulum"
    "lqr_prismatic"
    "lqr_sawyer"
    "pid_doublependulum"
    "pid_pendulum"
    "trackingLQR_triple_cartpole"
]

for file in files
    # println(file)
    include("examples/"*file*".jl")
end