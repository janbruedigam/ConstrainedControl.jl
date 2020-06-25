using BenchmarkTools

files = [
    "lqr_acrobot"
    "lqr_pendulum"
    "lqr_prismatic"
    "pid_doublependulum"
    "pid_pendulum"
    "sawyer"
]

for file in files
    # println(file)
    include("examples/"*file*".jl")
end