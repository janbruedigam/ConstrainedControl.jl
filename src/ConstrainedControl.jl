module ConstrainedControl

using ConstrainedDynamics
using ConstrainedDynamics: svcat, szeros
using LinearAlgebra
using StaticArrays
using Rotations

export PID,
    LQR,

    ilqr


include(joinpath("util","util.jl"))
include(joinpath("control","pid.jl"))
include(joinpath("control","lqr.jl"))
# include(joinpath("control","ilqr.jl"))

end
