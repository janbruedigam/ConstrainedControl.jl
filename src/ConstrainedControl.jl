module ConstrainedControl

using ConstrainedDynamics
using ConstrainedDynamics: svcat
using LinearAlgebra
using StaticArrays

export PID,
    LQR,

    ilqr


include(joinpath("util","util.jl"))
include(joinpath("control","pid.jl"))
include(joinpath("control","lqr.jl"))
# include(joinpath("control","ilqr.jl"))

end
