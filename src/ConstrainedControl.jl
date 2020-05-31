module ConstrainedControl

using ConstrainedDynamics
using ConstrainedDynamics: minimalCoordinates, geteqconstraint, Vmat
using LinearAlgebra
using StaticArrays

export PID,
    LQR,

    ilqr


include(joinpath("util","util.jl"))
include(joinpath("control","pid.jl"))
include(joinpath("control","lqr.jl"))
include(joinpath("control","ilqr.jl"))

end
