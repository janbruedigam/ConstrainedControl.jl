module ConstrainedControl

using ConstrainedDynamics
using ConstrainedDynamics: svcat, szeros, imag
using LinearAlgebra
using StaticArrays

export PID,
    LQR,
    TrackingLQR,

    ilqr


include(joinpath("util","util.jl"))
include(joinpath("control","pid.jl"))
include(joinpath("control","lqr.jl"))
include(joinpath("control","lqr_tracking.jl"))
# include(joinpath("control","ilqr.jl"))

end
