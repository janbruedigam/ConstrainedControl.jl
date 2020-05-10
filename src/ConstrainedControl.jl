module ConstrainedControl

using ConstrainedDynamics
using ConstrainedDynamics: minimalCoordinates, geteqconstraint
using StaticArrays

import ControlSystems: dlqr

export PID,
    LQR

include("pid.jl")
include("lqr.jl")

end
