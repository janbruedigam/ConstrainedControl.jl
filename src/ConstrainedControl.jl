module ConstrainedControl

using ConstrainedDynamics
using ConstrainedDynamics: minimalCoordinates, geteqconstraint, Vmat
using StaticArrays

import ControlSystems: dlqr

export PID,
    LQR,

    ilqr

include("pid.jl")
include("lqr.jl")
include("ilqr.jl")

end
