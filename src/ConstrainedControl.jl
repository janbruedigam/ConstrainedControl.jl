module ConstrainedControl

using ConstrainedDynamics
using ConstrainedDynamics: minimalCoordinates, geteqconstraint
using StaticArrays

export PID

include("pid.jl")

end
