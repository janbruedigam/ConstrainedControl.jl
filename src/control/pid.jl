# TODO
# Only for 1dof joints
mutable struct PID{T,N} <: Controller
    P::SVector{N,T}
    I::SVector{N,T}
    D::SVector{N,T}

    eqcids::SVector{N,Int64}
    goals::SVector{N,T}
    integratederrors::SVector{N,T}
    lasterrors::SVector{N,T}

    control!::Function


    function PID(mechanism::Mechanism{T}, eqcid::Int64, goal; P = zero(T), I = zero(T), D = zero(T), controlfunction::Function = control_pid!) where T
        eqc = geteqconstraint(mechanism, eqcid)
        Nb = 6 * length(unique(eqc.childids))
        Nc = ConstrainedDynamics.length(eqc)
        @assert Nb-Nc == 1 "Only 1 DOF joints are supported"
        
        new{T,1}([P], [I], [D], [eqcid], [goal], [0], [0], controlfunction)
    end

    function PID(mechanism, eqcids::AbstractVector{Int64}, goals::AVec; 
            P::AVec = zeros(T,length(eqcids)), I::AVec = zeros(T,length(eqcids)), D::AVec = zeros(T,length(eqcids)),
            controlfunction::Function = control_pid!
        ) where {T, AVec<:AbstractVector{T}}

        N = length(eqcids)
        for eqcid in eqcids
            eqc = geteqconstraint(mechanism, eqcid)
            Nb = 6 * length(unique(eqc.childids))
            Nc = 0
            Nc = ConstrainedDynamics.length(eqc)
            @assert Nb-Nc == 1 "Only 1 DOF joints are supported"
        end
        new{T,N}(P, I, D, eqcids, goals, zeros(T,N), zeros(T,N), controlfunction)
    end
end

# Angle correction only works for revolute joints
function stateError_pid(mechanism, eqc, goal; anglecorrection = true)
    if anglecorrection && (typeof(eqc.constraints[1]) <: ConstrainedDynamics.Translational3 && typeof(eqc.constraints[2]) <: ConstrainedDynamics.Rotational2) 
        mincoords = minimalCoordinates(mechanism, eqc)[1]
        tempdiff = goal - mincoords
        if tempdiff > pi
            return tempdiff - 2*pi
        elseif tempdiff < -pi
            return tempdiff + 2*pi
        else
            return tempdiff
        end
    else
        return goal - minimalCoordinates(mechanism, eqc)[1]
    end
end

# @generated function error_pid(mechanism, eqcids::SVector{N,Int64}, goals) where {N}
#     vec = [:(goals[$i] - minimalCoordinates(mechanism, geteqconstraint(mechanism, eqcids[$i]))[1]) for i = 1:N]
#     return :(svcat($(vec...)))
# end

function error_pid(mechanism, eqcids::SVector{N,Int64}, goals; anglecorrection = true) where {N}
    vec = [stateError_pid(mechanism, geteqconstraint(mechanism, eqcids[i]), goals[i]; anglecorrection = anglecorrection) for i = 1:N]
    return svcat(vec...)
end

function control_pid!(mechanism, pid::PID{T,N}, k) where {T,N}
    Δt = mechanism.Δt

    currenterrors = error_pid(mechanism, pid.eqcids, pid.goals, anglecorrection = true)
    k==1 && (pid.lasterrors = currenterrors)

    perrors = currenterrors
    pid.integratederrors += currenterrors * Δt
    differentialerrors = (currenterrors - pid.lasterrors) / Δt
    
    u = pid.P .* perrors + pid.I .* pid.integratederrors + pid.D .* differentialerrors

    pid.lasterrors = currenterrors

    for i=1:N
        setForce!(mechanism, geteqconstraint(mechanism, pid.eqcids[i]), SVector(u[i]))
    end

    return
end