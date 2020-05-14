# Only for revolute joints
mutable struct PID{T,N} <: Controller
    P::SVector{N,T}
    I::SVector{N,T}
    D::SVector{N,T}

    eqcid::SVector{N,Int64}
    goal::SVector{N,T}
    integratedError::SVector{N,T}

    control!::Function

    function PID(mechanism, eqcid::Int64, goal::T; P::T = zero(T), I::T = zero(T), D::T = zero(T)) where T
        @assert typeof(geteqconstraint(mechanism, eqcid)) <: EqualityConstraint{T,5,2,Tuple{ConstrainedDynamics.Translational3{T},ConstrainedDynamics.Rotational2{T}}} "Only revolute joint supported"
        new{T,1}([P], [I], [D], [eqcid], [goal], [0], control_pid!)
    end
    function PID(mechanism, eqcid::AbstractVector{Int64}, goal::AbstractVector{T}; 
        P::AbstractVector{T} = zeros(T,length(eqcid)), I::AbstractVector{T} = zeros(T,length(eqcid)), D::AbstractVector{T} = zeros(T,length(eqcid))) where T

        N = length(eqcid)
        for i=1:N
            @assert typeof(geteqconstraint(mechanism, eqcid[i])) <: EqualityConstraint{T,5,2,Tuple{ConstrainedDynamics.Translational3{T},ConstrainedDynamics.Rotational2{T}}} "Only revolute joint supported"
        end
        new{T,N}(P, I, D, eqcid, goal, zeros(T,N), control_pid!)
    end
end

@inline stateError_pid(mechanism, eqc, goal, N) = goal - minimalCoordinates(mechanism, eqc, K = N)[1]

function control_pid!(mechanism, pid::PID{T,N}, k) where {T,N}
    e0 = zeros(T,N)
    e1 = zeros(T,N)
    for i=1:N
        e0[i] = stateError_pid(mechanism, geteqconstraint(mechanism, pid.eqcid[i]), pid.goal[i], 1)
        e1[i] = stateError_pid(mechanism, geteqconstraint(mechanism, pid.eqcid[i]), pid.goal[i], 2)
    end

    pid.integratedError += e1 * mechanism.Δt
    differentialError = (e1 - e0) / mechanism.Δt
    perror = e1
    u = pid.P .* perror + pid.I .* pid.integratedError + pid.D .* differentialError

    for i=1:N
        setForce!(mechanism, geteqconstraint(mechanism, pid.eqcid[i]), SVector{1,T}(u[i]))
    end
end