mutable struct PID{T} <: Controller
    P::T
    I::T
    D::T

    eqcid::Int64
    goal::T
    integratedError::T

    function PID(eqcid::Int64, goal::T; P::T = zero(T), I::T = zero(T), D::T = zero(T)) where T
        @assert !(P == 0 && I == 0 && D == 0)

        new{T}(P, I, D, eqcid, goal, 0)
    end
end

@inline stateError(mechanism, eqc, goal, N) = goal - minimalCoordinates(eqc, mechanism, K = N)[1]

function control!(mechanism, pid::PID)
    e0 = stateError(mechanism, geteqconstraint(mechanism, pid.eqcid), pid.goal, 1)
    e1 = stateError(mechanism, geteqconstraint(mechanism, pid.eqcid), pid.goal, 2)

    pid.integratedError += e1 * mechanism.Δt
    differentialError = (e1 - e0) / mechanism.Δt
    perror = e1

    # return pid.P*perror + pid.I*pid.integratedError + pid.D*differentialError
    τ = SVector{3,Float64}(1, 0, 0) * (pid.P * perror + pid.I * pid.integratedError + pid.D * differentialError)
    setForce!(mechanism, mechanism.bodies[1], τ = τ)
end