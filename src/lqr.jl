# Only for revolute joints
mutable struct LQR{T} <: Controller
    K::SMatrix{3,6,T,18}

    bodyid::Int64

    control!::Function

    function LQR(mechanism, bodyid::Int64, xd::AbstractVector{T}, vd::AbstractVector{T}, ud::AbstractVector{T}, Q::AbstractMatrix{T}, R::AbstractMatrix{T};) where T
        Δt = mechanism.Δt
        body = getbody(mechanism, bodyid)
        # linearize:
        E = SMatrix{3,3,T,9}(1,0,0,0,1,0,0,0,1)
        Z = @SMatrix zeros(T,3,3)
        M = ConstrainedDynamics.getM(body)[1:3,1:3]
        A = [[E E*Δt];[Z E]]
        B = [inv(M)*Δt^2;inv(M)*Δt]
        
        # calculate K
        K = dlqr(A,B,Q,R)        
        
        new{T}(K, bodyid, control_lqr!)
    end
end

function control_lqr!(mechanism, lqr::LQR{T}) where {T}
    Δt = mechanism.Δt
    body = getbody(mechanism, lqr.bodyid)
    zp = [ConstrainedDynamics.getx3(body,Δt);ConstrainedDynamics.getvnew(body)]
    F = -lqr.K * zp

    setForce!(mechanism, body, F=F)
end