# Only for constant reference trajectories
mutable struct LQR{T} <: Controller
    K::SMatrix{3,6,T,18}
    xd::SVector{3,Float64}
    vd::SVector{3,Float64}
    ud::SVector{3,Float64}

    bodyid::Int64

    control!::Function

    function LQR(mechanism, bodyid::Int64, Q::AbstractMatrix{T}, R::AbstractMatrix{T}; xd::AbstractVector{T}=zeros(T,3), vd::AbstractVector{T}=zeros(T,3), Fd::AbstractVector{T}=zeros(T,3)) where T
        Δt = mechanism.Δt
        body = getbody(mechanism, bodyid)
        # linearize        
        Afull, Bfull = ConstrainedDynamics.∂zp1∂z(mechanism, body, xd, vd, Fd, Quaternion{T}(), zeros(T,3), zeros(T,3), Δt)
        A = Afull[1:6,1:6]
        B = Bfull[1:6,1:3]

        # calculate K
        K = dlqr(A,B,Q,R)        
        
        new{T}(K, xd, vd, ud, bodyid, control_lqr!)
    end
end

function control_lqr!(mechanism, lqr::LQR{T}) where {T}
    Δt = mechanism.Δt
    body = getbody(mechanism, lqr.bodyid)
    zp = [ConstrainedDynamics.getx3(body,Δt);ConstrainedDynamics.getvnew(body)]
    F = -lqr.K * (zp - [lqr.xd;lqr.vd])

    setForce!(mechanism, body, F=F+lqr.ud)
end