# Only for constant reference trajectories
mutable struct LQR{T} <: Controller
    K::SMatrix{3,6,T,18}
    xd::SVector{3,Float64}
    vd::SVector{3,Float64}
    ud::SVector{3,Float64}

    bodyid::Int64

    control!::Function

    function LQR(mechanism, bodyid::Int64, Q::AbstractMatrix{T}, R::AbstractMatrix{T}; xd::AbstractVector{T}=zeros(T,3), vd::AbstractVector{T}=zeros(T,3), ud::AbstractVector{T}=zeros(T,3)) where T
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