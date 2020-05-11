# Only for constant reference trajectories
mutable struct LQR{T} <: Controller
    KT::SMatrix{3,12,T,36}
    KR::SMatrix{3,12,T,36}

    xd::SVector{3,Float64}
    vd::SVector{3,Float64}
    Fd::SVector{3,Float64}
    qd::Quaternion{T}
    ωd::SVector{3,Float64}
    τd::SVector{3,Float64}

    bodyid::Int64

    control!::Function

    function LQR(mechanism, bodyid::Int64, Q::AbstractMatrix{T}, R::AbstractMatrix{T}; 
        xd::AbstractVector{T}=zeros(T,3), vd::AbstractVector{T}=zeros(T,3), Fd::AbstractVector{T}=zeros(T,3),
        qd::Quaternion{T}=Quaternion{T}(), ωd::AbstractVector{T}=zeros(T,3), τd::AbstractVector{T}=zeros(T,3)) where T

        Δt = mechanism.Δt
        body = getbody(mechanism, bodyid)
        # linearize        
        Afull, Bfull = ConstrainedDynamics.∂zp1∂z(mechanism, body, xd, vd, Fd, qd, ωd, τd, Δt)
        # A = Afull[1:6,1:6]
        # B = Bfull[1:6,1:3]

        # calculate K
        K = dlqr(Afull,Bfull,Q,R)    
        KT = K[1:3,:]    
        KR = K[4:6,:]
        
        new{T}(KT, KR, xd, vd, Fd, qd, ωd, τd, bodyid, control_lqr!), Afull, Bfull
    end
end

function control_lqr!(mechanism, lqr::LQR{T}) where {T}
    Δt = mechanism.Δt
    body = getbody(mechanism, lqr.bodyid)

    x3 = ConstrainedDynamics.getx3(body,Δt)
    v2 = ConstrainedDynamics.getvnew(body)
    q3 = ConstrainedDynamics.getq3(body,Δt)
    ω2 = ConstrainedDynamics.getωnew(body)

    ΔzpT = [x3-lqr.xd;v2-lqr.vd]
    ΔzpR = [ConstrainedDynamics.VLᵀmat(lqr.qd)*q3;ω2-lqr.ωd]
    Δzp = [ΔzpT;ΔzpR]
    F = -lqr.KT * Δzp
    τ = -lqr.KR * Δzp

    setForce!(mechanism, body, F=F+lqr.Fd, τ=τ+lqr.τd)
end