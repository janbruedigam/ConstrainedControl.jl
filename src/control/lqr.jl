# Only for constant reference trajectories
mutable struct LQR{T,N} <: Controller
    KT::Vector{SMatrix{3,12,T,36}}
    KR::Vector{SMatrix{3,12,T,36}}

    xd::SVector{3,Float64}
    vd::SVector{3,Float64}
    Fd::SVector{3,Float64}
    qd::Quaternion{T}
    ωd::SVector{3,Float64}
    τd::SVector{3,Float64}

    bodyid::Int64

    control!::Function

    function LQR(mechanism, bodyid::Int64, Q::AbstractMatrix{T}, R::AbstractMatrix{T}, horizon;
        xd::AbstractVector{T}=zeros(T,3), vd::AbstractVector{T}=zeros(T,3), Fd::AbstractVector{T}=zeros(T,3),
        qd::Quaternion{T}=Quaternion{T}(), ωd::AbstractVector{T}=zeros(T,3), τd::AbstractVector{T}=zeros(T,3)) where T


        Δt = mechanism.Δt
        N = horizon/Δt

        body = getbody(mechanism, bodyid)
        # linearize        
        A, B, G = linearizeMechanism(mechanism, xd, vd, Fd, qd, ωd, τd)

        # calculate K
        if size(G)[1] == 0
            KT, KR = dlqr(A,B,Q,R,N)
        else
            KT, KR = dlqr(A,B,G,Q,R,N)
        end

        if N<Inf
            N=Integer(ceil(horizon/Δt))
        end
        
        new{T, N}(KT, KR, xd, vd, Fd, qd, ωd, τd, bodyid, control_lqr!)
    end
end

function control_lqr!(mechanism, lqr::LQR{T,Inf}, k) where {T}
    Δt = mechanism.Δt
    body = getbody(mechanism, lqr.bodyid)
    state = body.state

    x2 = state.xsol[2]
    v2 = state.vsol[2]
    q2 = state.qsol[2]
    ω2 = state.ωsol[2]

    ΔzpT = [x2-lqr.xd;v2-lqr.vd]
    ΔzpR = [ConstrainedDynamics.VLᵀmat(lqr.qd)*q2;ω2-lqr.ωd]
    Δzp = [ΔzpT;ΔzpR]
    F = -lqr.KT[1] * Δzp
    τ = -lqr.KR[1] * Δzp

    setForce!(body, F=F+lqr.Fd, τ=τ+lqr.τd)
end

function control_lqr!(mechanism, lqr::LQR{T,N}, k) where {T,N}
    Δt = mechanism.Δt
    body = getbody(mechanism, lqr.bodyid)
    state = body.state

    x2 = state.xsol[2]
    v2 = state.vsol[2]
    q2 = state.qsol[2]
    ω2 = state.ωsol[2]

    ΔzpT = [x2-lqr.xd;v2-lqr.vd]
    ΔzpR = [ConstrainedDynamics.VLᵀmat(lqr.qd)*q2;ω2-lqr.ωd]
    Δzp = [ΔzpT;ΔzpR]
    if k<N
        F = -lqr.KT[k] * Δzp
        τ = -lqr.KR[k] * Δzp
    else
        F = @SVector zeros(T,3)
        τ = @SVector zeros(T,3)
    end

    setForce!(body, F=F+lqr.Fd, τ=τ+lqr.τd)
end

function dlqr(A,B,Q,R,N)
    if N==Inf
        P = dare(A,B,Q,R)
        K = (R+B'*P*B)\B'*P*A
        return [K[1:3,:]], [K[4:6,:]]
    else
        N = Integer(ceil(N))
        KT = [zeros(3,12) for i=1:N-1]
        KR = [zeros(3,12) for i=1:N-1]
        Pk = Q
        for k=N-1:-1:1
            Kk = (R+B'*Pk*B)\B'*Pk*A
            KT[k] = Kk[1:3,:]
            KR[k] = Kk[4:6,:]

            Pk = A'*Pk*A - A'*Pk*B/Kk + Q
        end
        return KT, KR
    end
end

function dlqr(A,B,G,Q,R,N)
    if N==Inf
        @error "not supported"
    else
        N = Integer(ceil(N))
        n = size(G)[1]
        m = size(B)[2]
        r = minimum([n;m])

        ZGU = zeros(n, m)
        KT = [zeros(3,12) for i=1:N-1]
        KR = [zeros(3,12) for i=1:N-1]

        Pk = Q
        Hk = G

        for k=N-1:-1:1
            Mxxk = Q + A'*Pk*A
            Muuk = R + B'*Pk*B
            Muxk = B'*Pk*A

            Cxk = [G;Hk*A]
            Cuk = [ZGU;Hk*B]

            VSVDk = svd(Cuk, full=true)
            Vck = VSVDk.V[:,1:r]
            Vuck = VSVDk.V[:,r+1:m]
            

            Kk = Vck*pinv(Cuk*Vck)*Cxk + Vuck/(Vuck'*Muuk*Vuck)*Vuck'*Muxk
            KT[k] = Kk[1:3,:]
            KR[k] = Kk[4:6,:]


            Pk = Mxxk - 2*Muxk'*Kk + Kk'*Muuk*Kk
            Hk = (I - Cuk*Vck*pinv(Cuk*Vck))*Cxk
        end
        return KT, KR
    end
end

