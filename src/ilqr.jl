mutable struct iLQR{T,N} <: Controller
    storage0::Storage{T,N}
    storage1::Storage{T,N}
    steps::Base.OneTo{Int64}
    F0::Vector{SVector{3,T}}
    τ0::Vector{SVector{3,T}}
    F1::Vector{SVector{3,T}}
    τ1::Vector{SVector{3,T}}
    α::T
    K::Vector{SMatrix{6,12,T,72}}
    d::Vector{SVector{6,T}}

    bodyid::Int64

    control!::Function

    function iLQR(mechanism, steps, bodyid::Int64, F0::Vector{<:AbstractVector{T}}, τ0::Vector{<:AbstractVector{T}}) where T

        storage0 = Storage{T}(steps,1)
        storage1 = Storage{T}(steps,1)
        F1 = F0
        τ1 = τ0
        α = 1.0
        K = [(@SMatrix zeros(T,6,12)) for i=steps]
        d = [(@SVector zeros(T,6)) for i=steps]
        
        new{T,length(steps)}(storage0, storage1, steps, F0, τ0, F1, τ1, α, K, d, bodyid, control_ilqr!)
    end
end

function control_ilqr!(mechanism, ilqr::iLQR, k)
    Δt = mechanism.Δt
    storage0 = ilqr.storage0
    F0 = ilqr.F0
    τ0 = ilqr.τ0
    K = ilqr.K
    d = ilqr.d
    if k==last(ilqr.steps)
        k-=1
    end
    body = getbody(mechanism,ilqr.bodyid)
    xbar = body.x[2]
    vbar = ConstrainedDynamics.getvnew(body)
    qbar = body.q[2]
    ωbar = ConstrainedDynamics.getωnew(body)

    Δx = xbar - storage0.x[1][k]
    Δv = vbar - storage0.v[1][k]
    Δϵ = ConstrainedDynamics.VLᵀmat(storage0.q[1][k])*qbar
    Δω = ωbar - storage0.ω[1][k]

    Δstate = [Δx;Δv;Δϵ;Δω]

    unew = [F0[k];τ0[k]] - K[k]*Δstate - ilqr.α*d[k]
    ilqr.F1[k] = unew[1:3]
    ilqr.τ1[k] = unew[4:6]

    setForce!(mechanism, body, F=ilqr.F1[k], τ = ilqr.τ1[k])
end

function ilqr(mechanism, bodyid::Int64, steps, Q::AbstractMatrix{T}, R::AbstractMatrix{T}, F0::Vector{<:AbstractVector{T}}, τ0::Vector{<:AbstractVector{T}};
    xd::AbstractVector{T}=zeros(T,3), vd::AbstractVector{T}=zeros(T,3), qd::Quaternion{T}=Quaternion{T}(), ωd::AbstractVector{T}=zeros(T,3)) where T

    Δt = mechanism.Δt
    ilqrobj = iLQR(mechanism, steps, bodyid, F0, τ0)

    # Initialize, set states and J
    ilqrobj.storage0 = simulate!(mechanism,steps,ilqrobj.storage0,ilqrobj,record = true)
    J0 = evalJ(ilqrobj.storage0,steps,Δt,Q,R,xd, vd, qd, ωd,ilqrobj.F0,ilqrobj.τ0)
    J1 = J0
    n = 0
    for outer n=1:100
        K, d, ΔV = backwardpass(mechanism, ilqrobj.storage0,steps, getbody(mechanism,bodyid), Q, R, xd, vd, qd, ωd, ilqrobj.F0, ilqrobj.τ0) # get controller
        ilqrobj.K = K
        ilqrobj.d = d
        # ilqr.ΔV = ΔV
        # Forwardpass
        for i=1:10 #linesearch
            ilqrobj.α = 1.0
            ilqrobj.storage1 = simulate!(mechanism,steps,ilqrobj.storage1,ilqrobj,record = true)
            # F1, τ1 = newF(mechanism, storage0, storage1, steps, F0,τ0, K,d,α)
            J1 = evalJ(ilqrobj.storage1,steps,Δt,Q,R,xd, vd, qd, ωd,ilqrobj.F1,ilqrobj.τ1)
            if J1<=J0
                break
            else
                ilqrobj.α = ilqrobj.α/2
            end
        end
        
        # End Forwardpass
        if abs(J1-J0)<1e-5
            break
        else
            ilqrobj.F0 = ilqrobj.F1
            ilqrobj.τ0 = ilqrobj.τ1
            J0 = J1
            ilqrobj.storage0 = ilqrobj.storage1
        end
    end

    return J0, n, ilqrobj
end

function backwardpass(mechanism, storage, steps, body, Q, R, xd, vd, qd, ωd, F, τ)
    Δt = mechanism.Δt
    p1 = zeros(12)
    P1 = Q

    K = [zeros(6,12) for i=steps]
    d = [zeros(6) for i=steps]
    ΔV = zeros(steps)

    for k=Iterators.drop(reverse(steps),1)
        x = storage.x[1][k]
        v = storage.v[1][k]
        q = storage.q[1][k]
        ω = storage.ω[1][k]

        state = [x;v;Vmat(q);ω]
        u = [F[k];τ[k]]

        A, B = ConstrainedDynamics.∂zp1∂z(mechanism, body, x, v, F[k], q, ω, τ[k], Δt)

        Qxx = Q + A'*P1*A
        Quu = R + B'*P1*B
        Qux = B'*P1*A
        Qxu = A'*P1*B
        Qx = Q*state + A'*p1
        Qu = R*u + B'*p1

        K[k] = Quu\Qux
        d[k] = Quu\Qu
        ΔV[k] = d[k]'*Qu + 1/2*d[k]'*Quu*d[k]

        P0 = Qxx + K[k]'*Quu*K[k] + K[k]'*Qux + Qxu*K[k]
        p0 = Qx + K[k]'*Quu*d[k] + K[k]'*Qu + Qxu*d[k]

        P1 = P0
        p1 = p0
    end


    return K, d, ΔV
end

function evalJ(storage,steps,Δt,Q,R,xd, vd, qd, ωd,F,τ)
    J = 0
    for k=Iterators.take(steps,length(steps)-1)
        Δx = storage.x[1][k]-xd
        Δv = (storage.x[1][k+1]-storage.x[1][k])/Δt - vd
        Δϵ = ConstrainedDynamics.VLᵀmat(qd)*storage.q[1][k]
        Δω = 2/Δt * ConstrainedDynamics.VLᵀmat(storage.q[1][k])*storage.q[1][k+1] - ωd

        Δstate = [Δx;Δv;Δϵ;Δω]
        u = [F[k];τ[k]]

        J += Δstate'*Q*Δstate + u'*R*u
    end
    return J
end

function newF(mechanism, storage0, storage1, steps, F0,τ0, K,d,α)
    Δt = mechanism.Δt
    for k=Iterators.take(steps,length(steps)-1)
        xbar = storage1.x[1][k]
        vbar = storage1.v[1][k]
        qbar = storage1.q[1][k]
        ωbar = storage1.ω[1][k]

        Δx = xbar - storage0.x[1][k]
        Δv = vbar - storage0.v[1][k]
        Δϵ = ConstrainedDynamics.VLᵀmat(storage0.q[1][k])*qbar
        Δω = ωbar - storage0.ω[1][k]

        Δstate = [Δx;Δv;Δϵ;Δω]

        unew = [F0[k];τ0[k]] - K[k]*Δstate - α*d[k]
        F0[k] = unew[1:3]
        τ0[k] = unew[4:6]
    end
    return F0, τ0
end