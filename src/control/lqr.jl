# TODO
# Only for constant reference trajectories
# Only for 1dof joints
mutable struct LQR{T,N,NK} <: Controller
    K::Vector{Vector{SMatrix{1,NK,T,NK}}} # for each time step and each eqc

    xd::Vector{SVector{3,Float64}}
    vd::Vector{SVector{3,Float64}}
    qd::Vector{UnitQuaternion{T}}
    ωd::Vector{SVector{3,Float64}}

    eqcids::Vector{Integer}
    Fτd::Vector{SVector{1,Float64}}

    control!::Function


    function LQR(mechanism::Mechanism{T,Nn,Nb}, bodyids::AbstractVector{<:Integer}, eqcids::AbstractVector{<:Integer},
            Q::Vector{<:AbstractMatrix{T}}, R::Vector{<:AbstractMatrix{T}}, horizon;
            xd::Vector{<:AbstractVector{T}} = [SA{T}[0; 0; 0] for i=1:Nb], 
            vd::Vector{<:AbstractVector{T}} = [SA{T}[0; 0; 0] for i=1:Nb],
            qd::Vector{UnitQuaternion{T}} = [one(UnitQuaternion{T}) for i=1:Nb], 
            ωd::Vector{<:AbstractVector{T}} = [SA{T}[0; 0; 0] for i=1:Nb],
            Fτd::Vector{<:AbstractVector{T}} = [SA{T}[0] for i=1:length(eqcids)]
        ) where {T, Nn, Nb}

        @assert length(bodyids) == length(Q) == length(xd) == length(vd) == length(qd) == length(ωd) == Nb "Missmatched length for bodies"
        @assert length(eqcids) == length(R) == length(Fτd) "Missmatched length for constraints"

        Δt = mechanism.Δt
        
        N = horizon/Δt
        if N<Inf
            N = Integer(ceil(horizon/Δt))
        end

        # linearize        
        A, Bu, Bλ, G = linearsystem(mechanism, xd, vd, qd, ωd, Fτd, bodyids, eqcids)

        Q = cat(Q...,dims=(1,2))
        R = cat(R...,dims=(1,2))

        # calculate K
        if size(G)[1] == 0
            @assert size(Bλ)[2] ==0
            Ku = dlqr(A, Bu, Q, R, N)
        else
            Ku = dlqr(A, Bu, Bλ, G, Q, R, N)
        end
        
        new{T, N, size(Ku[1][1])[2]}(Ku, xd, vd, qd, ωd, eqcids, Fτd, control_lqr!)
    end
end

function control_lqr!(mechanism::Mechanism{T,Nn,Nb}, lqr::LQR{T,N}, k) where {T,Nn,Nb,N}
    Δz = zeros(T,Nb*12)
    for (id,body) in pairs(mechanism.bodies)
        colx = (id-1)*12+1:(id-1)*12+3
        colv = (id-1)*12+4:(id-1)*12+6
        colq = (id-1)*12+7:(id-1)*12+9
        colω = (id-1)*12+10:(id-1)*12+12

        state = body.state
        Δz[colx] = state.xsol[2]-lqr.xd[id]
        Δz[colv] = state.vsol[2]-lqr.vd[id]
        Δz[colq] = ConstrainedDynamics.VLᵀmat(lqr.qd[id]) * Rotations.params(state.qsol[2])
        Δz[colω] = state.ωsol[2]-lqr.ωd[id]
    end

    if k<N
        for (i,id) in enumerate(lqr.eqcids)
            u = lqr.Fτd[i] - lqr.K[k][i]*Δz
            setForce!(mechanism, geteqconstraint(mechanism, id), u)
        end
    end

    return
end

function control_lqr!(mechanism::Mechanism{T,Nn,Nb}, lqr::LQR{T,Inf}, k) where {T,Nn,Nb}
    Δz = zeros(T,Nb*12)
    for (id,body) in pairs(mechanism.bodies)
        colx = (id-1)*12+1:(id-1)*12+3
        colv = (id-1)*12+4:(id-1)*12+6
        colq = (id-1)*12+7:(id-1)*12+9
        colω = (id-1)*12+10:(id-1)*12+12

        state = body.state
        Δz[colx] = state.xsol[2]-lqr.xd[id]
        Δz[colv] = state.vsol[2]-lqr.vd[id]
        Δz[colq] = ConstrainedDynamics.VLᵀmat(lqr.qd[id]) * Rotations.params(state.qsol[2])
        Δz[colω] = state.ωsol[2]-lqr.ωd[id]
    end

    for (i,id) in enumerate(lqr.eqcids)
        u = lqr.Fτd[i] - lqr.K[1][i]*Δz
        setForce!(mechanism, geteqconstraint(mechanism, id), u)
    end

    return
end

function dlqr(A,B,Q,R,N)
    if N==Inf
        P = dare(A,B,Q,R)
        K = (R+B'*P*B)\B'*P*A
        return [[K[i:i,:] for i=1:size(K)[1]]]
    else
        return dlqr(A,B,zeros(size(A)[1],0),zeros(0,size(A)[1]),Q,R,N)
    end
end

function dlqr(A,Bu,Bλ,G,Q,R,N)
    infflag = false
    if N == Inf
        infflag = true
        N = 1000
    end

    mx = size(A)[2]
    mu = size(Bu)[2]
    mλ = size(Bλ)[2]
    Ku = [[zeros(1,size(Q)[1]) for j=1:mu] for i=1:N-1]
    Kλ = [[zeros(1,size(Q)[1]) for j=1:mλ] for i=1:N-1]
    Pk = Q

    k = 0
    for outer k=N-1:-1:1
        M11 = R + Bu'*Pk*Bu
        M12 = Bu'*Pk*Bλ
        M21 = G*Bu
        M22 = G*Bλ

        M = [M11 M12;M21 M22]
        b = [Bu'*Pk;G]*A

        Kk = M\b

        for i=1:mu
            Ku[k][i] = Kk[i:i,:]
        end

        Kuk = Kk[1:mu,:]
        Kλk = Kk[mu+1:mu+mλ,:]

        Abar = A-Bu*Kuk-Bλ*Kλk
        Pkp1 = Q + Kuk'*R*Kuk + Abar'*Pk*Abar

        if infflag && norm(Pk-Pkp1) < 1e-5
            break
        end

        Pk = Pkp1
    end

    if infflag
        if k==1
            @info "Riccati recursion did not converge."
        else
            Ku = [Ku[k]]
        end
    end

    return Ku
end
