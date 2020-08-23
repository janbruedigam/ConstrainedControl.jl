# TODO
# Only for 1dof joints
mutable struct TrackingLQR{T,N,NK} <: Controller
    K::Vector{Vector{SMatrix{1,NK,T,NK}}} # for each time step and each eqc

    xd::Vector{Vector{SVector{3,Float64}}} # for each time step and each eqc
    vd::Vector{Vector{SVector{3,Float64}}}# for each time step and each eqc
    qd::Vector{Vector{UnitQuaternion{T}}} # for each time step and each eqc
    ωd::Vector{Vector{SVector{3,Float64}}} # for each time step and each eqc

    eqcids::Vector{Integer}
    Fτd::Vector{Vector{SVector{1,Float64}}} # for each time step and each eqc

    control!::Function


    function TrackingLQR(mechanism::Mechanism{T,Nn,Nb}, storage::Storage{T,N}, Fτ::Vector{<:Vector{<:AbstractVector{T}}}, eqcids::AbstractVector{<:Integer},
            Q::Vector{<:AbstractMatrix{T}}, R::Vector{<:AbstractMatrix{T}}
        ) where {T, Nn, Nb, N}

        Q = cat(Q...,dims=(1,2))
        R = cat(R...,dims=(1,2))

        xd = [[SA{T}[0; 0; 0] for i=1:Nb] for j=1:N]
        vd = [[SA{T}[0; 0; 0] for i=1:Nb] for j=1:N]
        qd = [[one(UnitQuaternion{T}) for i=1:Nb] for j=1:N]        
        ωd = [[SA{T}[0; 0; 0] for i=1:Nb] for j=1:N]

        for k = 1:N
            for i = 1:Nb
                xd[k][i] = storage.x[i][k]
                vd[k][i] = storage.v[i][k]
                qd[k][i] = storage.q[i][k]
                ωd[k][i] = storage.ω[i][k]
            end
        end

        # calculate K
        Ku = dlqr(mechanism, xd, vd, qd, ωd, Fτ, eqcids,Q,R,N)

        new{T, N, size(Ku[1][1])[2]}(Ku, xd, vd, qd, ωd, eqcids, Fτ, control_trackinglqr!)
    end
end

function control_trackinglqr!(mechanism::Mechanism{T,Nn,Nb}, lqr::TrackingLQR{T,N}, k) where {T,Nn,Nb,N}
    Δz = zeros(T,Nb*12)
    qvm = QuatVecMap()
    for (id,body) in pairs(mechanism.bodies)
        colx = (id-1)*12+1:(id-1)*12+3
        colv = (id-1)*12+4:(id-1)*12+6
        colq = (id-1)*12+7:(id-1)*12+9
        colω = (id-1)*12+10:(id-1)*12+12

        state = body.state
        Δz[colx] = state.xsol[2]-lqr.xd[k][id]
        Δz[colv] = state.vsol[2]-lqr.vd[k][id]
        Δz[colq] = rotation_error(state.qsol[2],lqr.qd[k][id],qvm)
        Δz[colω] = state.ωsol[2]-lqr.ωd[k][id]
    end

    if k<N
        for (i,id) in enumerate(lqr.eqcids)
            u = lqr.Fτd[k][i] - lqr.K[k][i]*Δz
            setForce!(mechanism, geteqconstraint(mechanism, id), u)
        end
    end

    return
end

function dlqr(mechanism::Mechanism{T,Nn,Nb,Ne}, xd, vd, qd, ωd, Fτd, eqcids,Q,R,N) where {T,Nn,Nb,Ne}
    bodyids = getid.(mechanism.bodies)
    mx = Nb*6*2
    mu = size(eqcids)[1]
    mλ = 0
    for eqc in mechanism.eqconstraints
        ConstrainedDynamics.isinactive(eqc) && continue
        mλ += ConstrainedDynamics.length(eqc)
    end
    Ku = [[zeros(1,size(Q)[1]) for j=1:mu] for i=1:(N-1)]
    Kλ = [[zeros(1,size(Q)[1]) for j=1:mλ] for i=1:(N-1)]
    Pk = Q

    k = 0
    for outer k=N-1:-1:1
        A, Bu, Bλ, G = linearsystem(mechanism, xd[k], vd[k], qd[k], ωd[k], Fτd[k], bodyids, eqcids)
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

        if norm(Pk-Pkp1) < 1e-5
            break
        end

        Pk = Pkp1
    end

    for k2=k-1:-1:1
        Ku[k2] = Ku[k2+1]
    end

    return Ku
end
