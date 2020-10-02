using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra
using Rotations
using Rotations: rotation_error
using Statistics


# Parameters
ex = [1.0;0.0;0.0]
ey = [0.0;1.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
cartshape = Box(0.1, 0.5, 0.1, length1/2)
poleshape = Box(width, depth, length1, length1)

p2 = [0.0;0.0;length1/2] # joint connection point

# Desired orientation
ϕ = 0.

# Links
origin = Origin{Float64}()
cart = Body(cartshape)
pole1 = Body(poleshape)
pole2 = Body(poleshape)
pole3 = Body(poleshape)

# Constraints
joint1 = EqualityConstraint(Prismatic(origin, cart, ey))
joint2 = EqualityConstraint(Revolute(cart, pole1, ex; p2=p2))
joint3 = EqualityConstraint(Revolute(pole1, pole2, ex; p1 = -p2, p2=p2))
joint4 = EqualityConstraint(Revolute(pole2, pole3, ex; p1 = -p2, p2=p2))

links = [cart;pole1;pole2;pole3]
constraints = [joint1;joint2;joint3;joint4]
shapes = [cartshape;poleshape]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81,Δt = 0.01)
setPosition!(origin,cart,Δx = [0;0.0;0])
setPosition!(cart,pole1,p2 = p2,Δq = UnitQuaternion(RotX(0.)))
setPosition!(pole1,pole2,p1 = -p2,p2 = p2,Δq = UnitQuaternion(RotX(0.)))
setPosition!(pole2,pole3,p1 = -p2,p2 = p2,Δq = UnitQuaternion(RotX(0.)))

Q = [diagm(ones(12))*0.0 for i=1:4]
Q[1][2,2] = 10
Q[1][5,5] = 1
Q[2][7,7] = 40
Q[2][10,10] = 1
Q[3][7,7] = 40
Q[3][10,10] = 1
Q[4][7,7] = 40
Q[4][10,10] = 1
R = [ones(1,1)*0.1]

ucost = zeros(1000)

function owncontrol_trackinglqr!(mechanism::Mechanism{T,Nn,Nb}, lqr::TrackingLQR{T,N}, k) where {T,Nn,Nb,N}
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

    v1 = mechanism.bodies[1].state.vc[2]
    ω2 = mechanism.bodies[2].state.ωc[1]
    ω3 = mechanism.bodies[3].state.ωc[1] - ω2
    ω4 = mechanism.bodies[4].state.ωc[1] - ω2 - ω3

    ucart = -sign(v1)*0.1*abs(v1) + randn()*2
    up2 = -sign(ω2)*0.1*abs(ω2)
    up3 = -sign(ω3)*0.1*abs(ω3)
    up4 = -sign(ω4)*0.1*abs(ω4)

    if k<N
        for (i,id) in enumerate(lqr.eqcids)
            u = lqr.Fτd[k][i] - lqr.K[k][i]*Δz + ucart
            ucost[k] = ((lqr.K[k][i]*Δz)'*R*(lqr.K[k][i]*Δz))[1]
            setForce!(mechanism, geteqconstraint(mechanism, id), u)
        end

        setForce!(mechanism, geteqconstraint(mechanism,6), [up2])
        setForce!(mechanism, geteqconstraint(mechanism,7), [up3])
        setForce!(mechanism, geteqconstraint(mechanism,8), [up4])
    end

    return
end

lqr = TrackingLQR(mech, storage0, [[[U[k]]] for k=1:1000], [5], Q, R, controlfunction = owncontrol_trackinglqr!)

function uncontrol!(mechanism, k)
    v1 = mechanism.bodies[1].state.vc[2]
    ω2 = mechanism.bodies[2].state.ωc[1]
    ω3 = mechanism.bodies[3].state.ωc[1] - ω2
    ω4 = mechanism.bodies[4].state.ωc[1] - ω2 - ω3

    ucart = U[k] - sign(v1)*0.1*abs(v1) + randn()*2
    up2 = -sign(ω2)*0.1*abs(ω2)
    up3 = -sign(ω3)*0.1*abs(ω3)
    up4 = -sign(ω4)*0.1*abs(ω4)

    setForce!(mechanism, geteqconstraint(mechanism,5), [ucart])
    setForce!(mechanism, geteqconstraint(mechanism,6), [up2])
    setForce!(mechanism, geteqconstraint(mechanism,7), [up3])
    setForce!(mechanism, geteqconstraint(mechanism,8), [up4])
end

function control!(mechanism, k)
    mincoords = [minimalCoordinates(mechanism)[5][1];minimalCoordinates(mechanism)[6][1];minimalCoordinates(mechanism)[7][1];minimalCoordinates(mechanism)[8][1]]

    v1 = mechanism.bodies[1].state.vc[2]
    ω2 = mechanism.bodies[2].state.ωc[1]
    ω3 = mechanism.bodies[3].state.ωc[1] - ω2
    ω4 = mechanism.bodies[4].state.ωc[1] - ω2 - ω3
    vels = [v1;ω2;ω3;ω4]

    
    c = [mincoords;vels]
    if c[2]<0
        c[2] += 2pi
    end
    if c[3]<0
        c[3] += 2pi
    end
    if c[4]<0
        c[4] += 2pi
    end
    if X[2,k]<0
        X[2,k] += 2pi
    end
    if X[3,k]<0
        X[3,k] += 2pi
    end
    if X[4,k]<0
        X[4,k] += 2pi
    end

    diffval = (c-X[:,k]) 

    if diffval[2] > pi
        diffval[2] = 2pi - c[2] -X[2,k]
    elseif diffval[2] < -pi
        diffval[2] = c[2] - (2pi - X[2,k])
    end
    if diffval[3] > pi
        diffval[3] = 2pi - c[3] -X[3,k]
    elseif diffval[3] < -pi
        diffval[3] = c[3] - (2pi - X[3,k])
    end
    if diffval[4] > pi
        diffval[4] = 2pi - c[4] -X[4,k]
    elseif diffval[4] < -pi
        diffval[4] = c[4] - (2pi - X[4,k])
    end


    ucart = U[k]-K[k,:]'*diffval - sign(vels[1])*0.1*abs(vels[1]) + randn()*2
    ucost[k] = ((K[k,:]'*diffval)'*R[1]*(K[k,:]'*diffval))[1]
    up2 = -sign(vels[2])*0.1*abs(vels[2])
    up3 = -sign(vels[3])*0.1*abs(vels[3])
    up4 = -sign(vels[4])*0.1*abs(vels[4])

    setForce!(mechanism, geteqconstraint(mechanism,5), [ucart])
    setForce!(mechanism, geteqconstraint(mechanism,6), [up2])
    setForce!(mechanism, geteqconstraint(mechanism,7), [up3])
    setForce!(mechanism, geteqconstraint(mechanism,8), [up4])

    return
end




steps = Base.OneTo(1000)
storage = Storage{Float64}(steps,4)

nmax = 1000
cost = [zeros(nmax,1000) for i=1:3]
costsum = [zeros(nmax,1000) for i=1:3]
costmean = [zeros(1000) for i=1:3]
costsummean = [zeros(1000) for i=1:3]
costvar = [zeros(1000) for i=1:3]
costsumvar = [zeros(1000) for i=1:3]
qvm = QuatVecMap()

tmpcm = cart.m
tmpp1m = pole1.m
tmpp2m = pole2.m
tmpp3m = pole3.m
tmpcJ = cart.J
tmpp1J = pole1.J
tmpp2J = pole2.J
tmpp3J = pole3.J

for sel = 1:3
    for n=1:nmax
        cart.m = tmpcm * (1+(rand()*0.2 - 0.1))
        pole1.m = tmpp1m * (1+(rand()*0.2 - 0.1))
        pole2.m = tmpp2m * (1+(rand()*0.2 - 0.1))
        pole3.m = tmpp3m * (1+(rand()*0.2 - 0.1))
        cart.J = tmpcJ * (1+(rand()*0.2 - 0.1))
        pole1.J = tmpp1J * (1+(rand()*0.2 - 0.1))
        pole2.J = tmpp2J * (1+(rand()*0.2 - 0.1))
        pole3.J = tmpp3J * (1+(rand()*0.2 - 0.1))
        
        mod(n,50) == 0 && display(n)
        setPosition!(origin,cart,Δx = [0;0.1-rand()/5;0])
        setPosition!(cart,pole1,p2 = p2,Δq = UnitQuaternion(RotX(0.1-rand()/5)))
        setPosition!(pole1,pole2,p1 = -p2,p2 = p2,Δq = UnitQuaternion(RotX(0.1-rand()/5)))
        setPosition!(pole2,pole3,p1 = -p2,p2 = p2,Δq = UnitQuaternion(RotX(0.1-rand()/5)))
        setVelocity!(cart)
        setVelocity!(pole1)
        setVelocity!(pole2)
        setVelocity!(pole3)
        if sel == 1
            simulate!(mech,storage,lqr,record = true)
        elseif sel == 2
            simulate!(mech,storage,control!,record = true)
        elseif sel == 3
            simulate!(mech,storage,uncontrol!,record = true)
        end

        # visualize(mech,storage,shapes)

        for i = 1:4
            cost[sel][n,1] = (storage.x[i][1][2]-storage0.x[i][1][2])^2*10
            cost[sel][n,1] += (storage.v[i][1][2]-storage0.v[i][1][2])^2*1
            cost[sel][n,1] += Rotations.rotation_error(storage.q[i][1],storage0.q[i][1],qvm)[1]^2*40
            cost[sel][n,1] += (storage.ω[i][1][2]-storage0.ω[i][1][2])^2*1
        end
        cost[sel][n,1] += ucost[1]
        costsum[sel][n,1] = cost[sel][n,1]
        for k=2:1000
            for i = 1:4
                cost[sel][n,k] = (storage.x[i][k][2]-storage0.x[i][k][2])^2*10
                cost[sel][n,k] += (storage.v[i][k][2]-storage0.v[i][k][2])^2*1
                cost[sel][n,k] += Rotations.rotation_error(storage.q[i][k],storage0.q[i][k],qvm)[1]^2*40
                cost[sel][n,k] += (storage.ω[i][k][2]-storage0.ω[i][k][2])^2*1
            end
            cost[sel][n,k] += ucost[k]
            costsum[sel][n,k] = costsum[sel][n,k-1]+cost[sel][n,k]
        end

        costmean[sel] = mean(cost[sel]',dims=2)[:,1]
        costsummean[sel] = mean(costsum[sel]',dims=2)[:,1]
        costvar[sel] = sqrt.(var(cost[sel]',dims=2)[:,1])
        costsumvar[sel] = sqrt.(var(costsum[sel]',dims=2)[:,1])
    end
end


