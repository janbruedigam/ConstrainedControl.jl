using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra
using Rotations
using Plots

# Parameters
joint_axis = [1.0;0.0;0.0]

l1 = 1.0
l2 = 2.0
m1 = 1.0
m2 = 1.0
width, depth = 0.1, 0.1
box1 = Box(width, depth, l1, m1)
box2 = Box(width, depth, l2, m2)

p2a = [0.0;0.0;l1/2] # joint connection point
p2b = [0.0;0.0;l2/2] # joint connection point

# Desired orientation
ϕ1 = pi
ϕ2 = pi

# Links
origin = Origin{Float64}()
link1 = Body(box1)
link2 = Body(box2)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2a))
joint2 = EqualityConstraint(Revolute(link1, link2, joint_axis; p1=-p2a, p2=p2b))

links = [link1;link2]
constraints = [joint1;joint2]
shapes = [box1;box2]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81, Δt = 0.001)

xd = [[[0;0;l1/2]];[[0;0;l1+l2/2]]]
qd=[[UnitQuaternion(RotX(ϕ1))];[UnitQuaternion(RotX(ϕ2))]]

Q = [diagm(ones(12))*0.0 for i=1:2]
Q[1][7,7]=4.0
Q[1][10,10]=1.0
Q[2][7,7]=4.0
Q[2][10,10]=1.0
R = [ones(1,1)]

lqr = LQR(mech, getid.(links), [getid(constraints[2])], Q, R, Inf, xd=xd, qd=qd)

function contr!(mechanism::Mechanism{T,N,Nb}, k) where {T,N,Nb}
    Δz = zeros(T,Nb*12)

    mincoords = (minimalCoordinates(mechanism)[3][1],minimalCoordinates(mechanism)[4][1])
    if mincoords[1]<0
        mincoords = (mincoords[1]+2pi,mincoords[2])
    end
    vels = [mechanism.bodies[1].state.ωc[1];mechanism.bodies[2].state.ωc[1]-mechanism.bodies[1].state.ωc[1]]
    
    u = ([242.52*(mincoords[1]-pi)] + [96.33*mincoords[2]] + [104.59*vels[1]] + [49.05*vels[2]])
    setForce!(mechanism, geteqconstraint(mechanism, 4), u)

    return
end

steps = Base.OneTo(25000)

s1 = 101
s2 = 101
# storages = [Storage{Float64}(steps,2) for i=1:s1*s2]
storage = Storage{Float64}(steps,2)
successful = [false for i=1:s1*s2]
angles = [[0.0;0.0] for i=1:s1*s2]
counter = 1

for i=1:s1
    for j=1:s2
        global counter
        setPosition!(origin,link1,p2 = p2a,Δq = UnitQuaternion(RotX(ϕ1-pi+2pi*(i-1)/(s1-1))))
        setPosition!(link1,link2,p1=-p2a,p2 = p2b,Δq = UnitQuaternion(RotX(-pi+2pi*(j-1)/(s2-1))))
        setVelocity!(link1)
        setVelocity!(link2)
        angles[counter] = [ϕ1-pi+2pi*(i-1)/(s1-1);-pi+2pi*(j-1)/(s2-1)]
        display(string(i)*string(j))

        joint1.λsol[1] *= 0
        joint1.λsol[2] *= 0
        joint2.λsol[1] *= 0
        joint2.λsol[2] *= 0
        try
            # simulate!(mech,storages[counter],lqr,record = true)
            # simulate!(mech,storages[counter],contr!,record = true)
            simulate!(mech,steps,storage,lqr,record = true)
            # simulate!(mech,steps,storage,contr!)
        catch
            display("Failed for "*string(i)*" and "*string(j))
        end
        if !any(getindex.(storage.ω[1],1).>100*pi) && !any(getindex.(storage.ω[2],1).>100*pi) && norm([link1.state.xc;link2.state.xc;link1.state.vc;link2.state.vc] - [xd[1];xd[2];zeros(6)])<1e-1
            successful[counter] = true
        end
        counter += 1
    end
end

posres = findall(x->x==true, successful)
A = [getindex.(angles[posres],1) getindex.(angles[posres],2)]
scatter(A[:,1],A[:,2])

# visualize(mech,storages[43],shapes)