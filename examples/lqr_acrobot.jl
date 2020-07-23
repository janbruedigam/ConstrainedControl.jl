using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra
using Rotations


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box1 = Box(width, depth, length1, length1)
box2 = Box(width, depth, length1*2, length1)

p2a = [0.0;0.0;length1/2] # joint connection point
p2b = [0.0;0.0;length1] # joint connection point
ϕ = 0.

# Initial orientation

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


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81)
setPosition!(origin,link1,p2 = p2a,Δq = UnitQuaternion(RotX(ϕ+pi+0.1)))
# setPosition!(origin,link1,p2 = p2a,Δq = UnitQuaternion(RotX(ϕ-0.)))
setPosition!(link1,link2,p1=-p2a,p2 = p2b,Δq = UnitQuaternion(RotX(ϕ-0.1)))

xd = [[[0;0;0.5]];[[0;0;2.0]]]
# xd = [[[0;0;0.5]];[[0;0;0.]]]
# xd = [[-[0;0;0.5]];[-[0;0;2.0]]]
qd=[[UnitQuaternion(RotX(ϕ+pi))];[UnitQuaternion(RotX(ϕ+pi))]]
# qd=[[UnitQuaternion(RotX(ϕ+pi))];[UnitQuaternion(RotX(ϕ))]]
# qd=[[UnitQuaternion(RotX(ϕ))];[UnitQuaternion(RotX(ϕ))]]

Q = [diagm(ones(12))*0.0 for i=1:2]
Q[1][7,7]=1.0
Q[1][10,10]=1.0
Q[2][7,7]=1.0
Q[2][10,10]=1.0
R = [ones(1,1)]

A, Bu, Bλ, G = linearsystem(mech, xd, [zeros(3) for i=1:2], qd, [zeros(3) for i=1:2], [[0.]], getid.(links), [getid(constraints[2])])
lqr = LQR(mech, getid.(links), [getid(constraints[2])], Q, R, 10., xd=xd, qd=qd)

# function contr!(mechanism::Mechanism{T,N,Nb}, k) where {T,N,Nb}
#     Δz = zeros(T,Nb*12)

#     mincoords = (minimalCoordinates(mechanism)[3][1],minimalCoordinates(mechanism)[4][1])
#     if mincoords[1]<0
#         mincoords = (mincoords[1]+2pi,mincoords[2])
#     end
#     vels = [mechanism.bodies[1].state.ωc[1];mechanism.bodies[2].state.ωc[1]-mechanism.bodies[1].state.ωc[1]]
#     # display(mincoords)

    
#     u = ([242.52*(mincoords[1]-pi)] + [96.33*mincoords[2]] + [104.59*vels[1]] + [49.05*vels[2]])
#     setForce!(mechanism, geteqconstraint(mechanism, 4), u)

#     return
# end

steps = Base.OneTo(1000)
storage = Storage{Float64}(steps,length(mech.bodies))

simulate!(mech,storage,lqr,record = true)
# simulate!(mech,storage,record = true)
# simulate!(mech,storage,contr!,record = true)
visualize(mech,storage,shapes)