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

p2 = [0.0;0.0;length1/2] # joint connection point

# Desired orientation
ϕ1 = 0
ϕ2 = 0

# Links
origin = Origin{Float64}()
link1 = Body(box1)
link2 = Body(box1)
link3 = Body(box1)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))
joint2 = EqualityConstraint(Revolute(link1, link2, joint_axis; p1=-p2, p2=p2))
joint3 = EqualityConstraint(Revolute(link3, link2, joint_axis; p1=p2,p2=-p2),CylindricalFree(link3, link1, joint_axis; p1=-p2,p2=p2))
# joint3 = EqualityConstraint(Revolute(link3, link2, joint_axis; p1=p2,p2=-p2))

links = [link1;link2;link3]
constraints = [joint1;joint2;joint3]
shapes = [box1]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81)
setPosition!(origin,link1,p2 = p2,Δq = UnitQuaternion(RotX(-pi/6+pi)))
setPosition!(link1,link2,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(2pi/3)))
setPosition!(link2,link3,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(2pi/3)))

xd = [links[i].state.xc for i=1:3]
qd=[links[i].state.qc for i=1:3]

setPosition!(origin,link1,p2 = p2,Δq = UnitQuaternion(RotX(-pi/6+pi+0.4)))
setPosition!(link1,link2,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(2pi/3)))
setPosition!(link2,link3,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(2pi/3)))

Q = [diagm(ones(12))*0.0 for i=1:3]
Q[1][7,7]=1.0
Q[1][10,10]=1.0
# Q[2][7,7]=1.0
# Q[2][10,10]=1.0
# Q[3][7,7]=1.0
# Q[3][10,10]=1.0
R = [ones(1,1) for i=1:1]

lqr = LQR(mech, getid.(links), [getid(constraints[1])], Q, R, 10., xd=xd, qd=qd)
A, Bu, Bλ, G = linearsystem(mech, xd, [zeros(3) for i=1:3], qd, [zeros(3) for i=1:3], [[0] for i=1:1], getid.(links), [getid(constraints[1])])

steps = Base.OneTo(1000)
storage = Storage{Float64}(steps,3)
storage = simulate!(mech,storage,lqr,record = true)
# storage = simulate!(mech,storage,record = true)
visualize(mech,storage,shapes)