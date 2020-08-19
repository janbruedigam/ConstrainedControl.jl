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
link4 = Body(box1)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))
joint2 = EqualityConstraint(Revolute(origin, link2, joint_axis; p2=p2))
joint3 = EqualityConstraint(Revolute(link1, link3, joint_axis; p1=-p2,p2=p2))
joint4 = EqualityConstraint(Revolute(link4, link2, joint_axis; p1=p2,p2=-p2),CylindricalFree(link4, link3, joint_axis; p1=-p2,p2=-p2))

links = [link1;link2;link3;link4]
constraints = [joint1;joint2;joint3;joint4]
shapes = [box1]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=0)
setPosition!(origin,link1,p2 = p2,Δq = UnitQuaternion(RotX(0-pi/4+0.0)))
setPosition!(origin,link2,p2 = p2,Δq = UnitQuaternion(RotX(pi/4-0.0)))
setPosition!(link1,link3,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(pi/2-0)))
setPosition!(link2,link4,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(-pi/2+0)))

xd = [links[i].state.xc for i=1:4]
qd=[links[i].state.qc for i=1:4]

ang = 1.0
setPosition!(origin,link1,p2 = p2,Δq = UnitQuaternion(RotX(0-pi/4+ang)))
setPosition!(origin,link2,p2 = p2,Δq = UnitQuaternion(RotX(pi/4-ang)))
setPosition!(link1,link3,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(pi/2-2*ang)))
setPosition!(link2,link4,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(-pi/2+2*ang)))

Q = [diagm(ones(12))*0.0 for i=1:4]
Q[1][7,7]=1.0
# Q[1][10,10]=1.0
Q[2][7,7]=1.0
# Q[2][10,10]=1.0
# Q[3][7,7]=1.0
# Q[3][10,10]=1.0
# Q[4][7,7]=1.0
# Q[4][10,10]=1.0
R = [ones(1,1) for i=1:2]

lqr = LQR(mech, getid.(links), [getid(constraints[1]);getid(constraints[2])], Q, R, Inf, xd=xd, qd=qd)
A, Bu, Bλ, G = linearsystem(mech, xd, [zeros(3) for i=1:4], qd, [zeros(3) for i=1:4], [[0] for i=1:2], getid.(links), [getid(constraints[1]);getid(constraints[2])])

steps = Base.OneTo(1000)
storage = Storage{Float64}(steps,4)
storage = simulate!(mech,storage,lqr,record = true)
# storage = simulate!(mech,storage,record = true)
visualize(mech,storage,shapes)