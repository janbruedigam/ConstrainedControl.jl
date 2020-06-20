using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)

p2 = [0.0;0.0;length1/2] # joint connection point
ϕ = 0.

# Initial orientation

# Links
origin = Origin{Float64}()
link1 = Body(box)
link2 = Body(box)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))
joint2 = EqualityConstraint(Revolute(link1, link2, joint_axis; p1=-p2, p2=p2))

links = [link1;link2]
constraints = [joint1;joint2]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=0.)
setPosition!(origin,link1,p2 = p2,Δq = UnitQuaternion(RotX(ϕ+pi-0.2)))
setPosition!(link1,link2,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(ϕ+0.3)))

xd=[[p2];[p2+p2+p2]]
qd=[[UnitQuaternion(RotX(ϕ+pi))];[UnitQuaternion(RotX(ϕ+pi))]]

Q = [diagm(ones(12))*0.0 for i=1:2]
Q[1][7,7]=100.0
Q[1][10,10]=10.0
Q[2][7,7]=100.0
Q[2][10,10]=10.0
R = [ones(1,1)]

lqr = LQR(mech, getid.(links), [getid(constraints[2])], Q, R, 10.,  qd=qd)


storage = simulate!(mech,10.,lqr,record = true)
visualize(mech,storage,shapes)