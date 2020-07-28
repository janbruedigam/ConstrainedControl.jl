using ConstrainedDynamics
using ConstrainedControl
using LinearAlgebra
using Rotations


# Parameters
joint_axis = [1.0;0.0;0.0]
p2a = [0.0;0.0;length1/2] # joint connection point
p2b = [0.0;0.0;length1] # joint connection point

length1 = 1.0
width, depth = 0.1, 0.1
box1 = Box(width, depth, length1, length1)
box2 = Box(width, depth, length1*2, length1)

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


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81)
setPosition!(origin,link1,p2 = p2a,Δq = UnitQuaternion(RotX(ϕ1-0.1)))
setPosition!(link1,link2,p1=-p2a,p2 = p2b,Δq = UnitQuaternion(RotX(0.1)))

xd = [[[0;0;0.5]];[[0;0;2.0]]]
qd=[[UnitQuaternion(RotX(ϕ1))];[UnitQuaternion(RotX(ϕ2))]]

Q = [diagm(ones(12))*0.0 for i=1:2]
Q[1][7,7]=1.0
Q[1][10,10]=1.0
Q[2][7,7]=1.0
Q[2][10,10]=1.0
R = [ones(1,1)]

lqr = LQR(mech, getid.(links), [getid(constraints[2])], Q, R, 10., xd=xd, qd=qd)
@test true
