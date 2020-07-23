using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra


# Parameters
joint_axis = [0.0;1.0;0.0]

box1 = Box(0.1, 1.0, 0.1, 1.0)

p2 = [0.0;1/2;0] # joint connection point

# Initial orientation

# Links
origin = Origin{Float64}()
link1 = Body(box1)
link2 = Body(box1)

# Constraints
joint1 = EqualityConstraint(Prismatic(origin, link1, joint_axis; p2=-p2))
joint2 = EqualityConstraint(Prismatic(link1, link2, joint_axis; p1=p2, p2=-p2))

links = [link1;link2]
constraints = [joint1;joint2]
# links = [link1]
# constraints = [joint1]
shapes = [box1]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-0.0)
setPosition!(origin, link1, p2 = -p2)
setPosition!(link1, link2, p1=p2,p2 = -p2)

# xd = [zeros(3) for i=1:2]
# xd = [[zeros(3)];[[0.0;1;0]]]
xd = [[[0.0;1;0]];[[0.0;2;0]]]

Q = [diagm(ones(12))*0.0 for i=1:2]
Q[1][2,2]=1.0
Q[2][2,2]=1.0
R = [ones(1,1)]

A, Bu, Bλ, G = linearsystem(mech, xd, [zeros(3) for i=1:2], [one(UnitQuaternion) for i=1:2], [zeros(3) for i=1:2], [[0.]], getid.(links), [getid(constraints[2])])
lqr = LQR(mech, getid.(links), [getid(constraints[2])], Q, R, 10., xd=xd)

storage = simulate!(mech,10.,lqr,record = true)
visualize(mech,storage,shapes)