using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, width, length1)

p2 = [0.0;0.0;0.0] # joint connection point

# Initial orientation

# Links
origin = Origin{Float64}()
link1 = Body(box)

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(Prismatic(origin, link1, joint_axis; p2=p2))

links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=0.)
setPosition!(origin,link1,p2 = p2,Δx = [1.0;0;0])

xd=[zeros(3)]
vd=[zeros(3)]
qd=[UnitQuaternion(RotX(pi*0))]
ωd=[zeros(3)]

Q = [diagm(ones(12))]
R = [ones(1,1)]

A, B, G = linearsystem(mech, xd, vd, qd, ωd, [[0.]], getid.(links), getid.(constraints))

lqr = LQR(mech, getid.(links), getid.(constraints), Q, R, 10., xd=xd, qd=qd)


storage = simulate!(mech,10.,lqr,record = true)
visualize(mech,storage,shapes)