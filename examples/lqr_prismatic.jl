using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, width, length1)

# Links
origin = Origin{Float64}()
link1 = Body(box)

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(Prismatic(origin, link1, joint_axis))

links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=0.)
setPosition!(origin,link1,Δx = [1.0;0;0])

Q = [diagm(ones(12))]
R = [ones(1,1)]

lqr = LQR(mech, getid.(links), getid.(constraints), Q, R, 10.)


storage = simulate!(mech,10.,lqr,record = true)
visualize(mech,storage,shapes)