using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1

p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = 0
q1 = Quaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))

links = [link1]
constraints = [joint_between_origin_and_link1]


mech = Mechanism(origin, links, constraints)
setPosition!(origin,link1,p2 = p2,Δq = q1)

pid = PID(mech, joint_between_origin_and_link1.id, pi/2, P = 10., I = 10., D = 5.)


storage = simulate!(mech,10.,pid,record = true)
visualize(mech,storage)