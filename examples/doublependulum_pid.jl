using ConstrainedDynamics
using ConstrainedControl

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))

p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = 0
q1 = Quaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Body(box)
link2 = Body(box)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, zeros(3), p2, joint_axis))
joint2 = EqualityConstraint(Revolute(link1, link2, -p2, p2, joint_axis))

links = [link1;link2]
constraints = [joint1;joint2]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes)
setPosition!(mech,origin,link1,p2 = p2,Δq = q1)
setPosition!(mech,origin,link2,p1=-p2,p2 = p2,Δq = q1)

pid = PID(mech, getfield.(constraints,:id), [pi/2;-3pi/4], P = [10.;10.], I = [10.;10], D = [5.;5.])


storage = simulate!(mech,10.,pid,record = true)
visualize!(mech,storage,shapes)