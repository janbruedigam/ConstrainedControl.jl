using ConstrainedDynamics
using ConstrainedControl

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1

p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = 0
q1 = Quaternion(RotX(ϕ1))
ϕ2 = 0
q2 = Quaternion(RotX(ϕ2))

# Links
origin = Origin{Float64}()
link1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))
link2 = deepcopy(link1)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))
joint2 = EqualityConstraint(Revolute(link1, link2, joint_axis; p1=-p2, p2=p2))

links = [link1;link2]
constraints = [joint1;joint2]


mech = Mechanism(origin, links, constraints)
setPosition!(origin,link1,p2 = p2,Δq = q1)
setPosition!(link1,link2,p1=-p2,p2 = p2,Δq = q2)

pid = PID(mech, getfield.(constraints,:id), [pi/2;-pi/4], P = [10.;10.], I = [10.;10], D = [5.;5.])

simulate!(mech,1,pid)
@test true