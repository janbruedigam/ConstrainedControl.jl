using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra
using Rotations


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1

p2a = [0.0;0.0;length1/2] # joint connection point
p2b = [0.0;0.0;length1] # joint connection point

# Desired orientation
ϕ1 = pi
ϕ2 = pi

# Links
origin = Origin{Float64}()
link1 = Box(width, depth, length1, length1)
link2 = Box(width, depth, length1*2, length1)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2a))
joint2 = EqualityConstraint(Revolute(link1, link2, joint_axis; p1=-p2a, p2=p2b))

links = [link1;link2]
constraints = [joint1;joint2]


mech = Mechanism(origin, links, constraints, g=-9.81)
setPosition!(origin,link1,p2 = p2a,Δq = QuatRotation(RotX(ϕ1-0.1)))
setPosition!(link1,link2,p1=-p2a,p2 = p2b,Δq = QuatRotation(RotX(0.1)))

xd = [[[0;0;0.5]];[[0;0;2.0]]]
qd=[[QuatRotation(RotX(ϕ1))];[QuatRotation(RotX(ϕ2))]]

Q = [diagm(ones(12))*0.0 for i=1:2]
Q[1][7,7]=4.0
Q[1][10,10]=4.0
Q[2][7,7]=1.0
Q[2][10,10]=1.0
R = [ones(1,1)]

lqr = LQR(mech, getid.(links), [getid(constraints[2])], Q, R, 10., xd=xd, qd=qd)


storage = simulate!(mech,10,lqr,record = true)
visualize(mech,storage)