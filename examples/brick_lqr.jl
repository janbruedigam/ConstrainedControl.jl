using ConstrainedDynamics
using ConstrainedControl
using LinearAlgebra
using Rotations

length1 = 1.0
width, depth = 1.0, 1.0
box = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))

# Links
origin = Origin{Float64}()
link1 = Body(box)

# Constraints
joint1 = EqualityConstraint(OriginConnection(origin, link1))

links = [link1]
constraints = [joint1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes, tend = 20.)
setPosition!(mech,origin,link1,Δx = [.5;1.;-1.5], Δq = Quaternion(RotY(-1.))*Quaternion(RotX(-1.)))

Q = diagm(ones(12))
Q[2,2] = 2
Q[3,3] = 0.1
Q[7,7] = 2
R = diagm(ones(6))
R[1] = 100


lqr,A,B = LQR(mech, link1.id, Q, R, Fd=[0.;0;9.81*box.m])


simulate!(mech,lqr,save = true)
visualize!(mech)