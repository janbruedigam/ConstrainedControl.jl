using ConstrainedDynamics
using ConstrainedControl
using LinearAlgebra

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
setPosition!(mech,origin,link1,Δx = [1.;2;-3])

Q = diagm(ones(6))
Q[2,2] = 2
Q[3,3] = 0.1
R = diagm(ones(3))
R[1] = 100

lqr = LQR(mech, link1.id, Q, R, xd=[0.;0;0], vd=[0.;0;0.], ud=[0.;0;9.81*box.m])


simulate!(mech,lqr,save = true)
visualize!(mech)