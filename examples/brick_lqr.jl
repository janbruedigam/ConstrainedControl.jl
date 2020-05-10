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


mech = Mechanism(origin, links, constraints, shapes = shapes, tend = 10.,g=0.)
setPosition!(mech,origin,link1,Δx = [1.;2;3])

Q = diagm(ones(6))
R = diagm(ones(3))
R[1] = 100

lqr = LQR(mech, link1.id, [0.], [0.], [0.], Q, R)


simulate!(mech,lqr,save = true)
visualize!(mech)