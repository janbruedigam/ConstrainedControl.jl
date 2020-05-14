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


mech = Mechanism(origin, links, constraints, shapes = shapes, g=0.)
setPosition!(mech,origin,link1,Δx = [.1;1.0*0;-1.5*0], Δq = Quaternion(RotY(-1.0*0))*Quaternion(RotX(-1.0*0)))

Q = diagm(ones(12))
Q[7,7] = 0
Q[8,8] = 0
Q[9,9] = 0
Q[10,10] = 0
Q[11,11] = 0
Q[12,12] = 0
R = diagm(ones(6))

steps = Base.OneTo(10)

F0 = [ones(3)*0. for i=steps]
τ0 = [zeros(3) for i=steps]

J0,n,ilqrobj = ilqr(mech, link1.id, steps, Q, R, F0, τ0)

function control!(mechanism,k)
    setForce!(mechanism,mechanism.bodies[1],F=ilqrobj.F1[k],τ=ilqrobj.τ1[k])
end


storage = simulate!(mech,0.1,control!,record = true)
visualize!(mech,storage,shapes)