using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)

p2 = [0.0;0.0;length1/2] # joint connection point
ϕ = 0.

# Initial orientation

# Links
origin = Origin{Float64}()
link1 = Body(box)
link2 = Body(box)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))
joint2 = EqualityConstraint(Revolute(link1, link2, joint_axis; p1=-p2, p2=p2))

links = [link1;link2]
constraints = [joint1;joint2]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=0)
setPosition!(origin,link1,p2 = p2,Δq = UnitQuaternion(RotX(ϕ+pi+1.0)))
setPosition!(link1,link2,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(ϕ-0.5)))

xd=[[p2];[3*p2]]
qd=[[UnitQuaternion(RotX(ϕ+pi))];[UnitQuaternion(RotX(ϕ+pi))]]


Q = [diagm(ones(12))*10.0 for i=1:2]
R = [ones(1,1) for i=1:2]

for (id,_) in pairs(mech.eqconstraints)
    deactivate!(mech,id)
end
lqr = LQR(mech, getid.(links), getid.(constraints), Q, R, 10., xd=xd, qd=qd)
for (id,_) in pairs(mech.eqconstraints)
    activate!(mech,id)
end

storage = simulate!(mech,10.,lqr,record = true)
visualize(mech,storage,shapes)