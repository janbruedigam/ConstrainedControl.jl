using ConstrainedDynamics
using ConstrainedControl
using LinearAlgebra
using Rotations


# Parameters
ex = [1.0;0.0;0.0]
ey = [0.0;1.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1

p2 = [0.0;0.0;length1/2] # joint connection point

# Desired orientation
ϕ = 0.

# Links
origin = Origin{Float64}()
cart = Box(0.1, 0.5, 0.1, length1/2)
pole = Box(width, depth, length1, length1)

# Constraints
joint1 = EqualityConstraint(Prismatic(origin, cart, ey))
joint2 = EqualityConstraint(Revolute(cart, pole, ex; p2=-p2))

links = [cart;pole]
constraints = [joint1;joint2]


mech = Mechanism(origin, links, constraints, g=-9.81)
setPosition!(origin,cart,Δx = [0;0.5;0])
setPosition!(cart,pole,p2 = -p2,Δq = UnitQuaternion(RotX(ϕ+0.2)))

xd = [[[0;0;0.0]];[[0;0;0.5]]]

Q = [diagm(ones(12))*1.0 for i=1:2]
R = [ones(1,1)]

lqr = LQR(mech, getid.(links), [getid(constraints[1])], Q, R, 10., xd=xd)
@test true