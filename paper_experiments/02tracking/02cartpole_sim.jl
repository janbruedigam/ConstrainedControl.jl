using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra
using Rotations


# Parameters
ex = [1.0;0.0;0.0]
ey = [0.0;1.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
cartshape = Box(0.1, 0.5, 0.1, length1/2)
poleshape = Box(width, depth, length1, length1)

p2 = [0.0;0.0;length1/2] # joint connection point

# Desired orientation
ϕ = 0.

# Links
origin = Origin{Float64}()
cart = Body(cartshape)
pole1 = Body(poleshape)
pole2 = Body(poleshape)
pole3 = Body(poleshape)

# Constraints
joint1 = EqualityConstraint(Prismatic(origin, cart, ey))
joint2 = EqualityConstraint(Revolute(cart, pole1, ex; p2=p2))
joint3 = EqualityConstraint(Revolute(pole1, pole2, ex; p1 = -p2, p2=p2))
joint4 = EqualityConstraint(Revolute(pole2, pole3, ex; p1 = -p2, p2=p2))

links = [cart;pole1;pole2;pole3]
constraints = [joint1;joint2;joint3;joint4]
shapes = [cartshape;poleshape]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81,Δt = 0.01)
setPosition!(origin,cart,Δx = [0;0.0;0])
setPosition!(cart,pole1,p2 = p2,Δq = UnitQuaternion(RotX(ϕ+0)))
setPosition!(pole1,pole2,p1 = -p2,p2 = p2,Δq = UnitQuaternion(RotX(0.)))
setPosition!(pole2,pole3,p1 = -p2,p2 = p2,Δq = UnitQuaternion(RotX(0.)))

function control!(mechanism, k)
    setForce!(mechanism, geteqconstraint(mechanism,5), [U[k]])
end

steps = Base.OneTo(1000)
storage0 = Storage{Float64}(steps,4)

simulate!(mech,storage0,control!,record = true);
# visualize(mech,storage0,shapes)