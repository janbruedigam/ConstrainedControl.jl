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
joint2 = EqualityConstraint(Revolute(cart, pole1, ex; p2=-p2))
joint3 = EqualityConstraint(Revolute(pole1, pole2, ex; p1 = p2, p2=-p2))
joint4 = EqualityConstraint(Revolute(pole2, pole3, ex; p1 = p2, p2=-p2))

links = [cart;pole1;pole2;pole3]
constraints = [joint1;joint2;joint3;joint4]
shapes = [cartshape;poleshape]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81,Δt = dt/10)
setPosition!(origin,cart,Δx = [0;0.0;0])
setPosition!(cart,pole1,p2 = -p2,Δq = UnitQuaternion(RotX(ϕ+pi)))
setPosition!(pole1,pole2,p1 = p2,p2 = -p2,Δq = UnitQuaternion(RotX(0.)))
setPosition!(pole2,pole3,p1 = p2,p2 = -p2,Δq = UnitQuaternion(RotX(0.)))

function control!(mechanism, k)
    setForce!(mechanism, geteqconstraint(mechanism,5), U[k])
end

steps = Base.OneTo(10000)
storage0 = Storage{Float64}(steps,4)

simulate!(mech,storage0,control!,record = true)
visualize(mech,storage0,shapes)

steps = Base.OneTo(1000)
storage0red = Storage{Float64}(steps,4)
Ured = [SA_F64[0] for i=1:1000]


for n=1:4
    for k=1:1000
        storage0red.x[n][k] = storage0.x[n][10*k]
        storage0red.v[n][k] = storage0.v[n][10*k]
        storage0red.q[n][k] = storage0.q[n][10*k]
        storage0red.ω[n][k] = storage0.ω[n][10*k]
    end
end

for k=1:1000
    Ured[k] = U[10*k]
end