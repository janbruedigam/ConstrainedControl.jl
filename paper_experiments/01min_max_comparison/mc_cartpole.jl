using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra
using Rotations
using Plots


# Parameters
ex = [1.0;0.0;0.0]
ey = [0.0;1.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
cartshape = Box(0.1, 0.5, 0.1, length1/2)
poleshape = Box(width, depth, length1, length1)

p2 = [0.0;0.0;length1/2] # joint connection point

# Desired orientation
ϕ = pi

# Links
origin = Origin{Float64}()
cart = Body(cartshape)
pole = Body(poleshape)

# Constraints
joint1 = EqualityConstraint(Prismatic(origin, cart, ey))
joint2 = EqualityConstraint(Revolute(cart, pole, ex; p2=p2))

links = [cart;pole]
constraints = [joint1;joint2]
shapes = [cartshape;poleshape]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81, Δt = 0.001)
# setPosition!(origin,cart,Δx = [0;0.5;0])
# setPosition!(cart,pole,p2 = p2,Δq = UnitQuaternion(RotX(ϕ+0.2)))

xd = [[[0;0;0.0]];[[0;0;0.5]]]
qd=[[UnitQuaternion(RotX(0.0))];[UnitQuaternion(RotX(ϕ))]]

Q = [diagm(ones(12))*0.0 for i=1:2]
Q[1][2,2]=1.0
Q[1][5,5]=1.0
Q[2][7,7]=4.0
Q[2][10,10]=1.0
R = [ones(1,1)]

lqr = LQR(mech, getid.(links), [getid(constraints[1])], Q, R, Inf, xd=xd, qd=qd)

function contr!(mechanism::Mechanism{T,N,Nb}, k) where {T,N,Nb}
    Δz = zeros(T,Nb*12)

    mincoords = (minimalCoordinates(mechanism)[3][1],minimalCoordinates(mechanism)[4][1])
    if mincoords[2]<0
        mincoords = (mincoords[1],mincoords[2]+2pi)
    end
    vels = [mechanism.bodies[1].state.vc[2];mechanism.bodies[2].state.ωc[1]]
    # display(mincoords)
    # display(vels)

    
    u = ([1.0*(mincoords[1])] + [-38.07*(mincoords[2]-pi)] + [2.4*vels[1]] + [-7.83*vels[2]])
    setForce!(mechanism, geteqconstraint(mechanism, 3), u)

    return
end

steps = Base.OneTo(25000)

s1 = 101
s2 = 101
# storages = [Storage{Float64}(steps,2) for i=1:s1*s2]
storage = Storage{Float64}(steps,2)
successful = [false for i=1:s1*s2]
coords = [[0.0;pi] for i=1:s1*s2]
counter = 1

for i=1:s1
    for j=1:s2
        global counter
        setPosition!(origin,cart,Δx = [0;-1+2*(i-1)/(s1-1);0])
        setPosition!(cart,pole,p2 = p2,Δq = UnitQuaternion(RotX(2pi*(j-1)/(s2-1))))   
        setVelocity!(cart)
        setVelocity!(pole)
        coords[counter] = [-1+2*(i-1)/(s1-1);2pi*(j-1)/(s2-1)]
        display(string(i)*string(j))

        joint1.λsol[1] *= 0
        joint1.λsol[2] *= 0
        joint2.λsol[1] *= 0
        joint2.λsol[2] *= 0
        try
            # # simulate!(mech,storages[counter],lqr,record = true)
            # simulate!(mech,storages[counter],contr!,record = true, debug=false)
            # simulate!(mech,steps,storage,lqr)
            simulate!(mech,steps,storage,contr!)
        catch
            display("Failed for "*string(i)*" and "*string(j))
        end
        if norm(cart.state.xc)<1e-1 && norm(pole.state.xc - [0;0;0.5])<1e-1 && norm(cart.state.vc)<1e-1 && norm(pole.state.ωc)<1e-1
            successful[counter] = true
        end
        counter += 1
    end
end

posres = findall(x->x==true, successful)
A = [getindex.(coords[posres],1) getindex.(coords[posres],2)]
scatter(A[:,1],A[:,2])

# visualize(mech,storages[43],shapes)