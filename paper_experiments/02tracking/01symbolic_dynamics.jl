using RigidBodyDynamics
using StaticArrays
using SymPy
using LinearAlgebra

inertias = @syms mc mp IC positive = true
lengths = @syms lc lp cc cp real = true
I_p = 1/3*mp*lp^2

gravitational_acceleration = @syms g real = true
params = [inertias..., lengths..., gravitational_acceleration...]
transpose(params)

T = Sym # the 'scalar type' of the Mechanism we'll construct
axis = SA{T}[1;0;0] # axis of rotation for each of the joints
praxis = SA{T}[0;1;0] # prismatic axis
cartpole = Mechanism(RigidBody{T}("world"); gravity = SA{T}[0;0;g])
world = root_body(cartpole) # the fixed 'world' rigid body

inertia1 = SpatialInertia(CartesianFrame3D("cart_frame"),
    moment = IC * axis * transpose(axis),
    com = SA{T}[0;0;0],
    mass = mc)
body1 = RigidBody("cart",inertia1)
joint1 = Joint("slider_to_cart", Prismatic(praxis))
joint1_to_world = one(Transform3D{T}, frame_before(joint1), default_frame(world));
attach!(cartpole, world, body1, joint1, joint_pose = joint1_to_world)

inertia2 = SpatialInertia(CartesianFrame3D("cart_to_pole1_frame"),
    moment=I_p * axis * transpose(axis),
    com=SVector(zero(T), zero(T), cp),
    mass=mp)
body2 = RigidBody("pole1",inertia2)
joint2 = Joint("cart_to_pole1", Revolute(axis))
joint2_to_body1 = Transform3D(frame_before(joint2), default_frame(body1), SVector(zero(T), zero(T), zero(T)))
attach!(cartpole, body1, body2, joint2, joint_pose = joint2_to_body1)

inertia3 = SpatialInertia(CartesianFrame3D("pole_to_pole2_frame"),
    moment=I_p * axis * transpose(axis),
    com=SVector(zero(T), zero(T), cp),
    mass=mp)
body3 = RigidBody("pole2",inertia3)
joint3 = Joint("pole1_to_pole2", Revolute(axis))
joint3_to_body2 = Transform3D(frame_before(joint3), default_frame(body2), SVector(zero(T), zero(T), lp))
attach!(cartpole, body2, body3, joint3, joint_pose = joint3_to_body2)

inertia4 = SpatialInertia(CartesianFrame3D("pole2_to_pole3_frame"),
    moment=I_p * axis * transpose(axis),
    com=SVector(zero(T), zero(T), cp),
    mass=mp)
body4 = RigidBody("pole3",inertia4)
joint4 = Joint("pole2_to_pole3", Revolute(axis))
joint4_to_body3 = Transform3D(frame_before(joint4), default_frame(body3), SVector(zero(T), zero(T), lp))
attach!(cartpole, body3, body4, joint4, joint_pose = joint4_to_body3)

x = MechanismState(cartpole)
q = configuration(x)
for i in eachindex(q)
    q[i] = symbols("q$i", real = true)
end
v = velocity(x)
for i in eachindex(v)
    v[i] = symbols("v$i", real = true)
end

M = simplify.(mass_matrix(x))
M = SMatrix{size(M)...,Sym,length(M)}(M)
B = simplify.(dynamics_bias(x))
B = SVector{length(B),Sym}(B)


invMsub = inv(M).subs([(mc,0.5),(mp,1.0),(lc,0.5),(lp,-1.0),(cc,0.0),(cp,-0.5)])
Bsub = B.subs([(mc,0.5),(mp,1.0),(lc,0.5),(lp,-1.0),(cc,0.0),(cp,-0.5)])

u = symbols("u")

#qdd = -invMsub*Bsub + invMsub*[u;0;0;0]
Anl = -invMsub*Bsub

A1 = diff.(Anl,q[1])
A2 = diff.(Anl,q[2])
A3 = diff.(Anl,q[3])
A4 = diff.(Anl,q[4])
A5 = diff.(Anl,v[1])
A6 = diff.(Anl,v[2])
A7 = diff.(Anl,v[3])
A8 = diff.(Anl,v[4])

A1 = diff.(Anl,q[1])
A2 = diff.(Anl,q[2])
A3 = diff.(Anl,q[3])
A4 = diff.(Anl,q[4])
A5 = diff.(Anl,v[1])
A6 = diff.(Anl,v[2])
A7 = diff.(Anl,v[3])
A8 = diff.(Anl,v[4])



Z = zeros(4,4)
E = diagm(ones(4))
Alin = [
    Z E  
    A1 A2 A3 A4 A5 A6 A7 A8
]
Blin = [
    Z
    invMsub
]

Alin.subs([(q[1],0),(q[2],0),(q[3],0),(q[4],0),(v[1],0),(v[2],0),(v[3],0),(v[4],0)])
Blin.subs([(q[1],0),(q[2],0),(q[3],0),(q[4],0),(v[1],0),(v[2],0),(v[3],0),(v[4],0)])

# Z = zeros(2,2)
# E = diagm(ones(2))
# Alin = [
#     Z E  
#     A1 A2 A5 A6
# ]
