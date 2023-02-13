using ConstrainedDynamics
using ConstrainedControl
using LinearAlgebra


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0

pll = [0.0;0.0;length1/2] # joint connection point
pul = [0.0;0.0;length1/4]
pp = [0.0;0.0;length1/4*sqrt(2)]

# Links
origin = Origin{Float64}()
lowerlegl = Box(0.1, 0.1, length1, length1)
lowerlegr = deepcopy(lowerlegl)
upperlegl = Box(0.1, 0.1, length1/2, length1/2)
upperlegr = deepcopy(upperlegl)
platform = Box(0.1, 0.1, length1/2*sqrt(2), length1/2*sqrt(2))

# Constraints
floorlr = EqualityConstraint(Revolute(origin, lowerlegl, joint_axis; p2=-pll),Revolute(origin, lowerlegr, joint_axis; p2=-pll),FixedOrientation(origin,platform;qoffset = Quaternion(RotX(pi/2))))
kneel = EqualityConstraint(Revolute(lowerlegl, upperlegl, joint_axis; p1=pll, p2=-pul))
kneer = EqualityConstraint(Revolute(lowerlegr, upperlegr, joint_axis; p1=pll, p2=-pul))
platl = EqualityConstraint(Revolute(platform, upperlegl, joint_axis; p2=pul, p1=pp))
platr = EqualityConstraint(Revolute(platform, upperlegr, joint_axis; p2=pul, p1=-pp))


links = [lowerlegl;lowerlegr;upperlegl;upperlegr;platform]
constraints = [platl;platr;floorlr;kneel;kneer]


mech = Mechanism(origin, links, constraints, g=-9.81, Δt = 0.01)
setPosition!(origin,lowerlegl,p2 = -pll,Δq = Quaternion(RotX(pi/4)))
setPosition!(origin,lowerlegr,p2 = -pll,Δq = Quaternion(RotX(-pi/4)))
setPosition!(lowerlegl,upperlegl,p1 = pll, p2 = -pul,Δq = Quaternion(RotX(-pi/2)))
setPosition!(lowerlegr,upperlegr,p1 = pll, p2 = -pul,Δq = Quaternion(RotX(pi/2)))
setPosition!(upperlegl,platform,p1 = pul, p2 = pp,Δq = Quaternion(RotX(3pi/4)))

xd = [links[i].state.xc for i=1:5]
qd=[links[i].state.qc for i=1:5]

Q = [diagm(ones(12))*0.0 for i=1:5]
Q[5][2,2]=10.0
Q[5][3,3]=10.0
Q[5][5,5]=1.0
Q[5][6,6]=1.0
R = [ones(1,1)*0.1 for i=1:2]

lqr = LQR(mech, getid.(links), getid.(constraints[1:2]), Q, R, Inf, xd=xd, qd=qd, Fτd = [[[6.7879484]];[[-6.7879484]]])

simulate!(mech,1,lqr)
@test true