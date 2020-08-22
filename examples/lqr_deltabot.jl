using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra
using Rotations
using Plots


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
llshape = Box(0.1, 0.1, length1, length1)
ulshape = Box(0.1, 0.1, length1/2, length1/2)
pshape = Box(0.1, 0.1, length1/2*sqrt(2), length1/2*sqrt(2))

pll = [0.0;0.0;length1/2] # joint connection point
pul = [0.0;0.0;length1/4]
pp = [0.0;0.0;length1/4*sqrt(2)]

# Links
origin = Origin{Float64}()
lowerlegl = Body(llshape)
lowerlegr = Body(llshape)
upperlegl = Body(ulshape)
upperlegr = Body(ulshape)
platform = Body(pshape)

# Constraints
floorlr = EqualityConstraint(Revolute(origin, lowerlegl, joint_axis; p2=-pll),Revolute(origin, lowerlegr, joint_axis; p2=-pll),FixedOrientation(origin,platform;qoffset = UnitQuaternion(RotX(pi/2))))
kneel = EqualityConstraint(Revolute(lowerlegl, upperlegl, joint_axis; p1=pll, p2=-pul))
kneer = EqualityConstraint(Revolute(lowerlegr, upperlegr, joint_axis; p1=pll, p2=-pul))
platl = EqualityConstraint(Revolute(platform, upperlegl, joint_axis; p2=pul, p1=pp))
platr = EqualityConstraint(Revolute(platform, upperlegr, joint_axis; p2=pul, p1=-pp))


links = [lowerlegl;lowerlegr;upperlegl;upperlegr;platform]
constraints = [platl;platr;floorlr;kneel;kneer]
shapes = [llshape;ulshape;pshape]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81, Δt = 0.01)
setPosition!(origin,lowerlegl,p2 = -pll,Δq = UnitQuaternion(RotX(pi/4)))
setPosition!(origin,lowerlegr,p2 = -pll,Δq = UnitQuaternion(RotX(-pi/4)))
setPosition!(lowerlegl,upperlegl,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(-pi/2)))
setPosition!(lowerlegr,upperlegr,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(pi/2)))
setPosition!(upperlegl,platform,p1 = pul, p2 = pp,Δq = UnitQuaternion(RotX(3pi/4)))

xd = [links[i].state.xc for i=1:5]
qd=[links[i].state.qc for i=1:5]

Q = [diagm(ones(12))*0.0 for i=1:5]
Q[5][2,2]=10.0
Q[5][3,3]=10.0
Q[5][5,5]=1.0
Q[5][6,6]=1.0
R = [ones(1,1)*0.1 for i=1:2]

lqr = LQR(mech, getid.(links), getid.(constraints[1:2]), Q, R, Inf, xd=xd, qd=qd, Fτd = [[[6.7879484]];[[-6.7879484]]])


### Calculation of initial conditions
s1 = 101
s2 = 101
successful = [false for i=1:s1*s2]
coords = [[0.0;0.0] for i=1:s1*s2]
counter = 1

# Find valid initial conditions
for i = 1:s1
    ppz = -1.5*length1 + 3*length1*(i-1)/(s1-1)
    ppz < 0 && continue
    for j = 1:s2
        global counter

        ppy = -1.5*length1 + 3*length1*(j-1)/(s2-1)
        coords[counter] = [ppy;ppz]

        pp0 = [0;ppy;ppz]
        # display(pp)
        ppr = pp0 + [0;pp[3];0]
        ppl = pp0 - [0;pp[3];0]

        if norm(ppr)<=1.5*length1 && norm(ppl)<=1.5*length1
            if norm(ppr)>=0.5*length1 && norm(ppl)>=0.5*length1
                successful[counter] = true
            end
        end
        counter+=1
    end
end

posres = findall(x->x==true, successful)
A = [getindex.(coords[posres],1) getindex.(coords[posres],2)]

# Find valid initial angles
nc = size(A)[1]
anglesl = [zeros(2) for i=1:nc]
anglesr = [zeros(2) for i=1:nc]

# a: lowerleg
# b: upperleg
# c: platform
for i = 1:nc
    pp0 = A[i,:]
    ppr = pp0 + [pp[3];0]
    ppl = pp0 - [pp[3];0]

    a = length1
    b = length1/2
    
    cl = norm(ppl)
    αl = acos((b^2+cl^2-a^2)/(2*b*cl))
    βl = acos((a^2+cl^2-b^2)/(2*a*cl))
    γl = acos((a^2+b^2-cl^2)/(2*a*b))
    δl = abs(atan(ppl[1]/ppl[2]))

    cr = norm(ppr)
    αr = acos((b^2+cr^2-a^2)/(2*b*cr))
    βr = acos((a^2+cr^2-b^2)/(2*a*cr))
    γr = acos((a^2+b^2-cr^2)/(2*a*b))
    δr = abs(atan(ppr[1]/ppr[2]))

    if ppl[1]<=0 && ppl[2]>=0
        anglesl[i] = [δl+βl;-pi+γl]
    elseif ppl[1]>=0 && ppl[2]>=0
        anglesl[i] = [-δl+βl;-pi+γl]
    elseif ppl[1]>=0 && ppl[2]<=0
        anglesl[i] = [-pi+δl-βl;pi-γl]
    elseif ppl[1]<=0 && ppl[2]<=0
        anglesl[i] = [-pi-δl-βl;pi-γl]
    end

    if ppr[1]<=0 && ppr[2]>=0
        anglesr[i] = [δr-βr;pi-γr]
    elseif ppr[1]>=0 && ppr[2]>=0
        anglesr[i] = [-δr-βr;pi-γr]
    elseif ppr[1]>=0 && ppr[2]<=0
        anglesr[i] = [-pi+δr+βr;-pi+γr]
    elseif ppr[1]<=0 && ppr[2]<=0
        anglesr[i] = [-pi-δr+βr;-pi+γr]
    end
end
### End

i = 750
setPosition!(origin,lowerlegl,p2 = -pll,Δq = UnitQuaternion(RotX(anglesl[i][1])))
setPosition!(origin,lowerlegr,p2 = -pll,Δq = UnitQuaternion(RotX(anglesr[i][1])))
setPosition!(lowerlegl,upperlegl,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(anglesl[i][2])))
setPosition!(lowerlegr,upperlegr,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(anglesr[i][2])))
setPosition!(origin,platform,p1 = [0;A[i,:]],Δq = UnitQuaternion(RotX(pi/2)))

storage = simulate!(mech,10,lqr,record = true)
visualize(mech,storage,shapes)