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

# Desired orientation
ϕ1 = 0
ϕ2 = 0

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


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81, Δt = 0.001)
setPosition!(origin,lowerlegl,p2 = -pll,Δq = UnitQuaternion(RotX(pi/4)))
setPosition!(origin,lowerlegr,p2 = -pll,Δq = UnitQuaternion(RotX(-pi/4)))
setPosition!(lowerlegl,upperlegl,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(-pi/2)))
setPosition!(lowerlegr,upperlegr,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(pi/2)))
setPosition!(upperlegl,platform,p1 = pul, p2 = pp,Δq = UnitQuaternion(RotX(3pi/4)))

xd = [links[i].state.xc for i=1:5]
qd=[links[i].state.qc for i=1:5]

Q = [diagm(ones(12))*0.0 for i=1:5]
# Q[5][2,2]=10.0
# Q[5][3,3]=10.0
# Q[5][5,5]=1.0
# Q[5][6,6]=1.0
# R = [ones(1,1)*0.1 for i=1:2]
Q[5][2,2]=100.0
Q[5][3,3]=100.0
Q[5][5,5]=1.0
Q[5][6,6]=1.0
R = [ones(1,1)*0.01 for i=1:2]

lqr = LQR(mech, getid.(links), getid.(constraints[1:2]), Q, R, Inf, xd=xd, qd=qd, Fτd = [[[6.7879484]];[[-6.7879484]]])

function contr!(mechanism,k)
    setForce!(mechanism,geteqconstraint(mechanism,platl.id),[6.7879484])
    setForce!(mechanism,geteqconstraint(mechanism,platr.id),[-6.7879484])
    return
end

s1 = 101
s2 = 101
# storages = [Storage{Float64}(steps,2) for i=1:s1*s2]
storage = Storage{Float64}(steps,5)
successful = [false for i=1:s1*s2]
coords = [[0.0;0.0] for i=1:s1*s2]
counter = 1

for i = 1:s1
    ppz = -1.5*length1 + 3*length1*(i-1)/(s1-1)
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
# A = [A;[-1.5 -1.5];[-1.5 1.5];[1.5 -1.5];[1.5 1.5]]
# scatter(A[:,1],A[:,2])

nc = size(A)[1]
anglesl = [zeros(2) for i=1:nc]
anglesr = [zeros(2) for i=1:nc]
templ = [zeros(4) for i=1:nc]
tempr = [zeros(4) for i=1:nc]

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

    templ[i] = [αl;βl;γl;δl]
    tempr[i] = [αr;βr;γr;δr]
end


steps = Base.OneTo(25000)
storage = Storage{Float64}(steps,5)
successful2 = [false for i=1:nc]

# i = 30
# setPosition!(origin,lowerlegl,p2 = -pll,Δq = UnitQuaternion(RotX(anglesl[i][1])))
# setPosition!(origin,lowerlegr,p2 = -pll,Δq = UnitQuaternion(RotX(anglesr[i][1])))
# setPosition!(lowerlegl,upperlegl,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(anglesl[i][2])))
# setPosition!(lowerlegr,upperlegr,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(anglesr[i][2])))
# setPosition!(origin,platform,p1 = [0;A[i,:]],Δq = UnitQuaternion(RotX(pi/2)))

# simulate!(mech,steps,storage,lqr, record=true)
# visualize(mech,storage,shapes)


for i=1:nc
    setPosition!(origin,lowerlegl,p2 = -pll,Δq = UnitQuaternion(RotX(anglesl[i][1])))
    setPosition!(origin,lowerlegr,p2 = -pll,Δq = UnitQuaternion(RotX(anglesr[i][1])))
    setPosition!(lowerlegl,upperlegl,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(anglesl[i][2])))
    setPosition!(lowerlegr,upperlegr,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(anglesr[i][2])))
    setPosition!(origin,platform,p1 = [0;A[i,:]],Δq = UnitQuaternion(RotX(pi/2)))
    setVelocity!(lowerlegl)
    setVelocity!(lowerlegr)
    setVelocity!(upperlegl)
    setVelocity!(upperlegr)
    setVelocity!(platform)
    display(string(i)*" von "*string(nc))

    for con in constraints
        con.λsol[1] *= 0
        con.λsol[2] *= 0
    end

    try
        # simulate!(mech,storages[counter],lqr,record = true)
        # # simulate!(mech,storages[counter],contr!,record = true)
        simulate!(mech,steps,storage,lqr)
        # simulate!(mech,steps,storage,contr!)
    catch
        display("Failed for "*string(i))
    end
    if norm(platform.state.xc - [0;0;1.5*length1*sqrt(2)/2])<1e-1 && norm(platform.state.vc)<1e-1
        successful2[i] = true
    end
end

posres = findall(x->x==true, successful2)
A = [A[i,:] for i=1:nc]
A = [getindex.(A[posres],1) getindex.(A[posres],2)]
scatter(A[:,1],A[:,2])