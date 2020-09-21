using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra
using Rotations
using Plots


# Parameters
joint_axis = [1.0;0.0;0.0]

l1 = 1.0
l2 = l1/2
l3 = l1/2*sqrt(2)
m1 = l1
m2 = l2
m3 = l3
llshape = Box(0.1, 0.1, l1, m1)
ulshape = Box(0.1, 0.1, l2, m2)
bshape = Box(0.1, 0.1, l3, m3)

pll = [0.0;0.0;l1/2] # joint connection point
pul = [0.0;0.0;l2/2]
pp = [0.0;0.0;l3/2]

# Desired orientation
ϕ1 = 0
ϕ2 = 0

# Links
origin = Origin{Float64}()
lowerlegl = Body(llshape)
lowerlegr = Body(llshape)
upperlegl = Body(ulshape)
upperlegr = Body(ulshape)
base = Body(bshape)

# Constraints
floorlr = EqualityConstraint(Revolute(origin, lowerlegl, joint_axis; p2=-pll),Revolute(origin, lowerlegr, joint_axis; p2=-pll),FixedOrientation(origin,base;qoffset = UnitQuaternion(RotX(pi/2))))
kneel = EqualityConstraint(Revolute(lowerlegl, upperlegl, joint_axis; p1=pll, p2=-pul))
kneer = EqualityConstraint(Revolute(lowerlegr, upperlegr, joint_axis; p1=pll, p2=-pul))
platl = EqualityConstraint(Revolute(base, upperlegl, joint_axis; p2=pul, p1=pp))
platr = EqualityConstraint(Revolute(base, upperlegr, joint_axis; p2=pul, p1=-pp))


links = [lowerlegl;lowerlegr;upperlegl;upperlegr;base]
constraints = [platl;platr;floorlr;kneel;kneer]
shapes = [llshape;ulshape;bshape]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-9.81, Δt = 0.001)
setPosition!(origin,lowerlegl,p2 = -pll,Δq = UnitQuaternion(RotX(pi/4)))
setPosition!(origin,lowerlegr,p2 = -pll,Δq = UnitQuaternion(RotX(-pi/4)))
setPosition!(lowerlegl,upperlegl,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(-pi/2)))
setPosition!(lowerlegr,upperlegr,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(pi/2)))
setPosition!(upperlegl,base,p1 = pul, p2 = pp,Δq = UnitQuaternion(RotX(3pi/4)))

xd = [links[i].state.xc for i=1:5]
qd=[links[i].state.qc for i=1:5]

Q = [diagm(ones(12))*0.0 for i=1:5]
Q[5][2,2]=100.0
Q[5][3,3]=100.0
Q[5][5,5]=1.0
Q[5][6,6]=1.0
R = [ones(1,1)*0.01 for i=1:2]

lqr = LQR(mech, getid.(links), getid.(constraints[1:2]), Q, R, Inf, xd=xd, qd=qd, Fτd = [[[6.7879484]];[[-6.7879484]]])

function contr!(mechanism::Mechanism{T,N,Nb}, k) where {T,N,Nb}
    Δy = base.state.xc[2]
    Δz = base.state.xc[3]-3/4*sqrt(2)
    Δvy = base.state.vc[2]
    Δvz = base.state.vc[3]
    
    u1 = ([75.5398*Δy] + [-65.6468*Δz] + [22.2935*Δvy] + [-18.5727*Δvz])
    u2 = ([75.5398*Δy] + [65.6468*Δz] + [22.2935*Δvy] + [18.5727*Δvz])
    setForce!(mechanism, geteqconstraint(mechanism, 6), u1+[6.7879484])
    setForce!(mechanism, geteqconstraint(mechanism, 7), u2-[6.7879484])

    return
end

s1 = 101
s2 = 101
successful = [false for i=1:s1*s2]
coords = [[0.0;0.0] for i=1:s1*s2]
counter = 1

for i = 1:s1
    ppz = -1.5*l1 + 3*l1*(i-1)/(s1-1)
    for j = 1:s2
        global counter

        ppy = -1.5*l1 + 3*l1*(j-1)/(s2-1)
        coords[counter] = [ppy;ppz]

        pp0 = [0;ppy;ppz]
        # display(pp)
        ppr = pp0 + [0;pp[3];0]
        ppl = pp0 - [0;pp[3];0]

        if norm(ppr)<=1.5*l1 && norm(ppl)<=1.5*l1
            if norm(ppr)>=0.5*l1 && norm(ppl)>=0.5*l1
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
# c: base
for i = 1:nc
    pp0 = A[i,:]
    ppr = pp0 + [pp[3];0]
    ppl = pp0 - [pp[3];0]

    a = l1
    b = l1/2
    
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
# storages = [Storage{Float64}(steps,2) for i=1:s1*s2]
storage = Storage{Float64}(steps,5)
counter = 1
successful2 = [false for i=1:nc]


for i=1:nc
    setPosition!(origin,lowerlegl,p2 = -pll,Δq = UnitQuaternion(RotX(anglesl[i][1])))
    setPosition!(origin,lowerlegr,p2 = -pll,Δq = UnitQuaternion(RotX(anglesr[i][1])))
    setPosition!(lowerlegl,upperlegl,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(anglesl[i][2])))
    setPosition!(lowerlegr,upperlegr,p1 = pll, p2 = -pul,Δq = UnitQuaternion(RotX(anglesr[i][2])))
    setPosition!(origin,base,p1 = [0;A[i,:]],Δq = UnitQuaternion(RotX(pi/2)))
    setVelocity!(lowerlegl)
    setVelocity!(lowerlegr)
    setVelocity!(upperlegl)
    setVelocity!(upperlegr)
    setVelocity!(base)
    display(string(i)*" von "*string(nc))

    for con in constraints
        con.λsol[1] *= 0
        con.λsol[2] *= 0
    end

    try
        # simulate!(mech,storages[counter],lqr,record = true)
        # # simulate!(mech,storages[counter],contr!,record = true)
        simulate!(mech,steps,storage,lqr,record = true)
        # simulate!(mech,steps,storage,contr!,record = true)
    catch
        display("Failed for "*string(i))
    end
    if !any(getindex.(storage.ω[1],1).>100*pi) && !any(getindex.(storage.ω[2],1).>100*pi) && 
        !any(getindex.(storage.ω[3],1).>100*pi) && !any(getindex.(storage.ω[4],1).>100*pi) && 
        !any(getindex.(storage.ω[5],1).>100*pi) && 
        norm(vcat([links[i].state.xc for i=1:5]...,[links[i].state.vc for i=1:5]...) - vcat(xd...,zeros(15)))<1e-1
        successful2[i] = true
    end
end

posres = findall(x->x==true, successful2)
A = [A[i,:] for i=1:nc]
A = [getindex.(A[posres],1) getindex.(A[posres],2)]
scatter(A[:,1],A[:,2])

# visualize(mech,storages[43],shapes)