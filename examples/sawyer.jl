using ConstrainedDynamics
using ConstrainedDynamicsVis


path = "examples/examples_files/sawyer_arm.urdf"
mech, shapes = Mechanism(path, floating=false, g = -0.0)

setPosition!(mech,mech.eqconstraints[8],[0.0])
setPosition!(mech,mech.eqconstraints[9],[-pi/2])
setPosition!(mech,mech.eqconstraints[10],[-pi/2])
setPosition!(mech,mech.eqconstraints[11],[pi/2])
setPosition!(mech,mech.eqconstraints[12],[0.0])
setPosition!(mech,mech.eqconstraints[13],[0.0])
setPosition!(mech,mech.eqconstraints[14],[0.0])

xd=[mech.bodies[i].state.xc for i=1:length(mech.bodies)]
qd=[mech.bodies[i].state.qc for i=1:length(mech.bodies)]

setPosition!(mech,mech.eqconstraints[8],[0.0])
setPosition!(mech,mech.eqconstraints[9],[0.0])
setPosition!(mech,mech.eqconstraints[10],[0.0])
setPosition!(mech,mech.eqconstraints[11],[0.0])
setPosition!(mech,mech.eqconstraints[12],[0.0])
setPosition!(mech,mech.eqconstraints[13],[0.0])
setPosition!(mech,mech.eqconstraints[14],[0.0])


Q = [diagm(ones(12))*1000.0 for i=1:7]
R = [ones(1,1) for i=1:7]

for (id,_) in pairs(mech.eqconstraints)
    deactivate!(mech,id)
end
lqr = LQR(mech, getid.(mech.bodies), getid.(mech.eqconstraints), Q, R, 20., xd=xd, qd=qd)
for (id,_) in pairs(mech.eqconstraints)
    activate!(mech,id)
end


storage = simulate!(mech,20.,lqr,record = true)
visualize(mech,storage,shapes)
