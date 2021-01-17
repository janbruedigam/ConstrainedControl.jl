using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedControl
using LinearAlgebra


path = "examples/examples_files/sawyer_arm.urdf"
mech = Mechanism(path, floating=false, g = 0.0)

setPosition!(mech,mech.eqconstraints["right_j0"],[0.0])
setPosition!(mech,mech.eqconstraints["right_j1"],[-pi/2])
setPosition!(mech,mech.eqconstraints["right_j2"],[-pi/2])
setPosition!(mech,mech.eqconstraints["right_j3"],[pi/2])

xd=[mech.bodies[i].state.xc for i=1:length(mech.bodies)]
qd=[mech.bodies[i].state.qc for i=1:length(mech.bodies)]

setPosition!(mech,mech.eqconstraints["right_j0"],[0.0])
setPosition!(mech,mech.eqconstraints["right_j1"],[0.0])
setPosition!(mech,mech.eqconstraints["right_j2"],[0.0])
setPosition!(mech,mech.eqconstraints["right_j3"],[0.0])


Q = [diagm(ones(12))*1000.0 for i=1:7]
R = [ones(1,1) for i=1:7]

ConstrainedDynamics.deactivateConstraints!(mech)
lqr = LQR(mech, getid.(mech.bodies), getid.(mech.eqconstraints), Q, R, 20., xd=xd, qd=qd)
ConstrainedDynamics.activateConstraints!(mech)


storage = simulate!(mech,20,lqr,record = true)
visualize(mech,storage)
