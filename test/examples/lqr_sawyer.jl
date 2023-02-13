using ConstrainedDynamics
using ConstrainedControl
using LinearAlgebra


path = "examples/examples_files/sawyer_arm.urdf"
mech = Mechanism(path, floating=false, g = 0.0)

setPosition!(mech,mech.eqconstraints["right_j0"],[0.0])
setPosition!(mech,mech.eqconstraints["right_j1"],[0])
setPosition!(mech,mech.eqconstraints["right_j2"],[0])
setPosition!(mech,mech.eqconstraints["right_j3"],[0])

xd=[body.state.xc for body in mech.bodies]
qd=[body.state.qc for body in mech.bodies]

setPosition!(mech,mech.eqconstraints["right_j0"],[0.0])
setPosition!(mech,mech.eqconstraints["right_j1"],[0.0])
setPosition!(mech,mech.eqconstraints["right_j2"],[0.0])
setPosition!(mech,mech.eqconstraints["right_j3"],[0.0])


Q = [diagm(ones(12))*1000.0 for i=1:7]
R = [ones(1,1) for i=1:7]

# ConstrainedDynamics.deactivateConstraints!(mech)
lqr = LQR(mech, getid.(mech.bodies), getid.(mech.eqconstraints), Q, R, 20., xd=xd, qd=qd)
# ConstrainedDynamics.activateConstraints!(mech)

simulate!(mech,1,lqr)
@test true