using ConstrainedControl
using LinearAlgebra
using StaticArrays

X = zeros(8,1000)
N = length(U)
Δt = 0.01


for k=1:1000
    X[1,k] = storage0.x[1][k][2]
    X[2,k] = ConstrainedDynamics.rotation_vector(storage0.q[2][k])[1]
    if X[2,k] > pi
        X[2,k] -= 2pi
    elseif X[2,k] < -pi
        X[2,k] += 2pi
    end
    X[3,k] = ConstrainedDynamics.rotation_vector(storage0.q[3][k])[1] - X[2,k]
    if X[3,k] > pi
        X[3,k] -= 2pi
    elseif X[3,k] < -pi
        X[3,k] += 2pi
    end
    X[4,k] = ConstrainedDynamics.rotation_vector(storage0.q[4][k])[1] - X[2,k] - X[3,k]
    if X[4,k] > pi
        X[4,k] -= 2pi
    elseif X[4,k] < -pi
        X[4,k] += 2pi
    end
    X[5,k] = storage0.v[1][k][2]
    X[6,k] = storage0.ω[2][k][1]
    X[7,k] = storage0.ω[3][k][1] - X[6,k]
    X[8,k] = storage0.ω[4][k][1] - X[6,k] - X[7,k]
end

K = [Float64[0 0 0 0 0 0 0 0] for k=1:N]


Q = [
    10 0 0 0 0 0 0 0
    0 10 0 0 0 0 0 0
    0 0 10 0 0 0 0 0
    0 0 0 10 0 0 0 0
    0 0 0 0 1 0 0 0
    0 0 0 0 0 1 0 0
    0 0 0 0 0 0 1 0 
    0 0 0 0 0 0 0 1
]*Δt
R = ones(1,1)*0.1*Δt
E = SMatrix{8,8,Float64,64}(diagm(ones(8)))

Pk = Q
for k=1000:-1:1
    if mod(k,10)==0
        display(k)
    end
    global Pk
    Pkp1 = Pk

    A = SMatrix{8,8,Float64,64}(Alin.subs([(q[1],X[1,k]),(q[2],X[2,k]),(q[3],X[3,k]),(q[4],X[4,k]),(v[1],X[5,k]),(v[2],X[6,k]),(v[3],X[7,k]),(v[4],X[8,k]),(g,-9.81)]))
    B = SMatrix{8,1,Float64,8}(Blin[:,1].subs([(q[1],X[1,k]),(q[2],X[2,k]),(q[3],X[3,k]),(q[4],X[4,k]),(v[1],X[5,k]),(v[2],X[6,k]),(v[3],X[7,k]),(v[4],X[8,k]),(g,-9.81)]))
    A = A*Δt+E
    B = B*Δt

    K[k] = (R+B'*Pkp1*B)\B'*Pkp1*A
    Abar = A - B*K[k]
    Pk = Q + K[k]'*R*K[k] + Abar'*Pkp1*Abar
end