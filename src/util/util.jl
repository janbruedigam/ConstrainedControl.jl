function dare(A, B, Q, R)
    n = size(A, 1)
    T = promote_type(eltype(A), eltype(B), eltype(Q), eltype(R))
    
    E = [
        Matrix{T}(I, n, n)  B/R*B'
        zeros(size(A))      A'
    ]

    F = [
        A   zeros(size(A))
        -Q  Matrix{T}(I, n, n)
    ]
    
    QZ = schur(F, E)
    QZ = ordschur(QZ, abs.(QZ.alpha./QZ.beta) .< 1)
    
    return QZ.Z[(n+1):end, 1:n]/QZ.Z[1:n, 1:n]
end

function care(A, B, Q, R)
    G = B/R*B'

    Z = [A   -G
         -Q  -A']

    S = schur(Z)
    S = ordschur(S, real(S.values) .< 0)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]

    return U21/U11
end