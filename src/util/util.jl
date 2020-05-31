function dare(A,B,Q,R)
    Ait = inv(A)'
    G = B/R*B'

    Z = [A + G*Ait*Q    -G*Ait
         -Ait*Q         Ait]

    S = schur(Z)
    S = ordschur(S, abs.(S.values).<=1)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]

    return U21/U11
end