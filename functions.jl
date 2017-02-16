function f1(U::Vector)
    x = transpose(U) * S * U - transpose(B) * U
	return x[1,1]
end

function df1!(U::Vector, dfU::Vector)
    dfU[:] = (S + transpose(S)) * U - B
    return dfU
end

function approxdf1!(U::Vector, dfU::Vector)
    epsilon = 1e-4
    for i = 1:length(U)
        Ueps = copy(U)
        Ueps[i] *= (1+epsilon)
        dfU[i] = (f1(Ueps)-f1(U))/(epsilon*U[i])
    end
    return dfU
end

function f2(U::Vector)
  y = transpose(U)*S*U + transpose(U)*exp(U)
  return y[1,1]
end

function df2!(U::Vector, dfU::Vector)
    dfU[:] = (S+transpose(S)) * U + (1+U).*exp(U)
    return dfU
end

function approxdf2!(U::Vector, dfU::Vector)
    epsilon = 1e-4
    for i = 1:length(U)
        Ueps = copy(U)
        Ueps[i] *= (1+epsilon)
        dfU[i] = (f2(Ueps)-f2(U))/(epsilon*U[i])
    end
    return dfU
end

function gradient_rho_constant(f::Function, df!::Function, U0::Vector, rho::Float64, tol::Float64)
    # Initialisation
    Unow = copy(U0)
    fnow = f(U0)
    Unext = U0 - rho*df!(Unow, zeros(Float64,5))
    fnext = f(Unext)
    while sqrt(sum((fnext-fnow).*(fnext-fnow))) > tol
        # Actualisation des valeurs
        Unow = copy(Unext)
        fnow = (f(Unow))
        # Nouvelle itération
        Unext = Unow - rho*df!(Unow, zeros(Float64,5))
        fnext = f(Unext)
    end
    return Unext, fnext
end

function gradient_rho_adaptatif(f::Function, df!::Function, U0::Vector, tol::Float64)
    n_test = 1000
    # Initialisation
    Unow = copy(U0)
    fnow = f(U0)
    # Calcul du pas optimal
    h = zeros(Float64,n_test)
    for rho_test in 1:n_test
        h[rho_test] = f(Unow-df!(Unow, zeros(Float64,5)*rho_test/n_test))
    end
    fnext, rho_test = findmin(h)
    rho = rho_test/n_test
    # Première itération
    Unext = U0 - rho*df!(Unow, zeros(Float64,5))
    fnext = f(Unext)
    while sqrt(sum((fnext-fnow).*(fnext-fnow))) > tol
        # Actualisation des valeurs
        Unow = copy(Unext)
        fnow = (f(Unow))
        # Calcul du pas optimal et nouvelle itération
        h = zeros(Float64,n_test)
        for rho_test in 1:n_test
            h[rho_test] = f(Unow-df!(Unow, zeros(Float64,5)*rho_test/n_test))
        end
        fnext, rho_test = findmin(h)
        rho = rho_test/n_test
        Unext = Unow - rho*df!(Unow, zeros(Float64,5))
    end
    return Unext, fnext
end
