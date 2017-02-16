using Optim
include("initializations.jl")
include("functions.jl")

println("\n Start of the Program ... \n")

#--------------------------------------------------------
#--------------------------------------------------------

# Question 1 --------------------------------------------
# a. b. c. : cf functions.jl
# Quelques tests :
U = [1.0,2.0,3.0,4.0,5.0]

# f1 & f2
println("f1(U) = $(f1(U)) , f2(U) = $(f2(U))")

# df1 analytique et approximé
println("df1 : $(df1!(U, zeros(Float64, 5)))")
println("approxdf1 : $(approxdf1!(U, zeros(Float64, 5)))")

# df2 analytique et approximé
println("df2 : $(df2!(U, zeros(Float64, 5)))")
println("approxdf1 : $(approxdf2!(U, zeros(Float64, 5)))")


# Question 2 --------------------------------------------
# cf functions.jl, le vecteur rho est dans initializations.jl
# Quelques tests :

# Méthode du Gradient à pas constant :
# Pour f1 :
for i=1:5
    umin,fmin = gradient_rho_constant(f1, df1!, U0, rho[i], 1e-4)
    println("u1min = $umin, f1min = $fmin")
end

# Pour f2 :
for i=1:5
    umin,fmin = gradient_rho_constant(f2, df2!, U0, rho[i], 1e-4)
    println("u2min = $umin, f2min = $fmin")
end

# Méthode du gradient adaptatif
# pour f1
u1min, f1min = gradient_rho_adaptatif(f1, df1!, U0, 1e-4)
println("u1min : $u1min, f1min : $f1min")

# pour f2
u2min, f2min = gradient_rho_adaptatif(f2, df2!, U0, 1e-4)
println("u2min : $u2min, f2min : $f2min")

# Question 2 bis -----------------------------------------
@time optimize(f1, df1!, U0, BFGS())
@time resBFGS = optimize(f1, df1!, U0, BFGS())
@time optimize(f1, df1!, U0, ConjugateGradient())
@time resCG = optimize(f1, df1!, U0, ConjugateGradient())
@time optimize(f1, df1!, U0, GradientDescent())
@time resGD =optimize(f1, df1!, U0, GradientDescent())

# Les algorithmes BFGS et Conjugate gradient, la descente de gradient est ~100 fois plus lente,
# ce qui ne parait incohérent.

println("BFGS : $resBFGS")
println("Conjugate Gradient : $resCG")
println("Gradient Descent : $resGD")

# Question 3 --------------------------------------------
# Avec un algorithme de gradient stochastique ?
