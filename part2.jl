using Optim
include("initializations.jl")
include("functions.jl")
include("recuit.jl")

# Test de f4
U = [1.0,2.0,3.0,4.0,5.0]
println("f4(U) = $(f4(U))")

# Comparaison des gradients analytiques et approximés
println("df4 : $(df4!(U, zeros(Float64, 5)))")
println("approxdf4 : $(approxdf4!(U, zeros(Float64, 5)))")

## b) Méthode d'optimisation classique (BFGS) ----------
optimize(f4, df4!, U0, BFGS())
# Minimum : -0.108

optimize(f4, df4!, U, BFGS())
# Minimum : 1814.47

optimize(f4, df4!, B, BFGS())
# Minimum : 75.03

# On obtient bien différents minima locax

# c) Methode du recuit simulé ----------------
recuit_nm(f4, U, 0.01, 100.0, 0.001, 0.999, 10)
# Le minimum change à chaque appel de la fonction
recuit_nm(f4, U, 0.01, 100.0, 0.001, 0.999, 10)
recuit_nm(f4, U0, 0.01, 100.0, 0.001, 0.999, 10)
recuit_nm(f4, U0, 0.01, 100.0, 0.001, 0.999, 10)
# En ralentissant la descente, on tombe toujours sur le même minimum, qui semble donc être le minimum global
