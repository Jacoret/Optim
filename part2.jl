# Test de f4
U = [1.0,2.0,3.0,4.0,5.0]
println("f4(U) = $(f4(U))")

# Comparaison des gradients analytiques et approximés
println("df4 : $(df4!(U, zeros(Float64, 5)))")
println("approxdf4 : $(approxdf4!(U, zeros(Float64, 5)))")
