using Optim
include("initializations.jl")
include("functions.jl")

println("\n Start of the Program ... \n")


println("df1 : $(df1!(Umin, zeros(Float64, 5)))")
println("approxdf1 : $(approxdf1!(Umin, zeros(Float64, 5)))")

df2_Umin = df2!(Umin, zeros(Float64, 5))
println("df2 en Umin : $df2_Umin")
approxdf2_Umin = approxdf2!(Umin, zeros(Float64, 5))
println("approxdf2 en Umin : $approxdf2_Umin")

xmin, fmin = gradient_rho_constant(f1, df1!, U0, 0.001, 1e-4)
println("xmin : $xmin, fmin : $fmin")
println("f1(Umin) = $(f1(Umin)) , f2(Umin) = $(f2(Umin))")

xmin, fmin = gradient_rho_adaptatif(f1, df1!, U0, 1e-4)
println("xmin : $xmin, fmin : $fmin")
println("f1(Umin) = $(f1(Umin)) , f2(Umin) = $(f2(Umin))")

res = optimize(f1, df1!, U0, BFGS())
res = optimize(f1, df1!, U0, ConjugateGradient())
res = optimize(f1, df1!, U0, GradientDescent())
