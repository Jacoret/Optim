using Optim
include("initializations.jl")
include("functions.jl")

function posit(x::Float64)
    if x > 0
        return x
    else
        return 0
    end
end

posit(-0.5)

function penal(U::vector)
    if
