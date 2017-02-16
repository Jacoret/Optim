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

function penal(U::Vector)
    b = 0
    for i in 1:5
        if U[i]<0
            b = posit(-U[i])
        elseif U[i]>1
            b = posit(U[i]-1)
        end
    end
    return b
end


penal([0.,1.,10.,0.,.5])


function f1penal(U::Vector, epsilon::Float64)
    f = f1(U) + b(u)/epsilon
    return f
end

function f2penal(U::Vector, epsilon::Float64)
    f = f2(U) + b(u)/epsilon
    return f
end
