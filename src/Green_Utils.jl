function de_obj(x::AbstractVector{T}, pdata) where {T}
if false
        A = pdata.A
        b = pdata.b
        m = pdata.m
        TP = pdata.TP
        TPD = pdata.TPD
        regparm = pdata.regparm
        nonlinear = pdata.nonlinear
        mu = pdata.mu
        z = TP * x
        h1n = regparm*TPD * x
        nonlinear ? GV = A * nlf.(z) .- b : GV = A * z .- b
        FV = [GV;h1n]
else
        mu=pdata.mu
        FV=Fres(x,pdata)
end
    delta = pdata.ADel[1]
    if delta > 1.e-8
        deob = znorm(FV) + 2.0 * delta * sL1(FV, mu)
    else 
        deob = znorm(FV)
    end
    return deob
end

function nlf(x)
    nlf = (sin(pi * x) + x^3) / (1.0 + x^2)
    return nlf
end

function Fres(x,pdata)
        A = pdata.A
        b = pdata.b
        m = pdata.m
        TP = pdata.TP
        TPD = pdata.TPD
        regparm = pdata.regparm
        nonlinear = pdata.nonlinear
        mu = pdata.mu
        z = TP * x
        h1n = regparm*TPD * x
        nonlinear ? GV = A * nlf.(z) .- b : GV = A * z .- b
        FV = [GV;h1n]
return FV
end


function chebmat(m = 100, n = 10)
    h = 1.0 / (m + 1)
    x = h:h:1.0-h
    y = 2.0 * (x .- 0.5)
    indp = zeros(Int64, n)
    T = zeros(m, n)
    for i = 1:n
        indp[i] = 1
        p = ChebyshevT(indp)
        indp[i] = 0
        T[:, i] .= p.(y)
    end
    return T
end
