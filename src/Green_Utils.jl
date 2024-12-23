function de_obj(x::AbstractVector{T}, pdata) where {T}
        nonlinear = pdata.nonlinear
        useC = pdata.useC
        mu=pdata.mu
#        FVout=Fres(x,pdata)
#        FV = FVout.FV
#        GV = FVout.GV
#     FV = Xres(x,pdata)
     FVout = Fres(x,pdata)
     FV = FVout.FV
     GV = FVout.GV
    useC ? sl1n=sL1(GV,mu) : sl1n=sL1(FV,mu)
#    sl1n = sL1(FV,mu)
     sl1x = sL1(x,mu)
    delta = pdata.ADel[1]
    pdata.lasso[1] ? l1part = sl1x : l1part = sl1n
    if delta > 1.e-8
        deob = znorm(FV) + 2.0*delta *l1part
#        deob = znorm(FV) + 2.0 * delta * sl1n
    else 
        deob = znorm(FV)
    end
    return deob
end

function nlf(x)
    nlf = (sin(pi * x) + x^3) / (1.0 + x^2)
    return nlf
end

function Xres(x,pdata)
Fout=Fres(x,pdata)
FV = Fout.FV
return FV
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
        h1n = sqrt(regparm)*TPD * x
        lfact = pi
        nonlinear ? GV = A * nlf.(z) .- b : GV = A * (lfact*z) .- b
        FV = [GV;h1n]
        FVout = (FV=FV, GV=GV)
        return FVout
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
