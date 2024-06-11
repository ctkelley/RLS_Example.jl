function znorm(x::AbstractVector{T}) where {T}
    znorm = sum(x .* x)
    return znorm
end 

function sL1(x::AbstractVector{T}, mu = 1.e-7) where {T}
    snrm = sum(Cabs.(x, mu))
    return snrm
end 
    
function Cabs(x, mu)
    p = 4.0 * mu * mu
    Cabs = sqrt(x * x + p)
    return Cabs
end

function Gen_Y(m, k, delta = 1.0)
    T1 = delta*rand(m, k)
    T1
end

function XFres(x, pdata)
delta=pdata.ADel[1]
pdata.ADel[1]=0.0
FV = de_obj(x, pdata)
pdata.ADel[1]=delta
return FV
end


function do_opt(pdata)
    m = pdata.m
    A = pdata.A
#    TP = pdata.TP
    x0 = pdata.x0
    nl_obj=pdata.obj
    delta = pdata.ADel[1]
    optrobde = OptimizationFunction(nl_obj, Optimization.AutoZygote())
#    optrobde = OptimizationFunction(nl_obj, Optimization.AutoEnzyme())
    prob = OptimizationProblem(optrobde, x0, pdata)
    sol = solve(prob, BFGS(); g_tol = 1.e-14)
    #sol=solve(prob, Newton(); g_tol = 1.e-12)
    fval = sol.objective
    xval = sol.u
    gradtup = gradient(nl_obj, xval, pdata)
    grad = gradtup[1]
    gradnorm = norm(grad)
    success = SciMLBase.successful_retcode(sol)
    return (xval = xval, fval = fval, gradnorm = gradnorm, success = success, sol)
end

