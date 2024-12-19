"""
PlotNew()

Make the plots in the paper.

PlotNew() makes figure 1
PlotNew(;longplot=true) makes figure 2
"""
function PlotNew(
    m = 1000,
    k = 1;
    noiselev = 0.0,
    delta = 0.1,
    nonlinear = true,
    mu = 1.e-8,
# regparm > 0 -> H1 regularization
    regparm=0.01,
    useC = true,
# longplot = false -> figure 1
# longplot = true -> figure 2
    longplot = false,
    lasso = [false]
)
lassoflag = lasso[1]
(larray, dran) = PlotSet(lassoflag, nonlinear, longplot)
    pdata = SetDEData(m, noiselev, delta, nonlinear, 1.e-8, regparm; 
            useC=useC, lasso=lasso)
    DeltaRec=zeros(length(dran),length(larray))
# Setting delta = 0 means you're getting the least squares (non-robust)
# solution. 
    DLen=length(dran)
    LLen=length(larray)
if ~lassoflag
    pdata.ADel[1] = 0.0
    LSout = G2(pdata)
    LSRes = LSout.NLRes.FV
    m=length(LSRes); mx=Int(m/2)
    useC ? (CLSRes = LSRes[1:mx]) : (CLSRes = LSRes)
end
#    pdata.x0 .= LSout.xval
    for ilam = 1:LLen
    pdata.ADel[1] = larray[ilam]
if lassoflag
# Lasso solution?
(LSRes, CLSRes) = Lasso_Base(pdata, lassoflag)
end
#
# RLS solutioin
    pdata.lasso[1]=false
    RBout = G2(pdata)
#
    pdata.x0 .= RBout.xval
    RBRes = RBout.NLRes.FV
    m=length(LSRes); mx=Int(m/2)
    useC ? (CRBRes = RBRes[1:mx]) : (CRBRes = RBRes)
#    RBRes = RBout.NLRes
#
# The only way either of the conditions for a warning could be
# violated is if the robust solution and the LS solution are near
# different local LS minimia. This can happen if you turn regularization
# off, so do not do that.
#
 ( norm(RBRes) < norm(LSRes) ) && @warn("norm error")
 ( norm(CLSRes,1) < norm(CRBRes,1) ) && @warn("L1 norm error")
    for il = 1:DLen
    delta=dran[il]
    STEP = One_Shot2(k, LSRes, RBRes, delta, pdata)
    DeltaRec[il,ilam]=STEP.DRB
    end
    end
    plot(dran,DeltaRec[:,1],"k--")
    plot(dran,DeltaRec[:,2],"k-")
    plot(dran,DeltaRec[:,3],"k-.")
    plot(dran,DeltaRec[:,4],"k.")
    plot(dran,DeltaRec[:,5],"k:")
    delta1=larray[1]
    delta2=larray[2]
    delta3=larray[3]
    delta4=larray[4]
    delta5=larray[5]
    zref=zeros(DLen)
    plot(dran,zref,"k_")
#    legend([L"$\lambda = \$delta1$",L"$\lambda = .5$"])
    legend([L"$\lambda = $"*"$delta1",
            L"$\lambda = $"*"$delta2",
            L"$\lambda = $"*"$delta3",
            L"$\lambda = $"*"$delta4",
            L"$\lambda = $"*"$delta5" ])
    xlabel(L"$\delta$")
    ylabel(L"$\Delta$")
end



"""
SetDEData(m, noiselev, delta, nonlinear, mu, regparm; 
                   xinit = [], obj=de_obj, useC=false)

Initialize the data structure and preallocate as much as possible.

de_obj is the objective function for this problem.
"""
function SetDEData(m, noiselev, delta, nonlinear, mu, regparm; 
                   xinit = [], obj=de_obj, useC=false, lasso=[false])
    LD1 = D1(m)
    G = GmatRL(m + 2)
    AF = G[2:m+1, 2:m+1]
    dx = 1.0 / (m + 1)
    xg = collect(dx:dx:1.0-dx)
# Make the "exact" solution a challenge to approximate with polynomials.
    xa=(xg .- .25); xe=abs.(xa);
    if nonlinear
        b = AF * xe
    else
        xnl = nlf.(xe)
        b = AF * xnl
    end
# There is no noise in the example in the paper, so noiselev = 0.0
    if noiselev > 1.e-8
         b .+= norm(b,Inf)*noiselev*rand(m)
    end
# Build the Chebychev polynomial matrix
    TP = chebmat(m, 10)
    TPD = LD1*TP
    A = AF
# If x0 hasn't been set by the caller, do something reasonable. Setting
# x0 to the zero vector makes the Jacobian singular, and is not a good
# idea.
    linit = length(xinit)
    if linit == 0
        (ma, na) = size(AF)
        (mt, nt) = size(TP)
        (nt == 1) ? (x0 = ones(na) / sqrt(na)) : (x0 = ones(nt) / sqrt(nt))
    else
        x0 = xinit
    end
# We store delta in an array with one element so we can change it without
# rebuilding the entire structure.
# Tuples in Juila are immutable so I can't just put a scalar in there.
    ADel = zeros(1)
    ADel[1] = delta
    pdata = (
        m = m,
        A = A,
        b = b,
        mu = mu,
        TP = TP,
        TPD = TPD,
        nonlinear = nonlinear,
        x0 = x0,
        ADel = ADel,
        xe = xe,
        obj = obj,
        regparm = regparm,
        useC = useC,
        lasso = lasso
    )
    return pdata
end

"""
G2(pdata)

Take pdata and solve the optimization problem.
Return the output from the optimizer and the data you need
to make the plots and tables.
"""
function G2(pdata)
# When do_opt sees delta = 0 you get the least squares solution
    optout = do_opt(pdata)
    TP = pdata.TP
    xe = pdata.xe
    useC = pdata.useC
    (mt, nt) = size(TP)
    if mt > 1
        zval = TP * optout.xval
    else
        zval = optout.xval
    end
    fval = optout.fval
    sol = optout.sol
    gradnorm = optout.gradnorm
    success = optout.success
    delx = norm(zval - xe, Inf)
#    NLRes = Fres(optout.xval, pdata)
    GLRes = Fres(optout.xval, pdata)
#    NLRes = Xres(optout.xval,pdata)
#    println(norm(GLRes.FV - NLRes,Inf))
#    NLRes = GLRes.FV
    return (
        xval = optout.xval,
        fval = fval,
        gradnorm = gradnorm,
        delx = delx,
        success = success,
        sol = sol,
#        NLRes = NLRes,
        NLRes = GLRes,
    )
end


function One_Shot2(k, LSRes, RBRes, delta, pdata)
m=length(LSRes); mx=Int(m/2); 
useC = pdata.useC
useC ? (CLSRes=LSRes[1:mx]) : (CLSRes=LSRes)
useC ? (CRBRes=RBRes[1:mx]) : (CRBRes=RBRes)
LSErr = norm(LSRes).^2 + 2.0*delta*norm(CLSRes,1)
RBErr = norm(RBRes).^2 + 2.0*delta*norm(CRBRes,1)
DRB = LSErr - RBErr
DRL = LSErr
DRR = RBErr
    return (
        DRB = DRB,
        DRL = DRL,
        DRR = DRR,
    )
end

function GmatRL(n,T=Float64)
    onet=T(1.0)
    h = onet / (n-onet)
    Ns = collect(0:1:n-1)
    X = h*Ns
    X[n]=onet
    G = [greens(x, y, onet) for x in X, y in X]
    G .*= h
    return G
end

function greens(x, y, onet)
    if x > y
        gf = y * (onet - x)
    else
        gf = x * (onet - y)
    end
    return gf
end

function PlotSet(lassoflag, nonlinear, longplot)
if nonlinear
    larray=[0.0, .1, 1.0, 5.0, 10.0, 100.0]
    longplot ?  dran = collect(0.0:0.025:5.0) : dran = collect(0.0:.1:1.0)
end
if ~nonlinear
    larray=[0.0, .5, 1.0, 10.0, 100.0, 200.0]
    longplot ?  dran = collect(0.0:1.0:50.0) : dran = collect(0.0:1.0:20.0)
    if lassoflag
    larray=[0.0, .5, 1.0, 10.0, 100.0, 200.0]
    longplot ?  dran = collect(0.0:1.0:50.0) : dran = collect(0.0:.1:1.0)
    end
end
return (larray, dran)
end

function Lasso_Base(pdata, lassoflag)
lassoflag && (pdata.lasso[1]=true)
useC=pdata.useC
LSout=G2(pdata)
LSRes = LSout.NLRes.FV
m=length(LSRes); mx=Int(m/2)
useC ? (CLSRes = LSRes[1:mx]) : (CLSRes = LSRes)
return (LSRes, CLSRes)
end
