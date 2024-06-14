function T2(
    m = 1000,
    k = 10000;
    noiselev = 0.0,
    delta = 0.1,
    nonlinear = false,
    mu = 1.e-8,
    regparm=.01
)
println(delta)
    pdata = SetDEData(m, noiselev, delta, nonlinear, mu, regparm)
    pdata.ADel[1] = 0.0
    LSout = G2(pdata)
    LSRes = LSout.NLRes
    pdata.ADel[1] = delta
    RBout = G2(pdata)
    pdata.x0 .= LSout.xval
    RBout=G2(pdata)
    RBRes = RBout.NLRes
    #println(norm(LSout.xval - RBout.xval, 2))
    STEP = One_Shot2(k, LSRes, RBRes, delta, pdata)
    println(STEP.DRB)
#    RBRes = STEP.RBRes
    println(norm(LSRes - RBRes))
    println(norm(LSRes) - norm(RBRes))
#    TP = pdata.TP
#    DXout = TP * (LSout.xval - STEP.xfin)
#    plot(DXout)
    return (LSout = LSout, RBout=RBout, STEP)
end

function PlotNew(
    m = 1000,
    k = 1;
    larray=[.1;.2;.5;.7;1.0],
    noiselev = 0.0,
    delta = 0.1,
    nonlinear = false,
    mu = 1.e-8,
    regparm=0.0
)
    pdata = SetDEData(m, noiselev, delta, nonlinear, 1.e-8, regparm)
#    dran = collect(0.0:0.025:5.0)
    dran = collect(0.0:0.01:0.5)
    DeltaRec=zeros(length(dran),length(larray))
    pdata.ADel[1] = 0.0
    DLen=length(dran)
    LLen=length(larray)
    LSout = G2(pdata)
    LSRes = LSout.NLRes
    for ilam = 1:LLen
    pdata.ADel[1] = larray[ilam]
    RBout = G2(pdata)
    RBRes = RBout.NLRes
 ( norm(RBRes) < norm(LSRes) ) && error("norm error")
 ( norm(LSRes,1) < norm(RBRes,1) ) && println("L1 norm error")
    for il = 1:DLen
    delta=dran[il]
    STEP = One_Shot2(k, LSRes, RBRes, delta, pdata)
    DeltaRec[il,ilam]=STEP.DRB
println("lambda = ", larray[ilam], "norm = ", STEP.DRB )
    end
    end
    plot(dran,DeltaRec[:,1],"k--")
    plot(dran,DeltaRec[:,2],"k-")
    plot(dran,DeltaRec[:,3],"k-.")
    plot(dran,DeltaRec[:,4],"k_")
    plot(dran,DeltaRec[:,5],"k.")
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



function SetDEData(m, noiselev, delta, nonlinear, mu, regparm; 
                   xinit = [], obj=de_obj)
    LD1 = D1(m)
    G = GmatRL(m + 2)
    AF = G[2:m+1, 2:m+1]
    Z = zeros(m, m)
    AF = I + Z
    dx = 1.0 / (m + 1)
    xg = collect(dx:dx:1.0-dx)
#    xe=ones(m)
    xa=(xg .- .25); xe=abs.(xa);
#    xe = exp.(xg) .* sin.(4.0 * xg)
    if nonlinear
        b = AF * xe
    else
        xnl = nlf.(xe)
        b = AF * xnl
    end
    if noiselev > 1.e-8
#        perturb = ones(m) + noiselev*rand(m)
#        b .*= perturb
         b .+= norm(b,Inf)*noiselev*rand(m)
    end
    TP = chebmat(m, 10)
    TPD = LD1*TP
    A = AF
    linit = length(xinit)
    if linit == 0
        (ma, na) = size(AF)
        (mt, nt) = size(TP)
        (nt == 1) ? (x0 = ones(na) / sqrt(na)) : (x0 = ones(nt) / sqrt(nt))
    else
        x0 = xinit
    end
    robust = (delta > 1.e-8)
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
        regparm = regparm
    )
    return pdata
end

function G2(pdata)
    optout = do_opt(pdata)
    TP = pdata.TP
    xe = pdata.xe
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
    NLRes = Fres(optout.xval, pdata)
    return (
        xval = optout.xval,
        fval = fval,
        gradnorm = gradnorm,
        delx = delx,
        success = success,
        sol = sol,
        NLRes = NLRes,
    )
end


function One_Shot2(k, LSRes, RBRes, delta, pdata)
oldway = false
if oldway
#    RBout = G2(pdata)
    delta = pdata.ADel[1]
    m = pdata.m
#    RBRes = RBout.NLRes
    mres = length(RBRes)
    YA = Gen_Y(mres, k, delta)
    LSErr = zeros(k)
    RBErr = zeros(k)
    for ip = 1:k
        LSErr[ip] = norm(abs.(LSRes) .+ YA[:, ip])^2
        RBErr[ip] = norm(abs.(RBRes) .+ YA[:, ip])^2
    end
    DRA = LSErr - RBErr
    DRB = sum(LSErr .- RBErr) / Float64(k)
    DRL = sum(LSErr) / Float64(k)
    DRR = sum(RBErr) / Float64(k)
else
#LSErr = norm(abs.(LSRes) .+ delta)^2
#RBErr = norm(abs.(RBRes) .+ delta)^2
LSErr = norm(LSRes).^2 + 2.0*delta*norm(LSRes,1)
RBErr = norm(RBRes).^2 + 2.0*delta*norm(RBRes,1)
end
DRB = LSErr - RBErr
DRL = LSErr
DRR = RBErr
#println(DRB,"  ", DRR," ", DRL)
#    RBData = [RBout.fval, RBout.gradnorm, RBout.delx, delta]
    return (
        DRB = DRB,
        DRL = DRL,
        DRR = DRR,
#        RBData = RBData,
#        xfin = RBout.xval,
#        RBRes = RBRes,
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



