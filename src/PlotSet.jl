function PlotSet(lassoflag, nonlinear)
if nonlinear
    larray=[0.0, .1, 1.0, 5.0, 10.0, 100.0]
    longplot ?  dran = collect(0.0:0.025:5.0) : dran = collect(0.0:.1:1.0)
end
if ~nonlinear
    larray=[0.0, .5, 1.0, 10.0, 100.0, 200.0]
    longplot ?  dran = collect(0.0:1.0:50.0) : dran = collect(0.0:1.0:20.0)
    if lassoflag
    larray=[0.0, .5, 1.0, 10.0, 100.0, 200.0]
    longplot ?  dran = collect(0.0:1.0:50.0) : dran = collect(0.0:.2:2.0)
    end
end
return (larray, dran)
end

function Lasso_Base(pdata, lassoflag)
lassoflag && (pdata.lasso[1]=true)
LSout=G2(pdata)
LSRes = LSout.NLRes.FV
m=length(LSRes); mx=Int(m/2)
useC ? (CLSRes = LSRes[1:mx]) : (CLSRes = LSRes)
return (LSRes, CLSRes)
end

