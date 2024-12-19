function HistPrint(m=1000;
    noiselev = 0.0,
    lambda=.5,
    mu = 1.e-8,
    regparm=.1,
    nonlinear = true,
    useC = true,
    lasso=[false])
#    larray=[0.0, 0.0005, 0.001, 0.01, .1, 1.0]
if nonlinear
    larray=[0.0, .1, 1.0, 5.0, 10.0, 100.0]
else
    larray=[0.0, 1.0, .5, 10.0, 100.0, 200.0]
end
    rownum = length(larray)
    OutVals = zeros(rownum,4)
    headers = ["\$\\lambda\$", "\$\\Psi_\\lambda\$", 
            "\$\\| \\nabla \\Psi_\\lambda \\|\$","slope"]
    formats = "%7.2e & %7.2e & %7.2e & %7.2e"
    pdata = SetDEData(m, noiselev, 0.0, nonlinear, mu, regparm; 
              useC = useC, lasso=lasso)
    if lasso[1] 
    qdata = SetDEData(m, noiselev, 0.0, nonlinear, mu, regparm; 
              useC = useC, lasso=[false])
    end
       nrm10=0.0
    for ilam=1:rownum
       pdata.ADel[1] = larray[ilam] 
#       optout = do_opt(pdata)
       optout = G2(pdata)
       nrmlam = norm(optout.NLRes.FV,1)
       useC && (nrmlam = norm(optout.NLRes.GV,1))
       if lasso[1]
       nrm10=nrmlam
       qdata.ADel[1] = larray[ilam] 
       optout = G2(qdata)
       nrmlam = norm(optout.NLRes.FV,1)
       useC && (nrmlam = norm(optout.NLRes.GV,1))
       else
       (ilam == 1) && (nrm10 = nrmlam)
       end
       slope = 2.0*(nrm10 - nrmlam)
       OutVals[ilam,1] =larray[ilam]
       OutVals[ilam,2] = optout.fval
       OutVals[ilam,3] = optout.gradnorm
       OutVals[ilam,4] = slope
    end
fprintTeX(headers, formats, OutVals)
end




