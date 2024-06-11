function HistPrint(
    noiselev = 0.0,
    lambda=.5,
    mu = 1.e-8,
    nonlinear = true)
    rownum = 2
    OutVals = zeros(rownum, 3) 
    headers = ["\$\\delta\$", "\$\\Psi\$","\$\\| \\nabla \\Psi \\|\$"]
    formats = "%7.2e & %7.2e & %7.2e" 
delarray = [0.0, lambda]
#for idel=1:2
#delta=delarray[idel]
Optout= T2(1000,1 ;noiselev = noiselev, delta = lambda,
            mu = mu, nonlinear = nonlinear, regparm=0.0)
OutVals[1,1] = 0.0
OutVals[1, 2] = Optout.LSout.fval
OutVals[1, 3] = Optout.LSout.gradnorm
OutVals[2,1] = lambda
OutVals[2, 2] = Optout.RBout.fval
OutVals[2, 3] = Optout.RBout.gradnorm
fprintTeX(headers, formats, OutVals)
OutVals
end

