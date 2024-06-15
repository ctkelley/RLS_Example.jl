function mini_test(
    m = 1000,
    k = 1;
    noiselev = 0.0,
    delta = 0.1,
    nonlinear = false,
    mu = 1.e-8,
    regparm=.0
)
pdata = SetDEData(m, noiselev, delta, nonlinear, mu, regparm)
return(pdata)
end


