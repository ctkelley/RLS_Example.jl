function Conv_Test(m=1000; noiselev=0.0, nonlinear=true, mu=1.e-8, delta=0.0,
                   regparm=0.0)
pdata = SetDEData(m, noiselev, delta, nonlinear, mu, regparm)
delarr=collect(0.0:.5:10.0)
Dlen=length(delarr)
xvals=zeros(10,Dlen)
fvals=zeros(Dlen)
xdiff=zeros(Dlen-1)
Gout=G2(pdata)
xvals[:,1] .= Gout.xval
fvals[1] = Gout.fval
for il=2:Dlen
     pdata.ADel[1] = delarr[il]
     Gout=G2(pdata)
     xvals[:,il] .= Gout.xval
     fvals[il] = Gout.fval
     dif = norm(xvals[:,il] - xvals[:,il-1])
     xdiff[il-1] = dif
end
figure(1)
semilogy(xdiff)
figure(2)
semilogy(fvals)
end

