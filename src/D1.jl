function D1(n)
dx=1.0/(n+1)
dv=.5/dx
dl=-dv*ones(n-1)
du=dv*ones(n-1)
d=zeros(n)
D1=Tridiagonal(dl,d,du)
D1[1,1]=-1.0/dx
D1[1,2]=1.0/dx
D1[n,n]=1.0/dx
D1[n,n-1]=-1.0/dx
return D1
end

function testd1(n,f,fp)
L=D1(n)
dx=1.0/(n+1)
x=collect(dx:dx:1.0-dx)
fv=f.(x)
fpv=fp.(x)
dv=L*fv
println(norm(dv-fpv)/sqrt(n))
end

