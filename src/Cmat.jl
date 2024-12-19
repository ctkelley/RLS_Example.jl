function Cmat(m)
mx=Int(m/2)
C=spdiagm(m, mx, ones(mx))
return C
end

function eye(m)
C=spdiagm(m,m,ones(m))
return C
end
