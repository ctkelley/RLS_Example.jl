function tables()
println("Table 1")
HistPrint()
println("Table 2")
HistPrint(;nonlinear=false)
println("Table 3")
HistPrint(;nonlinear=false, lasso=[true])
end

function figures()
figure(1)
PlotNew()
title("Nonlinear, left")
figure(2)
PlotNew(; longplot=true)
title("Nonlinear, right")
figure(3)
PlotNew(;nonlinear=false)
title("Linear, left")
figure(4)
PlotNew(; nonlinear=false, longplot=true)
title("Linear, right")
figure(5)
figure(6)
end

