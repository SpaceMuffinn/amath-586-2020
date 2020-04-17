
using Plots, LinearAlgebra

## For fourth-order explicit Runge--Kutta
R = (z) -> 1 + z + z^2/2 + z^3/6 + z^4/24

xrange = [-4,4]; yrange = [-4,4]
contourf(xrange[1]:0.01(1+rand()/10):xrange[2],yrange[1]:0.01(1+rand()/10):yrange[2],(x,y)-> sign((abs(R(x+1im*y)/exp(x+1im*y)))-1),colorbar=false)

## For a fifth-order Taylor series method
R = (z) -> 1 + z + z^2/2 + z^3/6 + z^4/24 + z^5/(24*5)

xrange = [-4,4]; yrange = [-4,4]
contourf(xrange[1]:0.01(1+rand()/10):xrange[2],yrange[1]:0.01(1+rand()/10):yrange[2],(x,y)-> sign((abs(R(x+1im*y)/exp(x+1im*y)))-1),colorbar=false)
