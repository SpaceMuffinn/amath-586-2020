
using Plots

using Elliptic.Jacobi

β₁ = 0.
β₂ = 1.
β₃ = 10.
c = (β₁ + β₂ + β₃)/3
t = 0.:.01:10
v = t -> β₂ + (β₃ - β₂)*cn(sqrt((β₃-β₁)/12)*t, (β₃-β₂)/(β₃-β₁) )^2 # Julia uses the square of the elliptic modulus

plot(t, map(v,t))

f = u -> [u[2], u[3], u[2]*(c - u[1])]
Df = u -> [0. 1. 0.; 0. 0. 1.; -u[2] c-u[1] 0.]
u₀ = [β₃,0.,-1.0/6*(β₃-β₁)*(β₃-β₂)]

h = 0.0001
[v(0), (v(h)-v(-h))/(2h), (v(h)-2v(0)+v(-h))/(h^2)]

T = 10.# Final time.
k = 0.0001 # Step size

# Forward Euler
n = convert(Int64,ceil(T/k))# Number of time steps, converted to Int64
U = zeros(3,n+1) # To save the solution values
U[:,1] = u₀
t = zeros(n+1) # To save times
t[1] = 0.
for i = 2:n+1
    U[:,i] = U[:,i-1] + k*f(U[:,i-1])
    t[i] = t[i-1] + k
end

plot(t,U[1,:],label="Solution with Forward Euler")
plot!(t,map(v,t),label="True solution")

T = 10. # Final time.
k = .02
p = 7
data = zeros(p)
ks = zeros(p)
for i = 1:p
    k = k/2
    n = convert(Int64,ceil(T/k))
    println("Number of time steps = ", n)
    U = zeros(3,n+1) # To save the solution values
    U[:,1] = u₀
    t = zeros(n+1,1)
    t[1] = 0.
    for i = 2:n+1
        U[:,i] = U[:,i-1] + k*f(U[:,i-1])
        t[i] = t[i-1] + k
    end
    data[i] = abs(U[1,end] - v(t[end]))
    ks[i] = k
end
data_fe = data

ks

plot(ks,data,lw=2,ms=5,marker=:d, minorgrid = true, xaxis=:log, yaxis= :log,label="1st order (Fwd Euler)")

T = 10. # Final time.
k = .02
p = 7
data = zeros(p)
ks = zeros(p)
for i = 1:p
    k = k/2
    n = convert(Int64,ceil(T/k))
    println("Number of time steps = ", n)
    U = zeros(3,n+1) # To save the solution values
    U[:,1] = u₀
    t = zeros(n+1,1)
    t[1] = 0.
    U[:,2] = U[:,1] + k*f(U[:,1])  # Begin the method using
    t[2] = t[1] + k                # forward Euler
    for i = 3:n+1
        U[:,i] = U[:,i-2] + (2*k)*f(U[:,i-1]) #Leapfrog
        t[i] = t[i-1] + k
    end
    data[i] = abs(U[1,end] - v(t[end]))
    ks[i] = k
end
data_leap = data

plot!(ks,data,lw=2,ms=5,marker=:d, minorgrid = true, xaxis=:log, yaxis= :log,label="2nd order (Leapfrog)")

using LinearAlgebra

A = randn(3,3)

I

A + I # The "size" of I is inferred 

Matrix{Float64}(I,2,2) # If you REALLY need to construct the identity matrix

A = randn(10000,10000)
@time I + A;
@time Matrix{Float64}(I,10000,10000) + A;

g = (u,Un) -> u - Un - (k/2)*(f(u)+f(Un))
Dg = u -> I - (k/2)*Df(u)

T = 10 # Final time.
k = 0.02
p = 7
data = zeros(p)
ks = zeros(p)
for i = 1:p
    k = k/2
    n = convert(Int64,ceil(T/k))
    println("Number of time steps = ", n)
    U = zeros(3,n+1) # To save the solution values
    U[:,1] = u₀
    t = zeros(n+1,1)
    t[1] = 0.
    max_iter = 10
    for i = 2:n+1
        t[i] = t[i-1] + k
        Unew = U[:,i-1]
        Uold = U[:,i-1]
        for j = 1:max_iter
            Uold = Unew
            Unew = Uold - (Dg(Uold)\g(Uold,U[:,i-1]))
            if maximum(abs.(Unew-Uold)) < k^2/10 # Newton's method until error tol.
                break                            # Use k^2 because the method is second order
            end
            
            if j == max_iter
                println("Newton didn't terminate")
            end
            
        end
        U[:,i] = Unew
    end
    data[i] = abs(U[1,end] - v(t[end]))
    ks[i] = k
end
data_trap = data

plot!(ks,data,lw=2,ms=5,marker=:d, minorgrid = true, xaxis=:log, yaxis= :log,label="2nd order (Trapezoid)")

g = (u,Un,Um) -> u - Un - (k/12)*(-f(Um)+8*f(Un)+5*f(u))
Dg = u -> I - (5k/12)*Df(u)

T = 10 # Final time.
k = .02
p = 7
data = zeros(p)
ks = zeros(p)
for i = 1:p
    k = k/2
    n = convert(Int64,ceil(T/k))
    println("Number of time steps = ", n)
    U = zeros(3,n+1) # To save the solution values
    U[:,1] = u₀
    t = zeros(n+1,1)
    t[1] = 0.
    max_iter = 10
    # Runge-Kutta second order here
    Us = U[:,1] + (k/2)*f(U[:,1])
    U[:,2] = U[:,1] + k*f(Us)
    t[2] = t[1] + k
    for i = 3:n+1
        t[i] = t[i-1] + k
        Unew = U[:,i-1]
        Uold = U[:,i-1]
        for j = 1:max_iter
            Uold = Unew
            Unew = Uold - (Dg(Uold)\g(Uold,U[:,i-1],U[:,i-2]))
            if maximum(abs.(Unew-Uold)) < k^3/100 # Newton's method until error tol.
                break                            # Use k^3 because the method is third order
            end
            
            if j == max_iter
                println("Newton didn't terminate")
            end
            
        end
        U[:,i] = Unew
    end
    data[i] = abs(U[1,end] - v(t[end]))
    ks[i] = k
end
data_am = data

plot!(ks,data,lw=2,ms=5,marker=:d, minorgrid = true, xaxis=:log, yaxis= :log,label="3rd order (AM)")

T = 10 # Final time.
k = .02
p = 7
data = zeros(p)
ks = zeros(p)
for i = 1:p
    k = k/2
    n = convert(Int64,ceil(T/k))
    println("Number of time steps = ", n)
    U = zeros(3,n+1) # To save the solution values
    U[:,1] = u₀
    t = zeros(n+1,1)
    t[1] = 0.
    max_iter = 10
    for i = 2:n+1
        t[i] = t[i-1] + k
        Y1 = U[:,i-1]
        f1 = f(Y1)    
        Y2 = U[:,i-1] + (k/2)*f1
        f2 = f(Y2)
        Y3 = U[:,i-1] + (k/2)*f2
        f3 = f(Y3)
        Y4 = U[:,i-1] + k*f3
        f4 = f(Y4)
        U[:,i] = U[:,i-1] + (k/6)*(f1+2*f2+2*f3+f4)
    end
    data[i] = abs(U[1,end] - v(t[end]))
    ks[i] = k
end
data_rk = data

1e-11/64000

2.2e-16*64000

plot!(ks,data,lw=2,ms=5,marker=:d, minorgrid = true, xaxis=:log, yaxis= :log,label="4th order until saturation (RK)")

methods = ["Forward Euler","Leapfrog","Trapezoid","Adams-Moulton","Runge-Kutta"];
d = Dict([(methods[1],data_fe),(methods[2],data_leap),(methods[3],data_trap),(methods[4],data_am),(methods[5],data_rk)]) # Use a dictionary, because we can

data_table = Array{Any,2}(undef,7,6);
data_table[1,1] = "k"
data_table[1,2:end] = methods # create column headings
for j = 2:7
    data_table[j,1] = ks[j]
end
l = 1
for i in methods
    l += 1
    for j = 2:7
        data_table[j,l] = d[i][j-1]/d[i][j]
    end
end

using Printf

@printf("%s        | %s |  %s |  %s |  %s  | %s \n",data_table[1,:]...)
for j=2:7
    @printf("%f | %0.4f        |  %0.4f   |  %0.4f    |  %0.4f         | %0.4f \n",data_table[j,:]...)
end
