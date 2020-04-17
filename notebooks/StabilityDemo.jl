
using Plots, Printf

h = t -> sin(t)^2
dh = t -> 2*sin(t)*cos(t)
f = (u,t) -> λ*(u-h(t))+dh(t)
Df = u -> λ

T = 10.# Final time.
t = 0:.01:T
plot(t,map(h,t),label=@sprintf("Attractor"))

# Forward Euler
λ = -2;

k = 0.001 # Step size
u₀ = 2.
n = convert(Int64,ceil(T/k))# Number of time steps, converted to Int64
U = zeros(n+1) # To save the solution values
U[1] = u₀
t = zeros(n+1) # To save times
t[1] = 0.
for i = 2:n+1
    U[i] = U[i-1] + k*f(U[i-1],t[i-1])
    t[i] = t[i-1] + k
end
plot!(t,U,label=@sprintf("Fwd Euler, k = %0.4f, lam = %3.1f",k,λ))

# Forward Euler
λ = -20;
T = 10.# Final time.
k = 0.001 # Step size
u₀ = 2.
n = convert(Int64,ceil(T/k))# Number of time steps, converted to Int64
U = zeros(n+1) # To save the solution values
U[1] = u₀
t = zeros(n+1) # To save times
t[1] = 0.
for i = 2:n+1
    U[i] = U[i-1] + k*f(U[i-1],t[i-1])
    t[i] = t[i-1] + k
end
plot!(t,U,label=@sprintf("Fwd Euler, k = %0.4f, lam = %3.1f",k,λ))

# Forward Euler
λ = -2001;
T = 10.# Final time.
k = 0.001 # Step size
u₀ = 2.
n = convert(Int64,ceil(T/k))# Number of time steps, converted to Int64
U = zeros(n+1) # To save the solution values
U[1] = u₀
t = zeros(n+1) # To save times
t[1] = 0.
for i = 2:n+1
    U[i] = U[i-1] + k*f(U[i-1],t[i-1])
    t[i] = t[i-1] + k
end
plot!(t,U,label=@sprintf("Fwd Euler, k = %0.4f, lam = %3.1f",k,λ))

T = 10.# Final time.
t = 0:.01:T
plot(t,map(h,t),label=@sprintf("Attractor"));

λ = -4.;
T = 4.# Final time.
k = 0.0001 # Step size
u₀ = 2.
n = convert(Int64,ceil(T/k))# Number of time steps, converted to Int64
U = zeros(n+1) # To save the solution values
U[1] = u₀
t = zeros(n+1)
t[1] = 0.
U[2] = U[1] + k*f(U[1],t[1])  # Begin the method using
t[2] = t[1] + k                # forward Euler
for i = 3:n+1
    U[i] = U[i-2] + (2*k)*f(U[i-1],t[i-1]) #Leapfrog
    t[i] = t[i-1] + k
end
plot!(t,U,label=@sprintf("Leapfrog, k = %0.4f, lam = %3.1f, T = %f",k,λ,T))

g = (u,Un,t,tn) -> u - Un - (k/2)*(f(u,t)+f(Un,tn))
Dg = u -> 1 - (k/2)*Df(u)

T = 10.# Final time.
t = 0:.01:T
plot(t,map(h,t),label=@sprintf("Attractor"));

λ = -4.;
T = 10.# Final time.
k = 0.001 # Step size
u₀ = 2.
n = convert(Int64,ceil(T/k))
U = zeros(n+1) # To save the solution values
U[1] = u₀
t = zeros(n+1)
t[1] = 0.
max_iter = 10
for i = 2:n+1
    t[i] = t[i-1] + k
    Unew = U[i-1]
    Uold = U[i-1]
    for j = 1:max_iter
        Uold = Unew      
        Unew = Uold - (Dg(Uold)\g(Uold,U[i-1],t[i],t[i-1]))  # Linear equation, will converge in 1 step, always
        if maximum(abs.(Unew-Uold)) < k^2/10 # Newton's method until error tol.
            break                            # Use k^2 because the method is second order
        end
            
        if j == max_iter
            println("Newton didn't terminate")
        end
            
        end
    U[i] = Unew
end
plot!(t,U,label=@sprintf("Trapezoid, k = %0.4f, lam = %3.1f",k,λ))

λ = -40.;
T = 10.# Final time.
k = 0.001 # Step size
u₀ = 2.
n = convert(Int64,ceil(T/k))
U = zeros(n+1) # To save the solution values
U[1] = u₀
t = zeros(n+1)
t[1] = 0.
max_iter = 10
for i = 2:n+1
    t[i] = t[i-1] + k
    Unew = U[i-1]
    Uold = U[i-1]
    for j = 1:max_iter
        Uold = Unew      
        Unew = Uold - (Dg(Uold)\g(Uold,U[i-1],t[i],t[i-1]))  # Linear equation, will converge in 1 step, always
        if maximum(abs.(Unew-Uold)) < k^2/10 # Newton's method until error tol.
            break                            # Use k^2 because the method is second order
        end
            
        if j == max_iter
            println("Newton didn't terminate")
        end
            
        end
    U[i] = Unew
end
plot!(t,U,label=@sprintf("Trapezoid, k = %0.4f, lam = %3.1f",k,λ))

λ = -2001.;
T = 10.# Final time.
k = 0.001 # Step size
u₀ = 2.
n = convert(Int64,ceil(T/k))
U = zeros(n+1) # To save the solution values
U[1] = u₀
t = zeros(n+1)
t[1] = 0.
max_iter = 10
for i = 2:n+1
    t[i] = t[i-1] + k
    Unew = U[i-1]
    Uold = U[i-1]
    for j = 1:max_iter
        Uold = Unew      
        Unew = Uold - (Dg(Uold)\g(Uold,U[i-1],t[i],t[i-1]))  # Linear equation, will converge in 1 step, always
        if maximum(abs.(Unew-Uold)) < k^2/10 # Newton's method until error tol.
            break                            # Use k^2 because the method is second order
        end
            
        if j == max_iter
            println("Newton didn't terminate")
        end
            
        end
    U[i] = Unew
end
plot!(t,U,label=@sprintf("Trapezoid, k = %0.4f, lam = %3.1f",k,λ))

λ = -40000.;
T = 10.# Final time.
k = 0.001 # Step size
u₀ = 2.
n = convert(Int64,ceil(T/k))
U = zeros(n+1) # To save the solution values
U[1] = u₀
t = zeros(n+1)
t[1] = 0.
max_iter = 10
for i = 2:n+1
    t[i] = t[i-1] + k
    Unew = U[i-1]
    Uold = U[i-1]
    for j = 1:max_iter
        Uold = Unew      
        Unew = Uold - (Dg(Uold)\g(Uold,U[i-1],t[i],t[i-1]))  # Linear equation, will converge in 1 step, always
        if maximum(abs.(Unew-Uold)) < k^2/10 # Newton's method until error tol.
            break                            # Use k^2 because the method is second order
        end
            
        if j == max_iter
            println("Newton didn't terminate")
        end
            
        end
    U[i] = Unew
end
plot!(t,U,label=@sprintf("Trapezoid, k = %0.4f, lam = %3.1f",k,λ))
