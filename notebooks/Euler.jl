
function test_fun(A)
    A[1,1] = 2*A[1,1]
    return A
end

A = [1 2 3; 4 5 6] # Integer array

test_fun(A)

A # A has changed

A = [1 2 3; 4 5 6] # Integer array
function test_fun2(A)
    A = 2*A
    return A
end
test_fun2(A)

A # A has not changed

A = [1 2 3; 4 5 6] # Integer array
function test_fun3(A)
    A[:,1] = 2*A[:,1]
    return A
end
test_fun2(A)

A # A has not changed

α = 1.

f = u -> [u[2], -u[3]*u[1]+2*u[1]^3, 1.0] # use commas to get a vector in Julia

T = 40.# Final time.
k = 0.0001 # Step size

# Forward Euler
n = convert(Int64,T/k)# Number of time steps, converted to Int64
U = zeros(3,n+1) # To save the solution values
U[:,1] = [1.,-1.,0.]
for i = 2:n+1
    U[:,i] = U[:,i-1] + k*f(U[:,i-1])
end

using Plots # Import plotting functionality

plot(U[3,:],U[1,:],label="Solution with Forward Euler")

g = (U,Un) -> U - Un - k*f(U)
Dg = (U) -> [1. -k 0.0; 
    k*U[3]-6*k*U[1]^2 1 k*U[1];
    0.0 0.0 1.0 ]

# Backward Euler
n = convert(Int64,T/k) # Number of time steps, converted to Int64
U = zeros(3,n+1) # To save the solution values
U[:,1] = [1.,-1.,0.]
max_iter = 10
for i = 2:n+1
    Unew = U[:,i-1]
    Uold = U[:,i-1]
    for j = 1:max_iter
        Uold = Unew
        Unew = Uold - (Dg(Uold)\g(Uold,U[:,i-1]))
        #println(maximum(abs.(Unew-Uold)))
        if maximum(abs.(Unew-Uold)) < k/10 # Newton's method until error tol.
            break
        end
        if j == max_iter
            println("Newton didn't terminate")
        end
    end
    U[:,i] = Unew
end

plot!(U[3,:],U[1,:],label="Solution with Backward Euler") # similar to "hold on"

f = u -> [u[2], u[3]*u[1]+2*u[1]^3, 1.0] 
T = 6.# Final time.
k = 0.0001 # Step size
n = convert(Int64,T/k) # Number of time steps, converted to Int64
U = zeros(3,n+1) 
U[:,1] = [0.3670615515480782,-0.2953721054475503,0.]
for i = 2:n+1
    U[:,i] = U[:,i-1] + k*f(U[:,i-1])
end

plot(U[3,:],U[1,:],label="Solution with Forward Euler, starts to blow up at t = 6")

f = u -> [u[2], u[3]*u[1]+2*u[1]^3, 1.0] 
g = (U,Un) -> U - Un - k*f(U)
Dg = (U) -> [1. -k 0.0; 
    -k*U[3]-6*k*U[1]^2 1 -k*U[1];
    0.0 0.0 1.0 ]

T = 6.# Final time.
k = 0.0001 # Step size
n = convert(Int64,T/k) # Number of time steps, converted to Int64
U = zeros(3,n+1) 
U[:,1] = [0.3670615515480782,-0.2953721054475503,0.]
max_iter = 10
for i = 2:n+1
    Unew = U[:,i-1]
    Uold = U[:,i-1]
    for j = 1:max_iter
        Uold = Unew
        Unew = Uold - (Dg(Uold)\g(Uold,U[:,i-1]))
        #println(maximum(abs.(Unew-Uold)))
        if maximum(abs.(Unew-Uold)) < k/10 # Newton's method until error tol.
            break
        end
        if j == max_iter
            println("Newton didn't terminate")
        end
    end
    U[:,i] = Unew
end

plot!(U[3,:],U[1,:],label="Solution with Backward Euler, starts to blow up in the other direction!",title="Solution should tend to zero exponentially")
