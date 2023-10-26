include("FKM.jl")
using .FKM

Nc = 64
q = 2.5
u = 0.0
gamma = 5.0
U = [unitary(Nc, 2.0 .* (rand(Float64,Nc^2) .- 0.5) .* pi) for i in 1:RANK]

S1 = action(U, Nc, q, u)
S2 = action(U, Nc, q)

tau = 0.1
Nmax = 10
check_force_balance(Nc, gamma, q, tau, Nmax)