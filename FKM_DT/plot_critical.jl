using CSV
using DataFrames
using Plots
using Glob
using LsqFit
using Printf
using LsqFit
using Statistics
using DelimitedFiles
using LinearAlgebra
using JLD
using StatsBase

F0(g,q) = -g^2*log(1/((2*q^2 - 1)*(q^2 + 1)^2*(q^2 - 1)^3*(2*q^4 + q^2 + 1)^3))/2
F1(g,q) = -g*log(1/((-1 + 2*q)*(q + 1)^2*(q - 1)^3*(2*q^2 + q + 1)^3)) + (3*log(g))/2 + log(32*q^9/((-1 + 2*q)^3*(2*q^2 + q + 1)^6*(q - 1)^3))/2 + 9/4

E0(g,q) = -g^2*log(1/((2*q^2 - 1)*(q^2 + 1)^2*(q^2 - 1)^3*(2*q^4 + q^2 + 1)^3))
E1(g,q) = -g*log(1/((-1 + 2*q)*(q + 1)^2*(q - 1)^3*(2*q^2 + q + 1)^3)) + 3/2

C0(g,q) = g^2*log(1/((2*q^2 - 1)*(q^2 + 1)^2*(q^2 - 1)^3*(2*q^4 + q^2 + 1)^3))
C1(g,q) = 3/2

## energy
g0 = 2^20
plt = plot()
file_d=glob("Obs_down/energy_N16g1048576u00.csv")
df = DataFrame(CSV.File(file_d))
x = df.q
y = df.mean_val .* g0
eerr = df.stderr_val .* g0
scatter!(plt, x, y, yerror=eerr, label="large to small", markershape=:circle)

file_u=glob("Obs_up/energy_N16g1048576u00.csv")
df = DataFrame(CSV.File(file_u))
x = df.q
y = df.mean_val .* g0
eerr = df.stderr_val .* g0
scatter!(plt, x, y, yerror=eerr, label="small to large")

plot!(plt, x, E0.(2^20,x), label="\$E_0\$")
plot!(plt, x, E1.(2^20,x), label="\$E_1\$")

savefig(plt,"energy_hystereses.png")

##########################
plt = plot()
file_d=glob("Obs_down/specificheat_N16g1048576u00.csv")
df = DataFrame(CSV.File(file_d))
x = df.q
y = df.mean_val
eerr = df.stderr_val
scatter!(plt, x, y, yerror=eerr, label="large to small", markershape=:circle)

file_u=glob("Obs_up/specificheat_N16g1048576u00.csv")
df = DataFrame(CSV.File(file_u))
x = df.q
y = df.mean_val
eerr = df.stderr_val
scatter!(plt, x, y, yerror=eerr, label="small to large")

plot!(plt, x, C0.(2^20,x), label="\$C_0\$")
plot!(plt, x, C1.(2^20,x), label="\$C_1\$",ylims=(0,2.5))

savefig(plt,"specificheat_hystereses.png")