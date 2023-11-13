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

F0(g,q) = -g^2*log(-1/(4*q^20 - q^16 - 4*q^14 - 4*q^12 + 2*q^8 + 4*q^6 - 1))/2
F1(g,q) = -g*log(-1/((q + 1)*(q^2 + 1)*(2*q^2 + q + 1)*(2*q^3 + q^2 - 1)*(q - 1)^2)) + log(g) + log(4*q^6*(q + 1)^2/((2*q^2 + q + 1)^2*(q - 1)^2*(q^2 + 1)*(2*q^3 + q^2 - 1)^2))/2 + 3/2

E0(g,q) = -g^2*log(-1/(4*q^20 - q^16 - 4*q^14 - 4*q^12 + 2*q^8 + 4*q^6 - 1))
E1(g,q) = -g*log(-1/(4*q^10 - q^8 - 4*q^7 - 4*q^6 + 2*q^4 + 4*q^3 - 1)) + 1

C0(g,q) = g^2*log(-1/(4*q^20 - q^16 - 4*q^14 - 4*q^12 + 2*q^8 + 4*q^6 - 1))
C1(g,q) = 1

## energy
g0 = 2^10
plt = plot()
file=glob("energy_N16g102400u00.csv")
df = DataFrame(CSV.File(file))
x = df.q
y = df.mean_val .* g0
eerr = df.stderr_val .* g0
scatter!(plt, x, y, yerror=eerr, label="energy", xlims=(0.05,0.1),ylims=(-10.5,5))

x = range(0.0,0.6,1000)
plot!(plt, x, E0.(g0,x), label="\$E_0\$")
plot!(plt, x, E1.(g0,x), label="\$E_1\$")

xlabel!(plt,"\$q\$")
title!(plt,"DT, \$\\gamma=1024\$, \$N_c=16\$")
savefig(plt,"energy_critical.png")

##########################
plt = plot()
file=glob("specificheat_N16g102400u00.csv")
df = DataFrame(CSV.File(file))
x = df.q
y = df.mean_val
eerr = df.stderr_val
scatter!(plt, x, y, yerror=eerr, label="specificheat", xlims=(0.05,0.1),ylims=(0,1.5))

x = range(0.0,0.6,1000)
plot!(plt, x, C0.(g0,x), label="\$C_0\$")
plot!(plt, x, C1.(g0,x), label="\$C_1\$")

xlabel!(plt,"\$q\$")
title!(plt,"DT, \$\\gamma=1024\$, \$N_c=16\$")
savefig(plt,"specificheat_critical.png")