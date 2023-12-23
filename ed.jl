using CSV
using DataFrames
using Plots
using Plots.PlotMeasures


file = "./FKM_Cn/Phases/phases_N16g102400q788u00.csv"

df = CSV.read(file, DataFrame)

data = df."Cycle 1"

plt = plot(xlabelfontsize=14, ylabelfontsize=14, size=(1000, 700), margin=20px)

xlims!(-3.14159, 3.14159)
ylims!(0, 0.04)
xlabel!("\$\\theta\$")
ylabel!("\$\\rho(\\theta)\$")

x = range(-3.14159, 3.14159, length=100)
bar!(plt,x,data,legend=nothing)

savefig(plt, "ed.pdf")