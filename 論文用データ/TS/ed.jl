using CSV
using DataFrames
using Plots
using Plots.PlotMeasures


#file = "./論文用データ/TS/Phases/phases_N16g16384q3137E-2u00.csv"
file = "./論文用データ/TS/Phases/phases_N16g16384q7433E-2u00.csv"

df = CSV.read(file, DataFrame)

cycle1 = df."Cycle 1" * 16
cycle2 = df."Cycle 2" * 16

plt = plot(xlabelfontsize=14, ylabelfontsize=14, legendfontsize=12, size=(500, 350), margin=20px)

xlims!(-3.14159, 3.14159)
ylims!(0, 0.5)
xlabel!("\$\\theta\$")
ylabel!("\$\\rho(\\theta)\$")

x = range(-3.14159, 3.14159, length=100)
bar!(plt,x,cycle1,alpha=0.5,label="\$U_1\$")
bar!(plt,x,cycle2,alpha=0.5,label="\$U_2\$")

savefig(plt, "TS-ED1.pdf")