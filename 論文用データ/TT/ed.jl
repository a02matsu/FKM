using CSV
using DataFrames
using Plots
using Plots.PlotMeasures


#file = "./論文用データ/TT/Phases/phases_N16g131072q3522E1u00.csv" # s=-5.72039
#file = "./論文用データ/TT/Phases/phases_N16g131072q2780E1u00.csv" # s=-5.3487
file = "./論文用データ/TT/Phases/phases_N16g131072q1563E-2u00.csv" # s=6.68011


df = CSV.read(file, DataFrame)

cycle1 = df."Cycle 1" * 16
cycle2 = df."Cycle 2" * 16
cycle3 = df."Cycle 3" * 16

plt = plot(xlabelfontsize=14, ylabelfontsize=14, legendfontsize=14, size=(500, 350), margin=20px)

xlims!(-3.14159, 3.14159)
ylims!(0, 0.5)
xlabel!("\$\\theta\$")
ylabel!("\$\\rho(\\theta)\$")

x = range(-3.14159, 3.14159, length=100)
bar!(plt,x,cycle1,alpha=0.5,label="\$U_1\$")
bar!(plt,x,cycle2,alpha=0.5,label="\$U_2\$")
bar!(plt,x,cycle3,alpha=0.5,label="\$U_3\$")

savefig(plt, "TT-ED3.pdf")