using DataFrames
using CSV
using Plots
using Plots.PlotMeasures
using Glob

obs = "specificheat"
N=16
maxbin = 20 

begin
  files0 = glob("論文用データ/TT/JN$(obs)_N$(N)*.csv")
  files = sort(files0, by=file -> parse(Int, match(r"g(\d+).csv", basename(file)).captures[1]))
  
  removed_q = []
  
  plt = plot(xlabelfontsize=14, ylabelfontsize=14, size=(1000, 700), margin=20px)
  for file in files
    df = DataFrame(CSV.File(file))
    # ファイル名からgammaを取得
    gamma = match(r"g(\d+).csv", file).captures[1]
    # binの値がmaxbin未満のデータのみをフィルタリング
    df2 = filter(row -> row[:bin] < maxbin, df)
    # binの値がmaxbinを超えるqの値を収集
    tmp = df[df.bin .> maxbin, :q] 
    pushfirst!(tmp,parse(Int,gamma))
    push!(removed_q,tmp )
    scatter!(plt, df2.s, df2.jmean, yerror=df2.jerr, markersize=2, alpha=0.8, label="\$\\gamma = $(gamma)\$")
  end
  
g = 16384
R = 0.53656

#C1(g,s)  = 2*2*g^2*R^(6*s)
#C1r(g,s) = 2*2*g^2*R^(6*(1-s))
#f(q) = 1/(4*q^3)-3/(16*q^5)+17/(64*q^7)-63/(256*q^9)
#C1d(g,s) = 2*2*g^2*f(R^s)^2

sc = log(1/(2*g)^(1/3))/log(R)

x = 1.0:0.001:sc
C(g,s)  = 1.5
plot!(plt, x, C.(g,x), linewidth=2, linecolor=:blue, label=nothing)

x = sc:0.001:7
C(g,s)  = 3*2*g^2*R^(6*s)
plot!(plt, x, C.(g,x), linewidth=2, linecolor=:blue, label=nothing)

x = 1-sc:0.001:0
C(g,s)  = 1.5
plot!(plt, x, C.(g,x), linewidth=2, linecolor=:red, linestyle=:dot, label="simple reflection")

x = -8:0.001:1-sc
C(g,s)  = 3*2*g^2*R^(6*(1-s))
plot!(plt, x, C.(g,x), linewidth=2, linecolor=:red, linestyle=:dot, label=nothing)

# Dual approximation
g1(q) = 1/(6*q^3)-11/(72*q^5)+233/(864*q^7)
g2(q) = 1/(12*q^3)-17/(144*q^5)+115/(576*q^7)

sce1 = -4.605927568158268
sce2 = -4.2325020832132205

x = sce2:0.001:0
C(g,s)  = 1.5
plot!(plt, x, C.(g,x), linewidth=2, linecolor=:blue, label="dual approximation")

x = sce1:0.001:sce2
C(g,s)  = 1.0 + 2*g^2*g2(R^s)^2
plot!(plt, x, C.(g,x), linewidth=2, linecolor=:blue, label=nothing)

x = -6:0.001:sce1
C(g,s)  = 2*2*g^2*g1(R^s)^2 + 2*g^2*g2(R^s)^2
plot!(plt, x, C.(g,x), linewidth=2, linecolor=:blue, label=nothing)

  # 臨界帯
plot!(x -> 0.0, 0.0, 1.0, fillrange=x -> 10.0, linealpha=0,
fillcolor=RGB(0.8, 0.8, 1), alpha=0.5, label="critical strip\n(unstable region)")

  xlabel!(plt,"\$s\$")
  ylabel!(plt,"specific heat")
  ylims!(plt,(0.0, 2.0))
  xlims!(plt,(-6,7))
  plot!(plt, legend=:bottom)
  display(plt)
end
savefig(plt,"./TT-SH.pdf")

removed_q

