using DataFrames
using CSV
using Plots
using Plots.PlotMeasures
using Glob

obs = "specificheat"
N=16
maxbin = 20 

begin
  files0 = glob("論文用データ/DT/JN$(obs)_N$(N)*.csv")
  files = sort(files0, by=file -> parse(Int, match(r"g(\d+).csv", basename(file)).captures[1]))
  
  removed_q = []
  
  plt = plot(xlabelfontsize=14, ylabelfontsize=14, size=(1000, 550), margin=20px)
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
    scatter!(plt, df2.s, df2.jmean, yerror=df2.jerr, markersize=4, markerstrokewidth=0.4, alpha=0.8, label="\$\\gamma = $(gamma)\$")
  end
  
g = 16384
R = 0.657298106138376
C1(g,s)  = 2*2*g^2*R^(6*s)
C1r(g,s) = 2*2*g^2*R^(6*(1-s))
f(q) = 1/(4*q^3)-3/(16*q^5)+17/(64*q^7)-63/(256*q^9)
C1d(g,s) = 2*2*g^2*f(R^s)^2
C2(g,s)  = 1.0

sc = log(1/(2*g)^(1/3))/log(R)
sce = -7.15657

x = sc:0.001:10.0
plot!(plt, x, C1.(g,x), linewidth=2, linecolor=:blue2, label="GWW approximation")
x = -10.0:0.001:1-sc
plot!(plt, x, C1r.(g,x), linewidth=2, linecolor=:red2, linestyle=:dot, label="flip of GWW in \$s>1\$")
x = -10.0:0.001:sce
plot!(plt, x, C1d.(g,x), linewidth=2, linecolor=:green2, label="dual approximation")
x = 1.0:0.001:sc
plot!(plt, x, C2.(g,x), linewidth=2, linecolor=:blue2, label=nothing)
x = 1-sc:0.001:0
plot!(plt, x, C2.(g,x), linewidth=2, linecolor=:red2, linestyle=:dot, label=nothing)
x = sce:0.001:0
plot!(plt, x, C2.(g,x), linewidth=2, linecolor=:green2, label=nothing)


  # 臨界帯
plot!(x -> 0.0, 0.0, 1.0, fillrange=x -> 10.0, linealpha=0,
fillcolor=RGB(0.8, 0.8, 1), alpha=0.5, label="critical strip\n(unstable region)")

  xlabel!(plt,"\$s\$")
  ylabel!(plt,"specific heat")
  ylims!(plt,(0, 1.2))
  xlims!(plt,(-10,11))
  plot!(plt, legend=:bottom)
  display(plt)
end
savefig(plt,"./DT-SH.pdf")

removed_q

