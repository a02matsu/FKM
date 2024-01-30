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

  # 臨界帯
plot!(x -> 0.0, 0.0, 1.0, fillrange=x -> 10.0, linealpha=0,
fillcolor=RGB(0.8, 0.8, 1), alpha=0.5, label="critical strip\n(unstable region)")

  xlabel!(plt,"\$s\$")
  ylabel!(plt,"specific heat")
  ylims!(plt,(0, 1.5))
  xlims!(plt,(-11,12))
  plot!(plt, legend=:top)
  display(plt)
end
savefig(plt,"./DT-SH.pdf")

removed_q

