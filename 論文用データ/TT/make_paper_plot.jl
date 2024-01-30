using DataFrames
using CSV
using Plots
using Glob

# 相転移点
sc128 = [2.968984391, -1.969255141, -1.548091682]
sc1024 = [4.082353538, -3.112939419, -2.727361185]
sc16384 =  [5.566845732, -4.605923237, -4.232502082]
sc131072 = [6.680214877, -5.720480015, -5.348782084]

obs = "specificheat"
N=16
maxbin = 20

begin 
  files0 = glob("JNData/JN$(obs)_N$(N)*.csv")
  files = sort(files0, by=file -> parse(Int, match(r"g(\d+).csv", basename(file)).captures[1]))

  removed_q = []

  plt = plot()
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
  vline!(plt,sc131072,linestyle=:dash,linecolor=:black,label="")
  vline!(plt,sc16384,linestyle=:dash,linecolor=:black,label="")
  vline!(plt,sc1024,linestyle=:dash,linecolor=:black,label="")
  vline!(plt,sc128,linestyle=:dash,linecolor=:black,label="")
  xlabel!(plt,"\$s\$")
  ylabel!(plt,obs)
  xlims!(plt,(-8,12))
  display(plt)
end

savefig(plt,"TT_"*obs*".png")



removed_q

