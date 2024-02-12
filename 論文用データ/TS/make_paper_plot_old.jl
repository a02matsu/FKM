using DataFrames
using CSV
using Plots
using Glob

# 相転移点
sc128 = [-3.88, -2.88, 5.22876020154383, 3.921570151157872]
sc1024 = [-5.87, -4.39, 7.1895452771227655, 5.392158957842074]
sc16384 = [-8.49, -6.37, 9.80392537789468, 7.35294403342101]

obs = "dC"
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
  vline!(plt,sc16384,linestyle=:dash,linecolor=:black,label="")
  vline!(plt,sc1024,linestyle=:dash,linecolor=:black,label="")
  vline!(plt,sc128,linestyle=:dash,linecolor=:black,label="")
  xlabel!(plt,"\$s\$")
  ylabel!(plt,"specific heat")
  xlims!(plt,(-14,15))
  display(plt)
end
savefig(plt,"TS_dC.png")



removed_q

