using DataFrames
using CSV
using Plots

obs = "specificheat"
N=16
maxbin = 150

files0 = glob("Obs/JN$(obs)_N$(N)*.csv")
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
xlabel!(plt,"\$s\$")
ylabel!(plt,"specific heat")
xlims!(plt,(-14,15))
display(plt)
savefig(plt,"TS_SH.png")

removed_q

