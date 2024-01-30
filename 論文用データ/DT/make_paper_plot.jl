using DataFrames
using CSV
using Plots

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
  xlabel!(plt,"\$s\$")
  ylabel!(plt,"specific heat")
  xlims!(plt,(-15,16))
  display(plt)
end
savefig(plt,"DT_dC.png")

removed_q

