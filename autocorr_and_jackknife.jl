include("Measurement.jl")
using .Measurement
using Glob
using DelimitedFiles
using DataFrames
using CSV
using Plots

obs = "dC"
N=16
gamma=16384

begin 
  bins = range(1,stop=2000,step=10)
  epsilon = 10^-3

  pattern="Obs/$(obs)_N$(N)g$(gamma)q*.txt"

  files = glob(pattern)
  Q = []
  B = []
  autocorr = []
  jmean = []
  jerr = []

  for file in files 
    data = vec(readdlm(file, Float64))
    # ファイル名からNcを取得
    q = read_q(file)
    push!(Q,q) 
    info = 1
    for bin in bins 
      tmp = autocorrelation(data,bin)
      if tmp < epsilon 
        push!(B, bin)
        push!(autocorr, tmp)
        jmean0, jerr0 = jackknife(data,bin)
        push!(jmean, jmean0)
        push!(jerr, jerr0)
        println(q,"\t",bin,"\t",tmp,"\t",jmean0,"\t",jerr0)
        info = 0
        break
      end
    end
    if info == 1 
      bin = 2000
      tmp = autocorrelation(data,bin)
      push!(B, bin)
      push!(autocorr, tmp)
      jmean0, jerr0 = jackknife(data,bin)
      push!(jmean, jmean0)
      push!(jerr, jerr0)
      println("!!OVER 2000!!", q,"\t",bin,"\t",tmp,"\t",jmean0,"\t",jerr0)
    end
  end
  # Qのソート順を取得
  perm = sortperm(Q)
  # 他の配列も同じ順序でソート
  Q_sorted = Q[perm]
  bin_sorted = B[perm]
  autocorr_sorted = autocorr[perm]
  jmean_sorted = jmean[perm]
  jerr_sorted = jerr[perm]
  
  # データフレームの作成
  df = DataFrame(q=Q_sorted, s=sof.(Q_sorted),  bin=bin_sorted, autocorr=autocorr_sorted, jmean=jmean_sorted, jerr=jerr_sorted)
  
  # CSVファイルの名前
  csvfile = "Obs/JN$(obs)_N$(N)g$(gamma).csv"
  # CSVファイルに保存
  CSV.write(csvfile, df)
end



scatter(sof.(Q_sorted),jmean_sorted,yerror=jerr_sorted)