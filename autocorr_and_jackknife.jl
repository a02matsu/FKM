using Distributed
NPROCS = 18 
addprocs(NPROCS)
nworkers()
@everywhere include("Measurement.jl")
@everywhere using .Measurement
@everywhere using Glob
@everywhere using DelimitedFiles
@everywhere using DataFrames
@everywhere using CSV
using Plots

#@everywhere begin 
#  obs = "specificheat"
#  N=16
#  gamma=16384
#
#  bins = range(1,stop=2000,step=1)
#  epsilon = 10^-3
#
#  pattern="Obs/$(obs)_N$(N)g$(gamma)q*.txt"
#
#  Q = []
#  B = []
#  autocorr = []
#  jmean = []
#  jerr = []
#end

obs = "specificheat"
N=16
gamma=16384

pattern="Obs/$(obs)_N$(N)g$(gamma)q*.txt"
files = glob(pattern)

# 測定ファイルから自己相関とjack knife dataを読み出す関数
@everywhere function find_jack_data(file)
  bins = range(1,stop=2000,step=1)
  epsilon = 10^-3
  data = vec(readdlm(file, Float64))
  q = read_q(file)
  info = 1
  B, autocorr, jmean, jerr = 0, 0.0, 0.0, 0.0
  for bin in bins
    tmp = autocorrelation(data,bin)
    if tmp < epsilon 
      B = bin
      autocorr = tmp
      jmean, jerr = jackknife(data,bin)
      println(q,"\t",bin,"\t",tmp,"\t",jmean,"\t",jerr)
      info = 0
      break
    end
  end
  if info == 1 
    bin = 2000
    tmp = autocorrelation(data,bin)
    B, autocorr = bin, tmp
    jmean, jerr = jackknife(data,bin)
    println("!!OVER 2000!!", q,"\t",bin,"\t",tmp,"\t",jmean,"\t",jerr)
  end
  return q, B, autocorr, jmean, jerr
end

###################################
# ファイルの処理
# pmap により、各workerでjackknifeが実行され、
# (q,B,autocorr,jmean,jerr)のタプルが返される。
# resultにはこのタプルの配列が格納される
results = pmap(find_jack_data, files)
# 結果の分解
Q, B, autocorr, jmean, jerr = [], [], [], [], []
for tuple in results
    push!(Q, tuple[1])
    push!(B, tuple[2])
    push!(autocorr, tuple[3])
    push!(jmean, tuple[4])
    push!(jerr, tuple[5])
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
csvfile = "JNData/JN$(obs)_N$(N)g$(gamma).csv"
# CSVファイルに保存
CSV.write(csvfile, df)



scatter(sof.(Q_sorted),jmean_sorted,yerror=jerr_sorted)