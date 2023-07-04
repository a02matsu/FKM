include("Measurement.jl")
using .Measurement
using StatsBase
using CSV
using DataFrames
using Plots
using Glob
using LsqFit
using Printf
using LsqFit
using Statistics
using DelimitedFiles
using LinearAlgebra
using JLD

######################################################
## 与えられたUconfとcycleに対応するWilson loopのヒストグラムを返す
function hist_Wilson(Uconf,Nc,cycle)
  phases = [] # phaseを格納するリスト
  for U in Uconf
    tmp = Array{ComplexF64}(I, Nc, Nc)
    for a in cycle
      if a <= RANK
        tmp = copy(tmp) * U[a]
      else
        tmp = copy(tmp) * adjoint(U[a-RANK])
      end
    end
    push!(phases, angle.(eigvals(tmp)))
  end
  # phaseのヒストグラムデータを書き出す
  hist = []
  # ビンのエッジを定義
  bin_edges = range(-3.1416, stop=3.1416, length=101) # length=101 to create 100 bins
  phase = reduce(vcat,phases)
  h = fit(Histogram, phase, bin_edges)
  return h.weights ./ length(phase) 
end

file = "Uconfig/Uconf_N4g1024q1145E0u00.jld"
Uconf = load(file, "Uconf")
h = hist_Wilson(Uconf,4,[1,2])
println(h)
m = match(r"q(\d+E-?[0-9]+)", file)
m.captures[1]
cycles = [[1,2,3],[1],[3],[4],[2],[2,8],[1,2],[2,3],[1,2,3,1,2,3],[1,3],[1,4],[3,4],[1,4,3]]
  Cnames = []
for i in eachindex(cycles)
  push!(Cnames,"C$(join(cycles[i],""))")
end
line1 = join(Cnames,",")

######################################################
## 与えられたNcとgammaに対応するUconfを読み取り、
## cycleのリストに対応するWilson loopの固有値のphaseのヒストグラムを
## CSVファイル(Cphase_N**g**q**u00.csv)にまとめる
function PhaseHist(Nc::Int, gamma::Int, cycles)
  # それぞれのサイクルに名前を付ける
  Cnames = []
  for i in eachindex(cycles)
    push!(Cnames,"C$(join(cycles[i],""))")
  end
  line1 = join(Cnames,",")

  # 該当するUconfファイルを見つける
  files_org = glob("Uconfig/Uconf_N$(Nc)g$(gamma)*u00.jld")
  # qを抽出し、ファイルを昇順にソート
  files = sort(files_org, by=file -> read_q(file))
  # ファイルごとにCSVを作成
  count = 0
  for file in files
    # qの値を読み取る
    q = read_q(file)
    m = match(r"q(\d+E-?[0-9]+)", file)
    qval = m.captures[1]
    # configurationをload
    Uconf = load(file, "Uconf")
    # cycleごとにUを読み取り、ヒストグラムを生成
    hists_tmp = []
    for cycle in cycles
      push!(hists_tmp, hist_Wilson(Uconf, Nc, cycle))
    end
    # 行と列を入れ替える
    hists = [[row[i] for row in hists_tmp] for i in eachindex(hists_tmp[1])]
    # CSVファイルを準備
    Hfile = "Phases/Cphases_N$(Nc)g$(gamma)q$(qval)u00.csv"
    open(Hfile, "w") do f
      write(f, line1, "\n") # 予め準備したindexを書き込み
      for hist in hists
        write(f, join(hist,","), "\n")
      end
    end
    count += 1
    println("$count / $(length(files)) $(Hfile)")
  end
end


######################################################
## PhaseHistで作ったCSVファイルからヒストグラムを読み取り、相転移の位置を読み取る
## 具体的には、-πからπまでを100分割したヒストグラムの最初と最後の度数の和が
## epsilon2以上、epsilon以下になる領域を相転移領域と定め、
## その領域のqの平均値を相転移の位置、その幅の半分をerrとする。
## 返り値は、
## Qc1,err1 : q<1領域での相転移位置とその誤差
## Qc2,err2 : q<1領域での相転移位置とその誤差
## Q : 全体のqの集合
function findPT2(Nc::Int,gamma::Int,epsilon=10^-3,epsilon2=10^-4)
  # 該当するUconfファイルを見つける
  files_org = glob("Phases/Cphases_N$(Nc)g$(gamma)*u00.csv")
  # qを抽出し、ファイルを昇順にソート
  files = sort(files_org, by=file -> read_q(file))
  # cycleの数を読み取る
  df = DataFrame(CSV.File(files[1]))
  _, Ncycles = size(df) 
  Q=[] # qの値
  A = [[] for _ in 1:Ncycles]
  for file in files
    # qの値を読み取ってQに格納
    push!(Q, read_q(file))
    # CSVファイルからデータを読み込む
    df = CSV.read(file, DataFrame) 
    # 各列について、最初と最後の行の値を加え、対応するAのサブ配列に追加
    for a in 1:Ncycles
      a_value = df[1, a] + df[end, a]
      push!(A[a], a_value)
    end
  end
  # 相転移点を見つける
  #e = 10^(-4)
  Qc1 = [] # q<1 での相転移点
  Qc2 = [] # q>1 での相転移点
  err1 = [] # q<1 での相転移点の誤差
  err2 = [] # q>1 での相転移点の誤差
  for a in 1:Ncycles
    q1, q2 = 0.0, 0.0
    q1b, q2b = 0.0, 0.0
    crossed_below, crossed_above = false, false
    crossed_below2, crossed_above2 = false, false
    for i in 1:length(A[a])
      if !crossed_below && A[a][i] < epsilon
          q1 = Q[i]
          crossed_below = true
      end
      if crossed_below && !crossed_below2 && A[a][i] < epsilon2
          q1b = Q[i]
          crossed_below2 = true
      end
      if crossed_below2 && !crossed_above2 && A[a][i] > epsilon2 && Q[i]>1.0
          q2b = Q[i-1]
          crossed_above2 = true
      end
      if crossed_above2 && !crossed_above && A[a][i] > epsilon && Q[i]>1.0
          q2 = Q[i-1]
          crossed_above = true
      end
    end
    push!(Qc1,(q1+q1b)/2)
    push!(err1,abs((q1-q1b))/2)
    push!(Qc2,(q2+q2b)/2)
    push!(err2,abs((q2-q2b))/2)
  end
  return Qc1, err1, Qc2, err2, Q
end


cycles = [[1,2,3],[1],[3],[4],[2],[2,8],[1,2],[2,3],[1,2,3,1,2,3],[1,3],[1,4],[3,4],[1,4,3]]
PhaseHist(4,1024,cycles)
findPT2(4,1024)

##########################
## 調べたいサイクル
cycles = [[1,2,3],[1],[3],[4],[2],[2,8],[1,2],[2,3],[1,2,3,1,2,3],[1,3],[1,4],[3,4],[1,4,3]]
Cnames = []
for i in eachindex(cycles)
  push!(Cnames,"C$(join(cycles[i],""))")
end

for gamma in [1024,16384]
  PhaseHist(16,gamma,cycles)
end

Qc1,err1,Qc2,err2,Q = findPT2(16,1024)
sc1=sof.(Qc1)
sc2=sof.(Qc2)
err_sc1 = abs.(sof.(Qc1 .+ err1) - sof.(Qc1))
err_sc2 = abs.(sof.(Qc2 .+ err2) - sof.(Qc2))
sc1 .+ sc2

plt = plot()
scatter!(plt, Cnames, sc1, yerror=err_sc1, label="" )
scatter!(plt, Cnames, sc2, yerror=err_sc2, label="" )
ylabel!("\$ s \$")
savefig(plt,"C_vs_sc.png")


function plot_critical(Nc,gamma,cycle, cycles, Qc1, Qc2,Q)
  ind = findfirst(x->x==cycle, cycles)
  target1 = Qc1[ind]
  target2 = Qc2[ind]
  _, ind1 = findmin(abs.(Q .- target1)) 
  _, ind2 = findmin(abs.(Q .- target2)) 
  q1 = Q[ind1]
  q2 = Q[ind2]
  exp1 = floor(Int, log10(abs(q1))) #-leading_zeros(q)
  val1 = Int(round(q1*10^(-exp1+3)))
  q1txt = "$(val1)E$(exp1)"
  exp2 = floor(Int, log10(abs(q2))) #-leading_zeros(q)
  val2 = Int(round(q2*10^(-exp2+3)))
  q2txt = "$(val2)E$(exp2)"

  x = range(-3.1416, stop=3.1416, length=100)
  Cname = "cycle$(join(cycle,"_"))"
  plt = plot()

  ###
  file = "Phases/Cphases_N$(Nc)g$(gamma)q$(q1txt)u00.csv"
  df = CSV.read(file, DataFrame)
  y = df[:, Cname]
  bar!(plt, x, y, label="q=$q1txt",alpha=0.5)
  ###
  file = "Phases/Cphases_N$(Nc)g$(gamma)q$(q2txt)u00.csv"
  df = CSV.read(file, DataFrame)
  y = df[:, Cname]
  bar!(plt, x, y, label="q=$q2txt",alpha=0.5)
  ylabel!("s")
  title!("\$ U_{$(join(cycle,""))} \$")
end

plot_critical(16,1024,[2,8],cycles,Qc1,Qc2,Q)




