using CSV
using DataFrames
using Plots
using Glob
using LsqFit
using Printf
using LsqFit

include("Graph.jl")
using .Graph
################################
# (q,u)からsへの変換
sof(q,u=0.0) = -log(q)/log((1-u)*(OMEGA-u))
# (s,u)からqへの変換
qof(s,u=0.0) = ( (1-u)*(OMEGA-u) )^(-s) 

function findPT(Nc, gamma_int, epsilon=2*10^(-3))
  # CSVファイルを見つける
  files_org = glob("Phases/phases_N$(Nc)g$(gamma_int)*u00.csv")
  # qを抽出し、ファイルを昇順にソート
  files = sort(files_org, by=file -> parse(Int, match(r"q(\d+)u", basename(file)).captures[1]))
  # RANK を読み取る
  df = DataFrame(CSV.File(files[1]))
  Ndata, RANK = size(df) 
  # 
  Q=[]
  A = [[] for _ in 1:RANK]
  for file in files
    # qの値を読み取ってQに格納
    qobj = match(r"q(\d+)", file)
    qint = parse(Int, qobj.captures[1])
    push!(Q, qint/10000)
    # CSVファイルからデータを読み込む
    df = DataFrame(CSV.File(file))
    # 各列について、最初と最後の行の値を加え、対応するAのサブ配列に追加
    for a in 1:RANK
      a_value = df[1, a] + df[end, a]
      push!(A[a], a_value)
    end
  end
  # 相転移点を見つける
  #e = 10^(-4)
  Qc = []
  for a in 1:RANK
    q1, q2 = 0.0, 0.0
    crossed_below, crossed_above = false, false
    for i in 1:length(A[a])
      if !crossed_below && A[a][i] < epsilon
          q1 = Q[i]
          crossed_below = true
      end
      if crossed_below && !crossed_above && A[a][i] > epsilon
          q2 = Q[i-1]
          crossed_above = true
      end
    end
    push!(Qc,[q1,q2])
  end
  return(Qc)
end
##################
## 相転移の位置をプロット
gammas = [400,800, 1600]
Qcs = [[] for _ in 1:RANK]
dQcs = [[] for _ in 1:RANK]
for g in gammas
  Qc = findPT(16,g,10^(-4))
  #Qc = findPT(16,g,10^(-4))
  for a in 1:RANK
    push!(Qcs[a],Qc[a][1])
    push!(dQcs[a],Qc[a][2])
  end
end

LG=log.(gammas/100 .- 0.5) # log(γ-1/2)
#LG=log.(gammas/100 ) # log(γ-1/2)
## fitting
# モデル関数を定義
model(lg, p) = p[1]*lg .+ p[2]
# 初期パラメーター値を設定
# フィッティングを実行
# (1) q<1
p1=[]
erra1=[]
errb1=[]
weights = [1.0/0.05 for _ in 1:length(LG)]
for a in 1:RANK 
  S1=sof.(Qcs[a],0.0)
  p0 = [-0.2, 0.0]
  fit1 = curve_fit(model, LG, S1, weights, p0)
  # 結果を表示
  push!(p1, fit1.param)
  # 誤差共分散行列
  covar = estimate_covar(fit1)
  # 傾きの標準誤差
  push!(erra1, sqrt.(covar[1,1]))
  push!(errb1, sqrt.(covar[2,2]))
end
# (2) q>1
p2=[]
erra2=[]
errb2=[]
for a in 1:RANK 
  S1=sof.(dQcs[a],0.0)
  p0 = [0.2, 0.0]
  fit1 = curve_fit(model, LG, S1, weights, p0)
  # 結果を表示
  push!(p2, fit1.param)
  # 誤差共分散行列
  covar = estimate_covar(fit1)
  # 傾きの標準誤差
  push!(erra2, sqrt.(covar[1,1]))
  push!(errb2, sqrt.(covar[2,2]))
end
# プロット作成
pltq = plot()
for a in 1:RANK
  scatter!(pltq, LG, yerror=1.0./weights, sof.(Qcs[a],0.0) ,label="")
  scatter!(pltq, LG, yerror=1.0./weights, sof.(dQcs[a],0.0),label="" )
end
for a in 1:RANK
  F1(x) = model(x,p1[a])
  F2(x) = model(x,p2[a])
  plot!(pltq,LG,F1,linestyle=:dash,linecolor=:black,label="")
  plot!(pltq,LG,F2,linestyle=:dash,linecolor=:black,label="")
end
xlabel!("\$\\log (\\gamma - 1/2)\$")
ylabel!("\$s_c\$")
display(pltq)
savefig(pltq,"sc_vs_gamma_small.png")
println("----------")
println(Qcs)
println(dQcs)
for a in 1:RANK
  errasum = sqrt(erra1[a]^2+erra2[a]^2)
  errbsum = sqrt(errb1[a]^2+errb2[a]^2)
  println("----------")
  println("cycle $(a)(+) : (",@sprintf("%.4f",p1[a][1]), " +/- ", @sprintf("%.4f",erra1[a]),") log(γ-1/2) + (", @sprintf("%.4f",p1[a][2]), " +/- ", @sprintf("%.4f",errb1[a]), ")")
  println("cycle $(a)(-) : (",@sprintf("%.4f",p2[a][1]), " +/- ", @sprintf("%.7f",erra2[a]),") log(γ-1/2) + (", @sprintf("%.4f",p2[a][2]), " +/- ", @sprintf("%.7f",errb2[a]), ")" )
  println("sum $(a) : ",@sprintf("%.4f",p1[a][1]+p2[a][1]), " +/- ", @sprintf("%.4f",errasum),"  log(γ-1/2) + ", @sprintf("%.4f",p1[a][2]+p2[a][2]), " +/- ", @sprintf("%.4f",errbsum) )
  println("----------")
  println("q^3チェック")
  n = 3
  fac = OMEGA^(-n*p1[a][2])
  expn = -n*p1[a][1]*log(OMEGA)
  fac_err = abs( -n*log(OMEGA)*OMEGA^(-n*p1[a][2]) ) * errb1[a]
  expn_err = n*log(OMEGA)*erra1[a]
  println("n=$(n)(+): q_c^{$(n)} = ","(",(@sprintf("%.3f", fac))," +/- ",@sprintf("%.4f",fac_err),")", "(γ-1/2)^(", @sprintf("%.3f",expn), " +/- ", @sprintf("%.4f",expn_err),")")
  ##
  fac = OMEGA^(n*(p2[a][2]-1))
  expn = n*p2[a][1]*log(OMEGA)
  fac_err = abs( n*log(OMEGA)*OMEGA^(n*(p2[a][2])-1) ) * errb2[a]
  expn_err = n*log(OMEGA)*erra2[a]
  println("n=$(n)(-): q'_c^{-$(n)}/Omega^$n = ","(",(@sprintf("%.3f", fac))," +/- ",@sprintf("%.7f",fac_err),")", "(γ-1/2)^(", @sprintf("%.3f",expn), " +/- ", @sprintf("%.7f",expn_err),")")
end
