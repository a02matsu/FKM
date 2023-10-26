module Measurement
using CSV
using DataFrames
using Plots
using Glob
using LsqFit
using Printf
using LsqFit
using Statistics
using DelimitedFiles
include("Graph.jl")
using .Graph

export calc_obs, findPT, plot_gamma, plot_N, analysePT, plot_gamma_small_q, sof, qof

################################
# qとsの対応
################################
# (q,u)からsへの変換
sof(q,u=0.0) = -log(q)/log((1-u)*(OMEGA-u))
# (s,u)からqへの変換
qof(s,u=0.0) = ( (1-u)*(OMEGA-u) )^(-s) 

################################
# カラーリストを定義
ColorList = [:blue, :red, :green, :cyan, :purple, :yellow]


###################################
# 物理量のプロット
# 同じパラメータのファイルを集めて、統計処理をしてcsvに描き込む
function calc_obs(phys,Nc, gamma_int)
    pattern = "Obs/$(phys)_N$(Nc)g$(gamma_int)q*u00*.txt"
    file_paths = glob(pattern)
    
    results = Dict{Float64, Vector{Float64}}()
    
    for file_path in file_paths
        q_int = parse(Int, match(r"q(\d+)", file_path).captures[1])
        data = vec(readdlm(file_path, Float64))
        
        if haskey(results, q_int/10000)
            append!(results[q_int/10000], data)
        else
            results[q_int/10000] = data
        end
    end
    
    summary = DataFrame(q=Float64[], mean_val=Float64[], stderr_val=Float64[])
    
    for (q, values) in results
        mean_val = mean(values)
        stderr_val = std(values) / sqrt(length(values))
        
        push!(summary, (q, mean_val, stderr_val))
    end
    
    sort!(summary, :q)
    
    CSV.write("$(phys)_N$(Nc)g$(gamma_int)u00.csv", summary)
end


###################################
# 相転移点を見つける
function findPT(Nc, gamma_int, epsilon=10^(-4),epsilon2=10^(-5))
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
  Qcb = []
  for a in 1:RANK
    q1, q2 = 0.0, 0.0
    q1b, q2b = 0.0, 0.0
    crossed_below, crossed_above = false, false
    crossed_below2, crossed_above2 = false, false
    for i in 1:length(A[a])
      if !crossed_below && A[a][i] < epsilon
          q1 = Q[i]
          crossed_below = true
      end
      if !crossed_below2 && A[a][i] < epsilon2
          q1b = Q[i]
          crossed_below2 = true
      end
      if crossed_below2 && !crossed_above2 && A[a][i] > epsilon2
          q2b = Q[i-1]
          crossed_above2 = true
      end
      if crossed_below2 && !crossed_above && A[a][i] > epsilon
          q2 = Q[i-1]
          crossed_above = true
      end
    end
    push!(Qc,[q1,q2])
    push!(Qcb,[q1b,q2b])
  end
  return(Qc,Qcb)
end


###################################
# gammaを固定してプロット
function plot_gamma(phys, gamma_int)
    # CSVファイルを見つける
    files_org = glob("$(phys)_N*g$(gamma_int)u00.csv")
    # Nを抽出し、ファイルを昇順にソート
    files = sort(files_org, by=file -> parse(Int, match(r"_N(\d+)g", basename(file)).captures[1]))
    # プロットを作成
    plt = plot()
    # 色のインデックスを追跡
    cindex = 1
    Nmax = 0
    for file in files
        # CSVファイルからデータを読み込む
        df = DataFrame(CSV.File(file))
        # ファイル名からgammaを取得
        gamma_str = match(r"g(\d+)u00", file).captures[1]
        gamma_int = parse(Int, gamma_str)
        gamma = gamma_int / 100
        Nc = match(r"_N(\d+)", file).captures[1]
        # 最後のNcが保存されるようにする
        Nmax = Nc
        # 横軸、縦軸、エラーバーのデータを計算
        #q,mean_val,stderr_val
        x = sof.(df.q,0.0) # 横軸は log(q)
        y = df.mean_val
        err = df.stderr_val
        # プロットを追加
        scatter!(plt, x, y, yerror=err, label="\$N_c\$ = $(Nc)", markershape=:circle, markersize=3, markercolor=ColorList[cindex])
        cindex += 1
    end
    # 安定領域の境界
    sb1 = sof(0.70222,0.0)
    sb2 = sof(1.0,0.0)
    vline!(plt, [sb1,sb2], linestyle=:dash, linecolor=:black, label="")
    # 相転移の位置
    # 一番大きいNcのデータから読み取る
    Qc, Qcb = findPT(Nmax,gamma_int)
    for a in 1:RANK
      vline!(plt, sof.((Qc[a]+Qcb[a])/2.0,0.0), linestyle=:dash, linecolor=ColorList[a], label="GWW $(a)")
    end
    # 相転移点を出力
    for a in 1:RANK
      println("PT",a,": q=",(Qc[a]+Qcb[a])/2.0," s=",sof.((Qc[a]+Qcb[a])/2.0,0.0))
    end
    # タイトルを設定
    xlabel!("\$s\$")
    ylabel!("$(phys)")
    title!("$(NAME), $(phys), \$ \\gamma = $(gamma_int/100) \$")
    #if phys == "specificheat"
      #ylims!(0,2.2)
    #end
    #if phys == "dC"
      #ylims!(-1,7)
    #end
    # プロットを表示
    display(plt)
    # プロットをファイルに保存
    savefig(plt, "$(phys)_g$(gamma_int).png")
end


###################################################
# Nを固定して異なるgammaごとにプロット
function plot_N(phys, Nc)
  # CSVファイルを見つける
  files_org = glob("$(phys)_N$(Nc)g*u00.csv")
  # gamma_intを抽出し、ファイルを昇順にソート
  files = sort(files_org, by=file -> parse(Int, match(r"g(\d+)u00", basename(file)).captures[1]))
  # プロットを作成
  plt = plot()
  for file in files
      # CSVファイルからデータを読み込む
      df = DataFrame(CSV.File(file))
      # ファイル名からgammaを取得
      gamma_str = match(r"g(\d+)u00", file).captures[1]
      gamma_int = parse(Int, gamma_str)
      gamma = gamma_int / 100
      # 横軸、縦軸、エラーバーのデータを計算
      #q_int,mean_val,stderr_val
      x = sof.(df.q,0.0)
      y = df.mean_val
      err = df.stderr_val
      # プロットを追加
      scatter!(plt, x, y, yerror=err, label="\$\\gamma\$ = $(gamma)", markershape=:circle, markersize=3)
  end
  # 安定領域の境界
  sb1 = sof(1/OMEGA,0.0)
  sb2 = sof(1.0,0.0)
  vline!(plt, [sb1,sb2], linestyle=:dash, linecolor=:black, label="")
  # タイトルを設定
  xlabel!("\$s\$")
  ylabel!("$(phys)")
  title!("$(NAME), $(phys), \$ N_c = $(Nc) \$")
  #if phys == "specificheat"
    #ylims!(0,2.2)
  #end
  #if phys == "dC"
    #ylims!(-1,7)
    #hline!(plt, [0.0,1.0,2.0], linestyle=:dash, linecolor=:black, label="")
  #end
  # プロットを表示
  display(plt)
  # プロットをファイルに保存
  savefig(plt, "$(phys)_N$(Nc).png")
end


#################################################
## 相転移の位置を解析
function analysePT(Qmax=16, epsilon=10^(-4), epsilon2=10^(-5))
  gammas = [400,800, 1600,12800,102400]
  Qcs = [Float64[] for _ in 1:RANK]
  dQcs = [Float64[] for _ in 1:RANK]
  Qcsb = [Float64[] for _ in 1:RANK]
  dQcsb = [Float64[] for _ in 1:RANK]
  for g in gammas
    Qc, Qcb = findPT(Qmax,g,epsilon,epsilon2)
    for a in 1:RANK
      push!(Qcs[a],Qc[a][1])
      push!(dQcs[a],Qc[a][2])
      push!(Qcsb[a],Qcb[a][1])
      push!(dQcsb[a],Qcb[a][2])
    end
  end

  LG=log.(gammas/100 .- 0.5) # log(γ-1/2)
  ## fitting
  # モデル関数を定義
  model(lg, p) = p[1]*lg .+ p[2]
  # 初期パラメーター値を設定
  # フィッティングを実行
  # (1) q<1
  p1=[]
  erra1=[]
  errb1=[]
  #weights = [1.0/wt for _ in 1:length(LG)]
  weights_Sc = []
  weights_dSc = []
  for a in 1:RANK
    diffSc = abs.(sof.(Qcs[a])-sof.(Qcsb[a]))
    diffdSc = abs.(sof.(dQcs[a])-sof.(dQcsb[a]))
    push!(weights_Sc, 1.0 ./ diffSc)
    push!(weights_dSc, 1.0 ./ diffdSc)
  end
  for a in 1:RANK 
    S1=sof.(Qcs[a],0.0)
    p0 = [-0.2, 0.0]
    fit1 = curve_fit(model, LG, S1, weights_Sc[a], p0)
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
    fit1 = curve_fit(model, LG, S1, weights_dSc[a], p0)
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
    scatter!(pltq, LG, yerror=1.0./weights_Sc[a], sof.(Qcs[a],0.0) ,label="")
    scatter!(pltq, LG, yerror=1.0./weights_dSc[a], sof.(dQcs[a],0.0),label="" )
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
  savefig(pltq,"sc_vs_gamma.png")

  ## 結果の出力
  # Qcの値
  println("<position of PT>")
  for a in 1:RANK
    println("Cycle $(a)")
    Qc = [@sprintf("%.4f", val) for val in (Qcs[a])]
    dQc = [@sprintf("%.4f", val) for val in (dQcs[a])]
    Qcb = [@sprintf("%.4f", val) for val in (Qcsb[a])]
    dQcb = [@sprintf("%.4f", val) for val in (dQcsb[a])]
    Sc = [@sprintf("%.4f", val) for val in sof.((Qcs[a] .+ Qcsb[a])./2.0)]
    dSc = [@sprintf("%.4f", val) for val in sof.((dQcs[a] .+ dQcsb[a])./2.0)]
    Sum = [@sprintf("%.4f", val) for val in sof.((dQcs[a] .+ dQcsb[a])/2.0) + sof.((Qcs[a] .+ Qcsb[a]) ./2.0 )] 
    println("\t qc_bd1 = : [", join(Qc,","),"]" )  
    println("\t qc_bd2 = : [", join(Qcb,","),"]" )  
    println("\t dqc_bd1 = : [", join(dQc,","),"]" )
    println("\t dqc_bd2 = : [", join(dQcb,","),"]" )
    println("\t sc = : [", join(Sc,","),"]" )  
    println("\t dsc = : [", join(dSc,","),"]" )
    println("\t sum : [", join(Sum,","),"]")
  end
  for a in 1:RANK
    println("== Cycle $(a) ==")
    # 対称性
    println("   <symmetry>")
    errasum = sqrt(erra1[a]^2+erra2[a]^2)
    errbsum = sqrt(errb1[a]^2+errb2[a]^2)
    println("\t (+) : (",@sprintf("%.4f",p1[a][1]), " +/- ", @sprintf("%.4f",erra1[a]),") log(γ-1/2) + (", @sprintf("%.4f",p1[a][2]), " +/- ", @sprintf("%.4f",errb1[a]), ")")
    println("\t (-) : (",@sprintf("%.4f",p2[a][1]), " +/- ", @sprintf("%.4f",erra2[a]),") log(γ-1/2) + (", @sprintf("%.4f",p2[a][2]), " +/- ", @sprintf("%.4f",errb2[a]), ")" )
    println("\t sum : (",@sprintf("%.4f",p1[a][1]+p2[a][1]), " +/- ", @sprintf("%.4f",errasum),")  log(γ-1/2) + (", @sprintf("%.4f",p1[a][2]+p2[a][2]), " +/- ", @sprintf("%.4f",errbsum),")" )
  # GWWとの比較
    println("   <GWWとの比較>")
    fac = OMEGA^(-poly[a]*p1[a][2])
    expn = -poly[a]*p1[a][1]*log(OMEGA)
    fac_err = abs( -poly[a]*log(OMEGA)*OMEGA^(-poly[a]*p1[a][2]) ) * errb1[a]
    expn_err = poly[a]*log(OMEGA)*erra1[a]
    println("\t n=$(poly[a])(+): q_c^{$(poly[a])} = ","(",(@sprintf("%.3f", fac))," +/- ",@sprintf("%.4f",fac_err),")", "(γ-1/2)^(", @sprintf("%.3f",expn), " +/- ", @sprintf("%.4f",expn_err),")")
    #fac = exp(-n*p1[a][2])
    #expn = -n*p1[a][1]
    #fac_err = abs( -n*exp(-n*p1[a][2]) ) * errb1[a]
    #expn_err = n*erra1[a]
    #println("n=$(n)(+): q_c^{$(n)} = ","(",(@sprintf("%.3f", fac))," +/- ",@sprintf("%.4f",fac_err),")", "(γ-1/2)^(", @sprintf("%.3f",expn), " +/- ", @sprintf("%.4f",expn_err),")")
    ###
    fac = OMEGA^(poly[a]*(p2[a][2]-1))
    expn = poly[a]*p2[a][1]*log(OMEGA)
    fac_err = abs( poly[a]*log(OMEGA)*OMEGA^(poly[a]*(p2[a][2])-1) ) * errb2[a]
    expn_err = poly[a]*log(OMEGA)*erra2[a]
    println("\t n=$(poly[a])(-): q'_c^{-$(poly[a])}/Omega^$(poly[a]) = ","(",(@sprintf("%.3f", fac))," +/- ",@sprintf("%.4f",fac_err),")", "(γ-1/2)^(", @sprintf("%.3f",expn), " +/- ", @sprintf("%.4f",expn_err),")")
    #fac = exp(n*(p2[a][2]))
    #expn = n*p2[a][1]
    #fac_err = abs( n*exp(n*(p2[a][2])) ) * errb2[a]
    #expn_err = n*erra2[a]
    #println("n=$(n)(-): q'_c^{-$(n)} = ","(",(@sprintf("%.3f", fac))," +/- ",@sprintf("%.7f",fac_err),")", "(γ-1/2)^(", @sprintf("%.3f",expn), " +/- ", @sprintf("%.7f",expn_err),")")
  end
end
  
  


############################################################
# small q 領域だけのプロットを作成してCnと比較
function plot_gamma_small_q(phys, gamma_int)
  # Cnの比熱理論値
  C1(g,a) = -2*g^2*log(1 - a^2)
  C2(g,a) = -2*g^2*(2*log(2) + log(g*(g-1)/(2*g-1)^2))
  # gamma
  gamma = gamma_int / 100
  # CSVファイルを見つける
  files_org = glob("$(phys)_N*g$(gamma_int)u00.csv")
  # Ncを抽出し、ファイルを昇順にソート
  files = sort(files_org, by=file -> parse(Int, match(r"_N(\d+)", basename(file)).captures[1]))

  # プロットを作成
  plt = plot()
  # 色のインデックスを追跡
  cindex = 1
  Nmax = 0
  for file in files
      # CSVファイルからデータを読み込む
      df = DataFrame(CSV.File(file))

      # ファイル名からNcを取得
      N = parse(Int,match(r"_N(\d+)", file).captures[1])
      if N > Nmax
        Nmax = N
      end
      # 横軸、縦軸、エラーバーのデータを計算
      #q_int,mean_val,stderr_val
      x = sof.(df.q,0.0)
      y = df.mean_val
      err = df.stderr_val
      # プロットを追加
      #plot!(x, y, yerror=err, label="\$\\gamma\$ = $(gamma)", markershape=:circle, markersize=4)
      scatter!(plt, x, y, yerror=err, label="\$N\$ = $(N)", markershape=:circle, markersize=3, linecolor=ColorList[cindex], markercolor=ColorList[cindex])
      # color indexを更新
      cindex += 1
  end
  # 相転移の位置
  # 一番大きいNcのデータから読み取る
  Qc,Qcb = findPT(Nmax,gamma_int)
  # small qの領域のみプロット
  for a in 1:RANK
    vline!(plt, [sof(Qc[a][1],0.0)], linestyle=:dash, linecolor=:black, label=nothing)
  end

  # タイトルを設定
  xlabel!("\$s\$")
  ylabel!("$(phys)")
  title!("$(phys), \$ \\gamma = $(gamma_int/100),\$ small \$q\$")
  minS=1.0
  maxS=sof(Qc[1][1],0.0)+2.0
  xlims!(minS,maxS)
  if phys == "specificheat"
    S1 = sof(Qc[2][1],0.0)
    S2 = sof(Qc[1][1],0.0)
    y1 = minS:0.01:S1
    plot!(plt,y1,2*C2.(gamma,qof.(y1,0.0).^3), linestyle=:dash, linecolor=:black, label=nothing)
    y2 = S1:0.01:S2
    plot!(plt, y2,C2.(gamma,qof.(y2,0.0).^3) .+ C1.(gamma,qof.(y2,0.0).^4), linestyle=:dash, linecolor=:black, label=nothing)
    y3 = S2:0.01:maxS
    plot!(plt, y3,C1.(gamma,qof.(y3,0.0).^3) .+ C1.(gamma,qof.(y3,0.0).^(4)), linestyle=:dash, linecolor=:black,label=nothing)
  end
  # プロットを表示
  display(plt)
  # プロットをファイルに保存
  savefig(plt, "$(phys)_g$(gamma_int)_small_q.png")
end

end