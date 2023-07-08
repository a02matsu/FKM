## change log
# ver.01 （2023/05/30）
# ・dCの出力にγ^3をかけることにした
# ver.02　（2023/05/31）
# ・グラフの名前をGRAPH --> NAMEに変更。（module Graphの住民）
# ・熱化のiteration数を10000から20000に増やした
# ver.03 （2023/06/01）
# ・熱化のiteration数を40000にして、Ntauを決めたらcold startすることにした。
# ・dCのhistoryをplotしておくことにした。
# ver.04 （2023/06/04）
# ・最初の1000個のデータを自動的に削除するようにした
# ver.05 （2023/06/05）
# ・sofとqofをmodule Graphに統合した
# ・Ntauを決める前にstep_sizeを少しずつ大きくして多少熱化することにした
# ver.06 （2023/06/13）
# ・module化した

module Simulation

include("FKM.jl")
using .FKM
using Plots
using Statistics
using Random
using StatsBase
using CSV
using DataFrames
using Glob
using DelimitedFiles
using Printf
using LinearAlgebra
using JLD
using Distributed

export NV, NE, RANK, minl, OMEGA, NAME, e_free, S, T, poly, NL
# 外から見える関数
export 
Metropolis, 
HMC,
backup_file,
check_B,
check_gauge_inv, 
Phases, 
check_force_balance,
LFreverse, 
action, 
unitary,
sof,
qof,
init_rng,
Q_regular,
Q_fine,
start_simulation,
chunk

#########################################
# 乱数を初期化する関数。毎回違うseedを使うためにtime()を使う。
# 並列化するときは、
# id = myid()
# とするとよい。
function init_rng(id=1)
  seed = round(Int, time()) % 6599 * id
  return MersenneTwister(seed)
end

# realの観測量を書き出す
function write_realvalues(filename, RV)
  open(filename, "w") do f
      for val in RV
          write(f, join(real(val), ", "))
          write(f,"\n")
      end
  end
end

#########################################
# qの値を有効数字4桁の整数部分と桁に分離する
function format_q(q, significant_figures=4)
  expo = floor(Int, log10(abs(q)))
  #val = Int(round(q / 10.0^(expo), digits=significant_figures - 1) * 10^(significant_figures-1))
  val = Int( round(q * 10.0^(-expo+significant_figures-1)) )
  return val, expo
end

#########################################
## HMCでのシミュレーションを実行
function start_simulation(Nc, gamma, q, niter=200000, step_size=0.10)
  ## Control Simulation
  MaxNtau = 30 # limit of Ntau
  MinAcc = 0.75 # minimal acceptance ratio
  ## theory parameters
  u = 0e-1
  gamma_int = Int(round(gamma))
  # gammaを100倍して整数にし、ファイル名に利用
  #gamma_int = Int(round(gamma*100))
  #Nf = gamma * Nc 
  # シミュレーションのパラメータ
  SSint = Int(round(step_size*100)) 


  #Nc = 8   
  #gamma = 16.0
  #niter = 200000
  #step_size = 0.10
  if !isdir("Conf")
      mkdir("Conf")
  end
  # qを100000倍して整数にする。ファイル名に利用
  #q_int = Int(round(q * 100000))
  # qをメイン部分と桁に分離(有効数字4桁）
  qval, qexpo = format_q(q)
  # ほどよいacceptanceがわかっていたらNtauファイルから読み込み、
  # さもなければ、初期のNtauは1に指定。acceptanceを見ながら調整する
  if !isdir("Ntau")
      mkdir("Ntau")
  end
  Ntaufile = "Ntau/Ntau_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo)u00.txt"
  if isfile(Ntaufile) && filesize(Ntaufile) != 0
      open(Ntaufile, "r") do f
          Ntau = parse(Int, readline(f))
      end
      Trig = 0
  # Ntauファイルが存在しない場合、Ntauを調整するので初期設定
  else
      Ntau = 1
      Trig = 1
      acc = 0.0
  end
  # configがない時は熱化処理（Ntauも再設定）
  filebody = "Conf/config_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo)u00"
  filename = "$(filebody).txt"
  if !isfile(filename) 
    # 最初にNtauを調整
    acc = 0.0
    Ntau = 1
    Trig = 1
    # 小さいstep_sizeで初期配位を作成
    HMC(0,Nc,gamma,q,1000,0.0001,Ntau,filebody)
    HMC(1,Nc,gamma,q,1000,0.001,Ntau,filebody)
    HMC(1,Nc,gamma,q,1000,0.01,Ntau,filebody)
    while acc < MinAcc && Ntau <= MaxNtau
      niter0 = 1000
      ## path表示を利用(uを指定しないとpath表示を利用する仕様)
      acc, Uconf = HMC(1,Nc,gamma,q,niter0,step_size,Ntau,filebody)
      if acc < 0.10
        Ntau += 3
        println("try again: Ntau:",Ntau," q:",@sprintf("%.4f",q), " acc was ", @sprintf("%4f",acc))
      elseif acc < MinAcc 
        Ntau += 1
        println("try again: Ntau:",Ntau," q:",@sprintf("%.4f",q), " acc was ", @sprintf("%4f",acc))
      end
    end
    if Ntau > MaxNtau
      open("Ex_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo).txt", "w") do f
        # Write any necessary content to the file here
      end
      return 
    else 
      # 熱化
      niter0 = 40000
      #acc, phases, A = HMC(1,Nc,gamma,q,niter0,step_size,Ntau,filebody)
      HMC(1,Nc,gamma,q,niter0,step_size,Ntau,filebody)
    end
  end
  if Ntau > MaxNtau
    open("Ex_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo).txt", "w") do f
      ## Write any necessary content to the file here
    end
    return
  else
    # 本番前のNtau調整
    if Trig == 1 
      acc = 0.0
      #while acc < MinAcc
      #  niter0 = 1000
      #  acc, Uconf = HMC(1,Nc,gamma,q,niter0,step_size,Ntau,filebody)
      #  if acc < MinAcc
      #    Ntau += 1
      #  end
      #end
      open(Ntaufile, "w") do f
        write(f, string(Ntau))
      end
    end
    # 本番
    acc, Uconf = HMC(1,Nc,gamma,q,niter+10000,step_size,Ntau,filebody)

    Uconf = Uconf[1001:end]
    #phases = phases[1001:end]
    phases = [[] for _ in 1:RANK] # Uの固有値の偏角
    Wilson = [] # Wilson loop
    for U in Uconf
      tmp = []
      for a in 1:RANK
        eigs = eigvals(U[a])
        push!(tmp, sum(eigs)/Nc)
        push!(phases[a], map(angle,eigs))
      end
      push!(Wilson,tmp)
    end

    #########################
    # phaseのヒストグラムを保存
    plt=histogram(vcat(phases[1]...),xlim=(-3.1416,3.1416),bins=100,label="Cycle 1",alpha=0.5,normalize=true)
    for a in 2:RANK
         histogram!(vcat(phases[a]...),xlim=(-3.1416,3.1416),bins=100,label="Cycle $(a)",alpha=0.5,normalize=true)
    end
    xlabel!("\$\\theta\$")
    title!("$(NAME), \$\\gamma = $(gamma), q=$(q)\$, \$u=0.0\$")
    figname = "phases_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo)u00"
    backup_file(figname,"png")
    if !isdir("Phases")
        mkdir("Phases")
    end
    savefig(plt, "Phases/$(figname).png")


    ####################################
    # phaseのヒストグラムデータを書き出す
    hist = []
    # ビンのエッジを定義
    bin_edges = range(-3.1416, stop=3.1416, length=101) # length=101 to create 100 bins
    for a in 1:RANK
        phase = reduce(vcat,phases[a])
        h = fit(Histogram, phase, bin_edges)
        push!(hist, h.weights ./ length(phase) )
    end
    # データフレームの作成
    colnames = ["Cycle " * string(a) for a in 1:RANK]
    df = DataFrame(Dict(colnames[a] => hist[a] for a in 1:RANK))
    # CSVファイルへの書き出し
    CSV.write("Phases/$(figname).csv", df)

    ####################################
    # 物理量
    # actionを計算
    A = Float64[]
    for U in Uconf
      push!(A, action(U, Nc, q))
    end
    ## averages
    aveS = mean(A)
    ## moment
    devS = A .- aveS
    devS2 = devS.^2 
    devS3 = devS.^3

    ## 測定量の書き出し
    if !isdir("Obs")
        mkdir("Obs")
    end
    # energy
    filename = "Obs/energy_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo)u00_S$(SSint)Ntau$(Ntau)"
    backup_file(filename,"txt")
    write_realvalues("$(filename).txt", A ./ Nc^2)
    # 比熱
    filename = "Obs/specificheat_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo)u00_S$(SSint)Ntau$(Ntau)"
    backup_file(filename,"txt")
    write_realvalues("$(filename).txt", devS2 .* (gamma^2 / Nc^2) )
    # 比熱の温度微分
    filename = "Obs/dC_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo)u00_S$(SSint)Ntau$(Ntau)"
    backup_file(filename,"txt")
    write_realvalues("$(filename).txt", devS3 .* (gamma^3 / Nc^2 ) )
    # Wilson loop
    filename = "Obs/Wilson_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo)u00_S$(SSint)Ntau$(Ntau)"
    backup_file(filename,"txt")
    open("$(filename).txt", "w") do f
      for W in Wilson
        components = [ "$(real(W[a])), $(imag(W[a]))" for a in 1:RANK ] # list comprehension to create a string for each complex number
        line = join(components, ", ")  # join all components with a comma
        write(f, line, "\n") # write the line to the file
        #for a in 1:RANK
          #write(f, real(W[a]), ", ", imag(W[a]), ", ")
        #end
        #write(f,"\n")
      end
    end

    #########################
    # dCのhistoryを書き出し
    X = [1:length(devS3)]
    Ph = devS3 .* (gamma^3 / Nc^2 )
    plt2 = scatter(X,Ph)
    title!("$(NAME), dC, \$\\gamma = $(gamma), q=$(q)\$, \$u=0.0\$")
    figname = "dC_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo)u00"
    backup_file(figname,"png")
    if !isdir("Obs")
        mkdir("Obs")
    end
    savefig(plt2, "Obs/$(figname).png")


    #########################
    # configurationを書き出す(JLDを利用。Juliaのみ対応なので注意を)
    if !isdir("Uconfig")
      mkdir("Uconfig")
    end
    config_file = "Uconfig/Uconf_N$(Nc)g$(gamma_int)q$(qval)E$(qexpo)u00"
    backup_file(config_file,"jld")
    jldopen("$(config_file).jld", "w"; compress = true) do f 
      f["Uconf"] = Uconf
      #f["large_array"] = zeros(10000)
      #save("$(config_file).jld", "Uconf", Uconf)
    end
  end
end

#########################################
## Qの範囲を設定
## 相転移点±0.3 : 0.05刻み
## それ以外：0.2刻み
## qc = [qc1, qc2, ...]
## dqc = [dqc1, dqc2, ...]
function addS(S, ini, fin, diff)
  append!(S, collect(ini:diff:fin))
end


#########################################
## 通常のQを設定
function Q_regular(qc,qcd)
  # 近似的な相転移点の位置
  sc = sort(sof.(qc))
  scd = sort(sof.(qcd))
  # データを取る q の値を設定
  Smin = scd[1]-1.0
  Smax = sc[end]+1.0
  S = Float64[]
  addS(S,Smin,scd[1]-0.3,0.2)
  addS(S,scd[1]-0.3,scd[1],0.05)
  for i in 1:length(scd)-1
    if scd[i] + 0.3 < scd[i+1] - 0.3 
      S = addS(S,scd[i],scd[i+1],0.05)
    else
      S = addS(S,scd[i],scd[i]+0.3,0.05)
      S = addS(S,scd[i]+0.3,scd[i+1]-0.3,0.2)
      S = addS(S,scd[i+1]-0.3,scd[i+1],0.05)
    end
  end
  addS(S,scd[end],scd[end]+0.3,0.05)
  addS(S,scd[end]+0.3, sc[1]-0.3,0.2)
  addS(S,sc[1]-0.3,sc[1],0.05)
  for i in 1:length(sc)-1
    if sc[i] + 0.3 < sc[i+1] - 0.3 
      S = addS(S,sc[i]+0.05,sc[i+1],0.05)
    else
      S = addS(S,sc[i]+0.05,sc[i]+0.3,0.05)
      S = addS(S,sc[i]+0.3,sc[i+1]-0.3,0.2)
      S = addS(S,sc[i+1]-0.3,sc[i+1],0.05)
    end
  end
  addS(S,scd[end],scd[end]+0.3,0.05)
  addS(S,scd[end]+0.3,Smax,0.2)
  S = unique(S)
  Q = qof.(S)
  Q = filter(x -> !(1/OMEGA-0.01 < x < 1.0+0.01), Q)
  return Q
end

#########################################
## 相転移周辺を細かく取るQを設定
function Q_fine(qc,dqc)
  sc = sort(sof.(qc))
  scd = sort(sof.(dqc))
  
  S = Float64[]
  #for s in sc
  for i in eachindex(sc)
    if i == 1 || ( i > 1 && sc[i] > sc[i-1]+0.05)
      s = sc[i]
      addS(S, s+0.01,s+0.041,0.01)
      addS(S, s+0.06,s+0.091,0.01)
      addS(S, s-0.09,s-0.060,0.01)
      addS(S, s-0.04,s-0.010,0.01)
    end
  end
  ##
  #for s in scd 
  for i in eachindex(scd)
    if i == 1 || ( i > 1 && scd[i] > scd[i-1]+0.05)
      s = scd[i]
      addS(S, s+0.01,s+0.041,0.01)
      addS(S, s+0.06,s+0.091,0.01)
      addS(S, s-0.09,s-0.059,0.01)
      addS(S, s-0.04,s-0.009,0.01)
    end
  end
  ##
  S = unique(S)
  Q = qof.(S)
  return Q
end

# realの観測量を書き出す
function write_realvalues(filename, RV)
  open(filename, "w") do f
      for val in RV
          write(f, join(real(val), ", "))
          write(f,"\n")
      end
  end
end

##############################################################
function chunk(arr, n)
  k, r = divrem(length(arr), n)
  [arr[(i-1)*n+1:i*n] for i=1:k]..., arr[k*n+1:end]
end

end





        