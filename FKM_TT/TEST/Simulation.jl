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

using Distributed
# ローカルマシンのプロセスにワーカープロセスを追加する
NPROCS = 18 
addprocs(NPROCS)
# リモートマシンにワーカーを追加
#remote_username = "matsu"
#remote_hostname = "192.168.1.2"
#dirname="/Users/matsu/Dropbox/JuliaPrograms/FKM_TS"
#exename = "/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia"
#NPROCS = 6
#for _ in 1:NPROCS
#    addprocs(["$remote_username@$remote_hostname"], exename="$exename", dir="$dirname", tunnel=true)
#end

# nworkersを定義する
nworkers()

@everywhere include("FKM.jl")
@everywhere begin
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

    # 乱数を初期化する関数。毎回違うseedを使うためにtime()を使う。
    function init_rng()
        seed = round(Int, time()) % 6599 * myid()
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

    # 
end
# それぞれのプロセスで乱数を初期化
rngs = @sync @distributed (vcat) for _ in workers()
    init_rng()
end
# 各種パラメータ
@everywhere begin
    ## Control Simulation
    MaxNtau = 30 # limit of Ntau
    MinAcc = 0.75 # minimal acceptance ratio
    ## theory parameters
    Nc = 4  
    u = 0e-1
    gamma = 16384.0
    # gammaを100倍して整数にし、ファイル名に利用
    gamma_int = Int(round(gamma*100))
    Nf = gamma * Nc 
    # シミュレーションのパラメータ
    niter = 20000
    step_size = 0.10
    SSint = Int(round(step_size*100)) 
end

# 安定領域の境界
qb1 = 1.0/OMEGA
qb2 = 1.0
# 近似的な相転移点の位置
#qc1 = (1.0/(2*gamma-1))^(1/3)
#qc2 = 0.842*gamma^(0.319)
#qc1 =  0.26/log(gamma-1/2)
qc1 =  (0.50/(gamma-1/2))^(1/3)
#qc2 =  (0.50/(gamma-1/2))^(1/4)
qcd1 = 13.385
qcd2 = 16.000
sc1 = sof(qc1)
scd1 = sof(qcd1)
scd2 = sof(qcd2)
# データを取る q の値を設定
# q = [ (1-u)(OMEGA-u) ]^(-s)
# として、sが等間隔になるようにデータを収集する
#sof(qb1-0.01,0.0)
#sof(qb2+0.01,0.0)
Smin = -5.0
Smax = 6.0
S = Float64[]
append!(S, collect(Smin:0.2:scd1-0.3) )
S = unique(S)
Q = qof.(S)
Q = filter(x -> !(qb1-0.05 < x < qb2+0.05), Q)

#############################################################
#Q=vcat(Q1,Q2)
# acceptanceとNTAUを格納する共有配列を作成
ACC_NTAU = RemoteChannel(()->Channel{Vector{Float64}}(Inf))

Ind = collect(1:length(Q))
zipped_data = collect(zip(Q, Ind))

##############################################################
function chunk(arr, n)
    k, r = divrem(length(arr), n)
    [arr[(i-1)*n+1:i*n] for i=1:k]..., arr[k*n+1:end]
end
##############################################################
## zip(Q, Ind)をNPROCS個ずつのブロックに分割
blocks = chunk(zipped_data, NPROCS) 

#@distributed for (q,i) in zipped_data
for block in blocks 
  @distributed for (q,i) in block
    if !isdir("Conf")
        mkdir("Conf")
    end
    # qを1000倍して整数にする。ファイル名に利用
    q_int = Int(round(q * 10000))
    # ほどよいacceptanceがわかっていたらNtauファイルから読み込み、
    # さもなければ、初期のNtauは1に指定。acceptanceを見ながら調整する
    if !isdir("Ntau")
        mkdir("Ntau")
    end
    Ntaufile = "Ntau/Ntau_N$(Nc)g$(gamma_int)q$(q_int)u00.txt"
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
    filebody = "Conf/config_N$(Nc)g$(gamma_int)q$(q_int)u00"
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
          open("Ex_N$(Nc)g$(gamma_int)q$(q_int).txt", "w") do f
              # Write any necessary content to the file here
          end
          return 
        else 
        # 熱化(改めてcold start)
        niter0 = 40000
        #acc, phases, A = HMC(1,Nc,gamma,q,niter0,step_size,Ntau,filebody)
        HMC(0,Nc,gamma,q,niter0,step_size,Ntau,filebody)
      end
    end
    if Ntau > MaxNtau
        open("Ex_N$(Nc)g$(gamma_int)q$(q_int).txt", "w") do f
            ## Write any necessary content to the file here
        end
      return
    end
    # 本番前のNtau調整
    if Trig == 1 
        acc = 0.0
        while acc < MinAcc
            niter0 = 1000
            acc, Uconf = HMC(1,Nc,gamma,q,niter0,step_size,Ntau,filebody)
            if acc < MinAcc
                Ntau += 1
            end
        end
        open(Ntaufile, "w") do f
            write(f, string(Ntau))
        end
    end
    # 本番
    acc, Uconf = HMC(1,Nc,gamma,q,niter+10000,step_size,Ntau,filebody)
    put!(ACC_NTAU, [q, Ntau, acc])

    Uconf = Uconf[1001:end]
    #phases = phases[1001:end]
    phases = [[] for _ in 1:RANK]
    for U in Uconf
      for a in 1:RANK
        push!(phases[a], map(angle,eigvals(U[a])))
      end
    end

    #########################
    # phaseのヒストグラムを保存
    plt=histogram(vcat(phases[1]...),xlim=(-3.1416,3.1416),bins=100,label="Cycle 1",alpha=0.5,normalize=true)
    for a in 2:RANK
         histogram!(vcat(phases[a]...),xlim=(-3.1416,3.1416),bins=100,label="Cycle $(a)",alpha=0.5,normalize=true)
    end
    xlabel!("\$\\theta\$")
    title!("$(NAME), \$\\gamma = $(gamma), q=$(q)\$, \$u=0.0\$")
    figname = "phases_N$(Nc)g$(gamma_int)q$(q_int)u00"
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
    filename = "Obs/energy_N$(Nc)g$(gamma_int)q$(q_int)u00_S$(SSint)Ntau$(Ntau)"
    backup_file(filename,"txt")
    write_realvalues("$(filename).txt", A ./ Nc^2)
    # 比熱
    filename = "Obs/specificheat_N$(Nc)g$(gamma_int)q$(q_int)u00_S$(SSint)Ntau$(Ntau)"
    backup_file(filename,"txt")
    write_realvalues("$(filename).txt", devS2 .* (gamma^2 / Nc^2) )
    # 比熱の温度微分
    filename = "Obs/dC_N$(Nc)g$(gamma_int)q$(q_int)u00_S$(SSint)Ntau$(Ntau)"
    backup_file(filename,"txt")
    write_realvalues("$(filename).txt", devS3 .* (gamma^3 / Nc^2 ) )

    #########################
    # dCのhistoryを書き出し
    X = [1:length(devS3)]
    Ph = devS3 .* (gamma^3 / Nc^2 )
    plt2 = scatter(X,Ph,markersize=2)
    title!("$(NAME), dC, \$\\gamma = $(gamma), q=$(q)\$, \$u=0.0\$")
    figname = "dC_N$(Nc)g$(gamma_int)q$(q_int)u00"
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
    config_file = "Uconfig/Uconf_N$(Nc)g$(gamma_int)q$(q_int)u00"
    backup_file(config_file,"jld")
    save("$(config_file).jld", "Uconf", Uconf)

  end
end
###################################
# Acceptanceを表示
# acceptanceとNtauの値をqの値でsortする
#ACC_data = sort([take!(ACC_NTAU) for _ in 1:length(Q)], by = x -> x[1])
#using Printf
#for data in ACC_data
#    q, Ntau, acc = data
#    if acc < 0.75
#        printstyled(@sprintf("q=%.3f : acc=%.2f : Ntau=%d\n", q, acc, Ntau), color=:red)
#    else
#        printstyled(@sprintf("q=%.3f : acc=%.2f : Ntau=%d\n", q, acc, Ntau))
#    end
#    #println("q=$q : Ntau = $(Int(Ntau)), acc = $acc")
#end