## change log
# ver.1 （2023/05/30）
# dCの出力にγ^3をかけることにした

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
end
# それぞれのプロセスで乱数を初期化
rngs = @sync @distributed (vcat) for _ in workers()
    init_rng()
end
# 各種パラメータ
@everywhere begin
    ## Control Simulation
    GRAPH = "\$C_3\$"
    MaxNtau = 30 # limit of Ntau
    MinAcc = 0.75 # minimal acceptance ratio
    ## theory parameters
    Nc = 8
    u = 0e-1
    gamma = 1024.0
    # gammaを100倍して整数にし、ファイル名に利用
    gamma_int = Int(round(gamma*100))
    Nf = gamma * Nc 
    # シミュレーションのパラメータ
    niter = 200000
    step_size = 0.10
    SSint = Int(round(step_size*100)) 
end
# (q,u)からsへの変換
sof(q) = -log(q)
# (s,u)からqへの変換
qof(s) = exp(-s)
# 近似的な相転移点の位置
qc1 = (1.0/(2*gamma-1))^(1/3)
sc1 = sof(qc1)
scd1 = -sc1
qof(scd1)
# 安定領域の境界
qb1 = 0.0
qb2 = 0.0
# データを取る q の値を設定
# q = [ (1-u)(OMEGA-u) ]^(-s)
# として、sが等間隔になるようにデータを収集する
#sof(qb1-0.01,0.0)
#sof(qb2+0.01,0.0)
Smin = -3.0
Smax = 3.0
S1 = collect(Smin:0.2:scd1-0.3)
S2 = collect(scd1-0.3:0.05:scd1+0.3)
S3 = collect(scd1+0.3:0.2:sc1-0.3)
S4 = collect(sc1-0.3:0.05:sc1+0.3)
S5 = collect(sc1+0.3:0.2:Smax)
S = unique(vcat(S1,S2,S3,S4,S5))
Q = qof.(S)
Q = filter(x -> !(qb1-0.01 < x < qb2+0.01), Q)


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
        while acc < MinAcc && Ntau <= MaxNtau
            niter0 = 1000
            ## path表示を利用(uを指定しないとpath表示を利用する仕様)
            acc, phases, A = HMC(0,Nc,gamma,q,niter0,step_size,Ntau,filebody)
            if acc < 0.10
                Ntau += 3
                println("try again: Ntau:",Ntau," q:",@sprintf("%.4f",q))
            elseif acc < MinAcc 
                Ntau += 1
                println("try again: Ntau:",Ntau," q:",@sprintf("%.4f",q))
            end
        end
        if Ntau <= MaxNtau
          # 熱化
          niter0 = 10000
          acc, phases, A = HMC(1,Nc,gamma,q,niter0,step_size,Ntau,filebody)
        end
    end
    if Ntau > MaxNtau
        open("Ex_N$(Nc)g$(gamma_int)q$(q_int).txt", "w") do f
            # Write any necessary content to the file here
        end
    else
        # 本番前のNtau調整
        if Trig == 1 
            acc = 0.0
            while acc < MinAcc
                niter0 = 1000
                acc, phases, A = HMC(1,Nc,gamma,q,niter0,step_size,Ntau,filebody)
                if acc < MinAcc
                    Ntau += 1
                end
            end
            open(Ntaufile, "w") do f
                write(f, string(Ntau))
            end
        end
        # 本番
        acc, phases, A = HMC(1,Nc,gamma,q,niter,step_size,Ntau,filebody)
        put!(ACC_NTAU, [q, Ntau, acc])

        #########################
        # phaseのヒストグラムを保存
        plt=histogram(vcat(phases[1]...),xlim=(-3.1416,3.1416),bins=100,label="Cycle 1",alpha=0.5,normalize=true)
        for a in 2:RANK
             histogram!(vcat(phases[a]...),xlim=(-3.1416,3.1416),bins=100,label="Cycle $(a)",alpha=0.5,normalize=true)
        end
        xlabel!("\$\\theta\$")
        title!("$(GRAPH), \$\\gamma = $(gamma), q=$(q)\$, \$u=0.0\$")
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
    end
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
###################################
# 物理量のプロット
# 同じパラメータのファイルを集めて、統計処理をしてcsvに描き込む
#function calc_obs(phys,Nc, gamma_int)
#    pattern = "Obs/$(phys)_N$(Nc)g$(gamma_int)q*u00*.txt"
#    file_paths = glob(pattern)
#    
#    results = Dict{Float64, Vector{Float64}}()
#    
#    for file_path in file_paths
#        q_int = parse(Int, match(r"q(\d+)", file_path).captures[1])
#        data = vec(readdlm(file_path, Float64))
#        
#        if haskey(results, q_int/10000)
#            append!(results[q_int/10000], data)
#        else
#            results[q_int/10000] = data
#        end
#    end
#    
#    summary = DataFrame(q=Float64[], mean_val=Float64[], stderr_val=Float64[])
#    
#    for (q, values) in results
#        mean_val = mean(values)
#        stderr_val = std(values) / sqrt(length(values))
#        
#        push!(summary, (q, mean_val, stderr_val))
#    end
#    
#    sort!(summary, :q)
#    
#    CSV.write("$(phys)_N$(Nc)g$(gamma_int)u00.csv", summary)
#end
## 物理量の統計データを作成
#for phys in ["energy", "specificheat", "dC"]
#    file_path="$(phys)_N$(Nc)g$(gamma_int)u00.csv"
#    fig="$(phys)_N$(Nc)g$(gamma_int)u00.png"
#    calc_obs(phys,Nc, gamma_int)
#end

