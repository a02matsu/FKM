## Ncとgammaを与えると、
# 1) 理論的な相転移点の情報に基づいて大雑把なsimulationを実行し、
# 2) その結果から相転移点を割り出し
# 3) その相転移点周辺のシミュレーションを細かく行う
using Distributed
include("Simulation.jl")
include("Measurement.jl")
using .Simulation
using .Measurement: findPT

# プロセス数
NPROCS = 18 
# シミュレーションのパラメータ
#Ns = [4,8,16]
#Gammas = [8.0,16.0,128.0,1024.0]
#Gammas = [16384.0]
Gammas = [Float64(2^30)]
Ns = [8,16]
#gamma = 2^30
#sof((2*gamma)^(-1/4))
#println(length(Q_regular([(2*gamma)^(-1/3),(2*gamma)^(-1/4)],[(2*gamma/8)^(1/3),(2*gamma/16)^(1/4)])))

function parallel_simulation(Nc::Int, gamma::Float64, Q::Vector{Float64}, niter=200000)
  # QをNPROCS個ずつのブロックに分割
  blocks = chunk(Q, nprocs())
  futures = [] # added
  ## Simulation Start
  for block in blocks 
    #@sync @distributed for q in block
    f = @distributed for q in block
      start_simulation(Nc, gamma, q, niter)
    end
    push!(futures, f)
  end
  return futures
end

for g in Gammas
  for N in Ns
    #################################################
    ## プロセスを追加してシミュレーションを準備
    addprocs(NPROCS)
    for pid in workers()
      @spawnat pid include("Simulation.jl")
    end
    @everywhere begin
      using .Simulation
      global Nc = $N
      global gamma = $g
    end
    #################################################
    # CHANGE HERE 
    qc = [(2*gamma)^(-1/3),(2*gamma)^(-1/4)]
    qcd = [(2*gamma/8)^(1/3),(2*gamma/16)^(1/4)]
    # まずは大雑把に
    Q = Q_regular(qc,qcd)
    futures = parallel_simulation(Nc,gamma,Q)
    for f in futures
      wait(f)
    end
    # メモリ管理のためにプロセスを消去
    rmprocs(workers())

    # Nc=16のときは細かく取る
    if Nc == 16 
      # 相転移点をみつける
      Qc1, Qc2 = findPT(Nc, gamma)
      qc = []
      qcd = []
      for a in 1:RANK
        push!(qc, (Qc1[a][1] .+ Qc2[a][1]) ./ 2 )
        push!(qcd, (Qc1[a][2] .+ Qc2[a][2]) ./ 2 )
      end
      println(sof.(qc))
      println(sof.(qcd))
      Q = Q_fine(qc,qcd)
      #################################################
      ## プロセスを追加してシミュレーションを準備
      addprocs(NPROCS)
      for pid in workers()
        @spawnat pid include("Simulation.jl")
      end
      @everywhere begin
        using .Simulation
        global Nc = $N
        global gamma = $g
      end
      #################################################
      futures = parallel_simulation(Nc,gamma,Q)
      for f in futures
        wait(f)
      end
    end
    # メモリ管理のためにプロセスを消去
    rmprocs(workers())
  end
end