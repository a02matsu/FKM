## Ncとgammaを与えると、
# 1) 理論的な相転移点の情報に基づいて大雑把なsimulationを実行し、
# 2) その結果から相転移点を割り出し
# 3) その相転移点周辺のシミュレーションを細かく行う
using Distributed

Ns = [4,8,16]
Gammas = [1024.0]

NPROCS = 6 
addprocs(NPROCS)
nworkers()
@everywhere include("Simulation.jl")
@everywhere include("Measurement.jl")
@everywhere function parallel_simulation(Nc::Int, gamma::Float64, Q::Vector{Float64}, niter=200000)
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
@everywhere begin
  using .Simulation
  using .Measurement: findPT
end

for g in Gammas
  for N in Ns
    @everywhere begin
      global Nc = $N
      global gamma = $g
    end
    @everywhere begin
      #  using .Simulation
      #  using .Measurement: findPT
      # パラメータ設定とデータを取るQの設定
      # qc, qcd : 大まかな相転移点の位置を表す配列
      #gamma = 4.0
      #Nc = 4
      qc = [(1.0/(2*gamma-1))^(1/3)]
      sc = sof.(qc)
      scd = 1.0 .- sc
      qcd = qof.(scd)

      # まずは大雑把に
      Q = Q_regular(qc,qcd)
      futures = parallel_simulation(Nc,gamma,Q)
      for f in futures
        wait(f)
      end
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
  
        # 相転移点まわりを取る
        Q = Q_fine(qc,qcd)
        parallel_simulation(Nc,gamma,Q)
      end
    end
  end
end