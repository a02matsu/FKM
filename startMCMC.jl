using Distributed
NPROCS = 18 
addprocs(NPROCS)
nworkers()
@everywhere include("Simulation.jl")
@everywhere using .Simulation
# それぞれのプロセスで乱数を初期化
rngs = @sync @distributed (vcat) for id in workers()
    init_rng(id)
end

for NN in [4,8]
  @everywhere begin
    global Nc = $NN 
  end
@everywhere begin
  # パラメータ設定とデータを取るQの設定
  # qc, qcd : 大まかな相転移点の位置を表す配列
  gamma = 128.0
  #Nc = 16
  qc = [(1.0/(2*gamma-1))^(1/3)]
  sc = sof.(qc)
  scd = 1.0 .- sc
  qcd = qof.(scd)
  # まずは大雑把に
  Q = Q_regular(qc,qcd)
  # 相転移点まわりを取る
  #qc = [0.3151]
  #51d = [2.074]
  #51= Q_fine(qc,qcd)
end
## QをNPROCS個ずつのブロックに分割
blocks = chunk(Q, NPROCS) 
println("$(length(Q)), $Q")
## Simulation Start
for block in blocks 
  @distributed for q in block
    start_simulation(Nc,gamma,q)
  end
end
end #NN