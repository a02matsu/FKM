using Distributed
NPROCS = 18 
addprocs(NPROCS)
nworkers()
@everywhere include("Simulation.jl")
@everywhere begin
  using .Simulation
  # パラメータ設定とデータを取るQの設定
  # qc, qcd : 大まかな相転移点の位置を表す配列
  gamma = 4.0
  Nc = 8
  qc = [(1.0/(2*gamma-1))^(1/3)]
  sc = sof.(qc)
  scd = 1.0 .- sc
  qcd = qof.(scd)
  # まずは大雑把に
  Q = Q_regular(qc,qcd)
  # 相転移点まわりを取る
  # qc = []
  # qcd = []
  #Q = Q_fine(qc,qcd)
end
# それぞれのプロセスで乱数を初期化
rngs = @sync @distributed (vcat) for id in workers()
    init_rng(id)
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