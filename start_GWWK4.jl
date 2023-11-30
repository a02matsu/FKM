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
  # パラメータ設定とデータを取るaの値を設定
  Nc = 16
  # まずは大雑把に
  A = range(0.35,0.45,18)
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