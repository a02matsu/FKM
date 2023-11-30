using Distributed
NPROCS = 12  
addprocs(NPROCS)
nworkers()
@everywhere include("Simulation.jl")
@everywhere using .Simulation
# それぞれのプロセスで乱数を初期化
rngs = @sync @distributed (vcat) for id in workers()
    init_rng(id)
end

for NN in [16]
  @everywhere begin
    global Nc = $NN 
  end
@everywhere begin
  # aの範囲を設定
  A = range(0.1,0.4,36)
end
## QをNPROCS個ずつのブロックに分割
blocks = chunk(A, NPROCS) 
println("$(length(A)), $A")
## Simulation Start
for block in blocks 
  @distributed for aa in block
    start_simulation_GWWK4(Nc,aa)
  end
end
end #NN