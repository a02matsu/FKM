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

begin 
  Nc = 16
  gamma = 128.0
  #qcd = qof.(scd)
  ss(q) = -log(q)/log(OMEGA) 
  qq(s) = OMEGA^(-s)
  #g* q^3 = 1/2 
  q1(g) = (2*g)^(-1/3)
  q2(g) = (2*g)^(-1/4)
  s1(g) = -log(q1(g))/log(OMEGA)
  s2(g) = -log(q2(g))/log(OMEGA)

  r = 0.4
  c2=s2(gamma)
  c1=s1(gamma)
  S1 = collect(range(1, stop=c2-r-0.02, step=0.1))
  S2 = collect(range(c2-r, stop=c1+r, step=0.02))
  S3 = collect(range(c1+r+0.02, stop=c1*2, step=0.1))

  # dual side
  # Mapleによると、dual sideの相転移点は
  # q1(16) := 1.965334605
  # q2(16) := 1.654675036
  # s1(16) := -1.610197267 (-1.935053609)
  # s2(16) := -1.200159538 (-1.576736900)
  #
  # q1(128) := 3.979867356
  # q2(128) := 2.807995890
  # s1(128) := -3.88 (-3.780395501)
  # s2(128) := -2.88 (-2.988822776)
  #
  # q1(1024) := 7.989666218
  # q2(1024) := 4.744011982
  # s1(1024) := -5.87 (-5.604362168)
  # s2(1024) := -4.39 (-4.384551980)
  #
  # q1(16384) := 20.15460804
  # q2(16384) := 9.507126407
  # s1(16384) := -8.49 (-8.032719080)
  # s2(16384) := -6.37 (-6.242770510)
  #
  # <以前の実験(dual side)>
  #4,8,16,128,1024
  # dsc = : [-0.7008,-1.2663,-1.8818,-3.7780,-5.6242]
  # dsc = : [-0.6007,-1.0663,-1.5818,-2.9779,-4.3934]
  # (-0.8757 +/- 0.0341758) log(γ-1/2) + (0.4651 +/- 0.1406483)
  # (-0.6701 +/- 0.0341758) log(γ-1/2) + (0.2599 +/- 0.1406483)


  c1 = -5.87 
  c2 = -4.39 
  SS1 = collect(range(2.5*c1, stop=c1-r-0.02, step=0.1))
  SS2 = collect(range(c1-r, stop=c2+r, step=0.02))
  SS3 = collect(range(c2+r+0.02, stop=0, step=0.1))

  S = sort(vcat(S1,S2,S3,SS1,SS2,SS3))
  Q = qq.(S)
end
#@everywhere func(q) = start_simulation(Nc,gamma,q)
#pmap(func, Q)
pmap(q -> start_simulation(Nc,gamma,q), Q)
## QをNPROCS個ずつのブロックに分割
#blocks = chunk(Q, NPROCS) 
#println("$(length(Q)), $Q")
## Simulation Start
#for block in blocks 
#  @distributed for q in block
#    start_simulation(Nc,gamma,q)
#  end
#end
#end #NN
rmprocs(workers())
