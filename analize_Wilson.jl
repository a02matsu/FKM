include("Measurement.jl")
using .Measurement

##########################
## 調べたいサイクル
cycles = [[1,2,3],[1],[3],[4],[2],[2,8],[1,2],[2,3],[1,2,3,1,2,3],[1,3],[1,4],[3,4],[1,4,3],[1,2,3,1,2,3,1,2,3],[1,2,3,1,2,3,1,2,3,1,2,3],[1,1],[1,1,1],[2,2],[2,2,2],[3,3],[3,3,3],[4,4],[4,4,4]]
Cnames = []
for i in eachindex(cycles)
  push!(Cnames,"C$(join(cycles[i],""))")
end

for gamma in [1024,16384]
  PhaseHist(16,gamma,cycles)
end

Qc1,err1,Qc2,err2,Q = findPT2(16,1024)
sc1=sof.(Qc1)
sc2=sof.(Qc2)
err_sc1 = abs.(sof.(Qc1 .+ err1) - sof.(Qc1))
err_sc2 = abs.(sof.(Qc2 .+ err2) - sof.(Qc2))
sc1 .+ sc2

plt = plot()
scatter!(plt, Cnames, sc1, yerror=err_sc1, label="" )
scatter!(plt, Cnames, sc2, yerror=err_sc2, label="" )
ylabel!("\$ s \$")
savefig(plt,"C_vs_sc.png")


function plot_critical(Nc,gamma,cycle, cycles, Qc1, Qc2,Q)
  ind = findfirst(x->x==cycle, cycles)
  target1 = Qc1[ind]
  target2 = Qc2[ind]
  _, ind1 = findmin(abs.(Q .- target1)) 
  _, ind2 = findmin(abs.(Q .- target2)) 
  q1 = Q[ind1]
  q2 = Q[ind2]
  exp1 = floor(Int, log10(abs(q1))) #-leading_zeros(q)
  val1 = Int(round(q1*10^(-exp1+3)))
  q1txt = "$(val1)E$(exp1)"
  exp2 = floor(Int, log10(abs(q2))) #-leading_zeros(q)
  val2 = Int(round(q2*10^(-exp2+3)))
  q2txt = "$(val2)E$(exp2)"

  x = range(-3.1416, stop=3.1416, length=100)
  Cname = "cycle$(join(cycle,"_"))"
  plt = plot()

  ###
  file = "Phases/Cphases_N$(Nc)g$(gamma)q$(q1txt)u00.csv"
  df = CSV.read(file, DataFrame)
  y = df[:, Cname]
  bar!(plt, x, y, label="q=$q1txt",alpha=0.5)
  ###
  file = "Phases/Cphases_N$(Nc)g$(gamma)q$(q2txt)u00.csv"
  df = CSV.read(file, DataFrame)
  y = df[:, Cname]
  bar!(plt, x, y, label="q=$q2txt",alpha=0.5)
  ylabel!("s")
  title!("\$ U_{$(join(cycle,""))} \$")
end

plot_critical(16,1024,[2,8],cycles,Qc1,Qc2,Q)




