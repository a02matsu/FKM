using Glob
using Plots
using DelimitedFiles

files = glob("dC_*.txt")

for file in files 
  Nobj = match(r"N(\d+)", file)
  gobj = match(r"g(\d+)", file)
  qobj = match(r"q(\d+)", file)
  Nc = parse(Int, Nobj.captures[1])
  gamma = parse(Int, gobj.captures[1])
  q = parse(Int, qobj.captures[1])
  figfile = "dC_N$(Nc)g$(gamma)q$(q).png"

  data = readdlm(file, Float64)
  plt = scatter(1:length(data), data)
  savefig(plt, figfile)
end


