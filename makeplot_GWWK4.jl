using Printf
using Glob
using DataFrames
using CSV
using DelimitedFiles
using StatsBase
using Plots

function calc_obs(phys,Nc::Int)
  pattern = "Obs/$(phys)_GWWK4_N$(Nc)a*.txt"
  file_paths = glob(pattern)
    
  results = Dict{String, Vector{Float64}}()
    
  for file_path in file_paths
    atxt = match(r"a([0-9].+)_", file_path).captures[1]
    data = vec(readdlm(file_path, Float64))

    if haskey(results, atxt)
      append!(results[atxt], data)
    else
      results[atxt] = data
    end
  end
    
  summary = DataFrame(a=String[], mean_val=Float64[], stderr_val=Float64[])
  
  for (atxt, values) in results
      mean_val = mean(values)
      stderr_val = std(values) / sqrt(length(values))
      
      push!(summary, (atxt, mean_val, stderr_val))
  end
    
  sort!(summary, :a)
    
  CSV.write("Obs/$(phys)_GWWK4_N$(Nc).csv", summary)
end

for phys in ["energy", "specificheat", "dC"]
  for Nc in [8]#,8,16]
    calc_obs(phys,Nc)
  end
end

file = "Obs/specificheat_GWWK4_N8.csv"
df = DataFrame(CSV.File(file))
x = df.a
y = df.mean_val
err = df.stderr_val
plt = plot()
xlabel!(plt,"\$x\$")
scatter!(plt, x, y, yerror=err,label="C")