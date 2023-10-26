###################################
# 物理量のプロット
# 同じパラメータのファイルを集めて、統計処理をしてcsvに描き込む
function calc_obs(phys,Nc, gamma_int)
    pattern = "Obs/$(phys)_N$(Nc)g$(gamma_int)q*u00*.txt"
    file_paths = glob(pattern)
    
    results = Dict{Float64, Vector{Float64}}()
    
    for file_path in file_paths
        q_int = parse(Int, match(r"q(\d+)", file_path).captures[1])
        data = vec(readdlm(file_path, Float64))
        
        if haskey(results, q_int/10000)
            append!(results[q_int/10000], data)
        else
            results[q_int/10000] = data
        end
    end
    
    summary = DataFrame(q=Float64[], mean_val=Float64[], stderr_val=Float64[])
    
    for (q, values) in results
        mean_val = mean(values)
        stderr_val = std(values) / sqrt(length(values))
        
        push!(summary, (q, mean_val, stderr_val))
    end
    
    sort!(summary, :q)
    
    CSV.write("$(phys)_N$(Nc)g$(gamma_int)u00.csv", summary)
end
# 物理量の統計データを作成
for phys in ["energy", "specificheat", "dC"]
    file_path="$(phys)_N$(Nc)g$(gamma_int)u00.csv"
    fig="$(phys)_N$(Nc)g$(gamma_int)u00.png"
    calc_obs(phys,Nc, gamma_int)
end

