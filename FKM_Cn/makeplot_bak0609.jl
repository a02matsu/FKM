include("Measurement.jl")
using .Measurement

#################################################
# 物理量の平均値をまとめたCSVファイルを作る
for phys in ["energy", "specificheat", "dC"]
  for Nc in [4,8,16]
    for gamma_int in [400,800,1600,12800,102400]
        file_path="$(phys)_N$(Nc)g$(gamma_int)u00.csv"
        fig="$(phys)_N$(Nc)g$(gamma_int)u00.png"
        calc_obs(phys,Nc, gamma_int)
    end
  end
end


#################################################
# gammaを固定したプロットを作成
for G in [400,800,1600,12800,102400]
  plot_gamma("energy", G)
  plot_gamma("specificheat", G)
  plot_gamma("dC", G)
end

#################################################
# Nを固定したプロットを作成
for N in [16]
  plot_N("specificheat", N)
  plot_N("dC", N)
end


#################################################
## 相転移の位置を解析
analysePT(16, 10^(-4), 0.025)


############################################################
# small q 領域だけのプロットを作成してCnと比較
for g in [400,1600,12800,102400]
  plot_gamma_small_q("specificheat",g)
end

