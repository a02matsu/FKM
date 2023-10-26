include("Measurement.jl")
using .Measurement

plot_ST(13107200)
plot_completed_energy(13107200)

#################################################
# 物理量の平均値をまとめたCSVファイルを作る
for phys in ["energy", "specificheat", "dC"]
  for Nc in [8]
    for gamma_int in [1638400]
        calc_obs(phys,Nc, gamma_int)
    end
  end
end


#################################################
# gammaを固定したプロットを作成
for G in [1638400]
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
analysePT(16, 10^(-3), 10^(-4))

############################################################
# small q 領域だけのプロットを作成してCnと比較
for g in [400,1600,12800,102400]
  plot_gamma_small_q("specificheat",g)
end

############################################################
# large q 領域だけのプロット
for g in [102400, 1638400]
  plot_gamma_large_q("specificheat",g,-6.0,-3.0)
  plot_gamma_large_q("dC",g,-6.0,-3.0)
end
