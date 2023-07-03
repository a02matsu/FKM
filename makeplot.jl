include("Measurement.jl")
using .Measurement



#################################################
# 物理量の平均値をまとめたCSVファイルを作る
# fundamental cycleのWilson loopの平均値をまとめたCSVファイルを作る
for phys in ["energy", "specificheat", "dC"]
  for Nc in [4,8]
    for gamma in [4,8,16,128,1024]
        calc_obs(phys,Nc, gamma)
        calc_Wilson(Nc, gamma)
    end
  end
end

#################################################
# gammaを固定したプロットを作成
#for G in [4,8,16,128,1024]
for G in [1024,16384]
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
for g in [4,16,128,1024]
  plot_gamma_small_q("specificheat",g)
end

############################################################
# large q 領域だけのプロット
for g in [1024, 16384]
  plot_gamma_large_q("specificheat",g)
  plot_gamma_large_q("dC",g)
end
