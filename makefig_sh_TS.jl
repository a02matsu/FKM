#
# 論文用の図を作成
#

using CSV
using DataFrames
using Plots
using Plots.PlotMeasures
using Glob
include("Graph.jl")
using .Graph

# カラーリストを定義
ColorList = [:blue, :red, :green, :cyan, :purple, :yellow]

Nc = 16
phys = "specificheat"
# 理論値(TS)
zeta(q) = -1 / ((q+1) * (2*q^4+3*q^3+3*q^2+2*q+1) * (2*q^5-q^4+2*q^3-q^2+q-1) * (q-1)^2)

plt = plot(xlabelfontsize=14, ylabelfontsize=14, size=(1000, 700), margin=20px)

#for gamma_int in [400, 800, 1600, 12800, 102400, 409600, 1638400, 13107200]
for gamma_int in [102400]

  # CSVファイルの読み込み
  file = "./FKM_$(NAME)/$(phys)_N$(Nc)g$(gamma_int)u00.csv"
  df = CSV.read(file, DataFrame)

  # プロットを作成
  # 色のインデックスを追跡
  cindex = 1
  # 横軸、縦軸、エラーバーのデータを計算
  #q,mean_val,stderr_val
  x = sof.(df.q, 0.0) # 横軸は log(q)
  y = df.mean_val
  err = df.stderr_val
  # プロットを追加
  g = Int(gamma_int / 100)
  scatter!(plt, x, y, yerror=err, label="HMC", markershape=:circle, markersize=3)

  omega = 1.42405
  q(s) = omega^(-s)
  C0(s) = g^2 * log(zeta(q(s)^2))
  C1(s) = g^2 * 2 * q(s)^8 + 1/2
  C2(s) = 1.0
  plot!(1.0:0.05:5.39216, C2, label="analytic", linecolor=:green, linewidth=3.0)
  plot!(5.39216:0.05:7.18955, C1, label=nothing, linecolor=:green, linewidth=3.0)
  plot!(7.18955:0.05:10.0, C0, label=nothing, linecolor=:green, linewidth=3.0)

  qd(s) = omega^(-(1-s))
  C0d(s) = g^2 * log(zeta(qd(s)^2))
  C1d(s) = g^2 * 2 * qd(s)^8 + 1/2
  C2d(s) = 1.0
  plot!(-4.39216:0.05:0, C2d, label=nothing, linecolor=:green, linewidth=3.0)
  plot!(-6.18955:0.01:-4.39216, C1d, label=nothing, linecolor=:green, linewidth=3.0)
  plot!(-9.18955:0.05:-6.18955, C0d, label=nothing, linecolor=:green, linewidth=3.0)
end

# 臨界帯
plot!(x -> 0.0, 0.0, 1.0, fillrange=x -> 10.0, linealpha=0,
  fillcolor=RGB(0.8, 0.8, 1), alpha=0.5, label="critical strip\n(unstable region)")

# 軸ラベルを設定
xlims!(-9.0, 10.0)
ylims!(0, 1.5 )
xlabel!("\$s\$")
ylabel!("$(phys)")
plot!(plt, legend=:topright)

# プロットを表示
display(plt)
# プロットをファイルに保存
savefig(plt, "$(NAME)_$(phys).pdf")
