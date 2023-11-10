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
phys = "energy"

plt = plot(xlabelfontsize = 14, ylabelfontsize = 14, size=(1000,700), margin=20px)

#for gamma_int in [800,1600,12800,102400,409600,1638400]
for gamma_int in [102400,409600,1638400]

# CSVファイルの読み込み
file = "./FKM_$(NAME)/$(phys)_N$(Nc)g$(gamma_int)u00.csv"
df = CSV.read(file, DataFrame)

# プロットを作成
# 色のインデックスを追跡
cindex = 1
# 横軸、縦軸、エラーバーのデータを計算
#q,mean_val,stderr_val
x = sof.(df.q,0.0) # 横軸は log(q)
y = df.mean_val
err = df.stderr_val
# プロットを追加
g = Int(gamma_int/100)
scatter!(plt, x, y, yerror=err, label="\$\\gamma\$ = $(g)", markershape=:circle, markersize=3)

end

# 臨界帯
plot!(x -> 0.0, 0.0, 1.0, fillrange=x -> 10.0, linealpha=0,
  fillcolor=RGB(0.8, 0.8, 1), alpha=0.5, label="critical strip\n(unstable region)")

# 軸ラベルを設定
xlims!(3.0, 6.5)
ylims!(-15, 1)
xlabel!("\$s\$")
ylabel!("$(phys)")
plot!(plt, legend=:bottomright)

# 理論値
zeta(q) = 1/((1 - q)*(1 - 2q)*(1 - q^2)^2*(1 + q + 2q^2)^3)

g = 1024
E0(s) = -g^2*log(zeta(2^(-2*s)))
E1(s) = -g*log(zeta(2^(-s))) + 3/2
E3(s) = (E0(s)+E1(s))/2
plot!(3.0:0.01:6.5, E0, label="E0(g=1024)")
plot!(3.0:0.01:6.5, E1, label="E1(g=1024)")
plot!(3.0:0.01:6.5, E3, label="(E0+E1)/2 (g=1024)")

g = 4096
E0(s) = -g^2*log(zeta(2^(-2*s)))
E1(s) = -g*log(zeta(2^(-s))) + 3/2
E3(s) = (E0(s)+E1(s))/2
plot!(3.0:0.01:6.5, E0, label="E0(g=4096)")
plot!(3.0:0.01:6.5, E1, label="E1(g=4096)")
plot!(3.0:0.01:6.5, E3, label="(E0+E1)/2 (g=1024)")

g = 16384
E0(s) = -g^2*log(zeta(2^(-2*s)))
E1(s) = -g*log(zeta(2^(-s))) + 3/2
E3(s) = (E0(s)+E1(s))/2
plot!(3.0:0.01:6.5, E0, label="E0(g=16384)")
plot!(3.0:0.01:6.5, E1, label="E1(g=16384)")
plot!(3.0:0.01:6.5, E3, label="(E0+E1)/2 (g=1024)")

# プロットを表示
display(plt)
# プロットをファイルに保存
savefig(plt, "$(NAME)_$(phys).pdf")
