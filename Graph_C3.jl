################################
## Cycle Graph n=3 ; C3
module Graph
# 外から見える変数
export NAME, OMEGA, minl, NV, NE, RANK, e_free, S, T, poly
# 外から見える関数
export Vmat, sof, qof

const NAME = "\$C_3\$"
const OMEGA = 1
const minl = 3
const n = minl
const NV = minl
const NE = minl
const RANK = 1 # gauge fixingで残るUの数。
const poly = [2]
const e_free = [1] # gauge fixingで残すedge

# cycle graphのsourceとtarget
S = []
T = []
for e=1:n
# Source of edges
    push!(S,e)
# Target of edges
    if e == n
        push!(T,1)
    else
        push!(T,e+1)
    end
end

# (q,u)からsへの変換
sof(q,u=0.0) = -log(q)
# (s,u)からqへの変換
qof(s,u=0.0) = exp(-s)

############################################################################
# Path Ihara Zeta の matrix V
function Vmat(q)
    V = zeros(ComplexF64,2,2)
    ###
    V[1,1] = q^n
    ###
    V[2,2] = q^n
    return V
end
end


