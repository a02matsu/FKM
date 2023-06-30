################################
## Triangle-Prism (TP)
module Graph
# 外から見える変数
export NAME, OMEGA, minl, NV, NE, RANK, e_free, S, T, poly, NL
# 外から見える関数
export Vmat, sof, qof

const NAME = "TP"
const OMEGA = 2.0
const minl = 3 # 最小サイクルの長さ
const NL = 3
const NV = 6 # vertexの数
const NE = 9 # edgeの数
const RANK = 4 # gauge fixingで残るUの数
const poly = [4,4,4,3] # funcamental cycleの角数
const e_free = [1,2,3,7] # gauge fixingで残すedge

# sourceとtarget
S = [1,2,3,1,2,3,4,5,6]
T = [2,3,1,4,5,6,5,6,4]

# (q,u)からsへの変換
sof(q,u=0.0) = -log(q)/log((1-u)*(OMEGA-u))
# (s,u)からqへの変換
qof(s,u=0.0) = ( (1-u)*(OMEGA-u) )^(-s) 
############################################################################
# Path Ihara Zeta の matrix V
function Vmat(q)
    V = zeros(ComplexF64,8,8)
    ###
    V[1,1] = q^5
    V[1,2] = q^1
    V[1,3] = q^4
    V[1,4] = q^4
    V[1,6] = q^4
    V[1,7] = q^5
    V[1,8] = q^2
    ###
    V[2,1] = q^4
    V[2,2] = q^4
    V[2,3] = q^1
    V[2,4] = q^3
    V[2,5] = q^4
    V[2,7] = q^4
    V[2,8] = q^3
    ###
    V[3,1] = q^1
    V[3,2] = q^5
    V[3,3] = q^4
    V[3,4] = q^2
    V[3,5] = q^5
    V[3,6] = q^4
    V[3,8] = q^4
    ###
    V[4,1] = q^4
    V[4,2] = q^2
    V[4,3] = q^3
    V[4,4] = q^3
    V[4,5] = q^2
    V[4,6] = q^3
    V[4,7] = q^4
    ###
    V[5,2] = q^5
    V[5,3] = q^4
    V[5,4] = q^2
    V[5,5] = q^5
    V[5,6] = q^4
    V[5,7] = q^1
    V[5,8] = q^4
    ###
    V[6,1] = q^5
    V[6,3] = q^4
    V[6,4] = q^4
    V[6,5] = q^1
    V[6,6] = q^4
    V[6,7] = q^5
    V[6,8] = q^2
    ###
    V[7,1] = q^4
    V[7,2] = q^4
    V[7,4] = q^3
    V[7,5] = q^4
    V[7,6] = q^1
    V[7,7] = q^4
    V[7,8] = q^3
    ###
    V[8,1] = q^2
    V[8,2] = q^4
    V[8,3] = q^3
    V[8,5] = q^4
    V[8,6] = q^3
    V[8,7] = q^2
    V[8,8] = q^3
    ###
    return V
end

end


