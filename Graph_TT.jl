################################
## Triangle-Square (TS)
module Graph
# 外から見える変数
export NAME, OMEGA, minl, NV, NE, RANK, e_free, S, T, poly, NL
# 外から見える関数
export Vmat, sof, qof

const NAME = "TT"
const OMEGA = 1.863706528
const minl = 3
const NL = 3
const NV = 5
const NE = 7
const RANK = 3 # gauge fixingで残るUの数。
const poly = [3,3,3] # funcamental cycleの角数
const e_free = [2,3,5] # gauge fixingで残すedge

# cycle graphのsourceとtarget
S = [1,2,3,4,5,3,3]
T = [2,3,4,5,1,1,5]

# (q,u)からsへの変換
sof(q,u=0.0) = -log(q)/log((1-u)*(OMEGA-u))
# (s,u)からqへの変換
qof(s,u=0.0) = ( (1-u)*(OMEGA-u) )^(-s) 
############################################################################
# Path Ihara Zeta の matrix V
function Vmat(q)
    V = zeros(ComplexF64,6,6)
    ###
    V[1,1] = q^3
    V[1,2] = q^1
    V[1,3] = q^2
    V[1,5] = q^3
    V[1,6] = q^2
    ###
    V[2,1] = q^5
    V[2,2] = q^3
    V[2,3] = q^2
    V[2,4] = q^3
    V[2,6] = q^4
    ###
    V[3,1] = q^2
    V[3,2] = q^2
    V[3,3] = q^3
    V[3,4] = q^2
    V[3,5] = q^4
    ###
    V[4,2] = q^3
    V[4,3] = q^4
    V[4,4] = q^3
    V[4,5] = q^5
    V[4,6] = q^2
    ###
    V[5,1] = q^3
    V[5,3] = q^2
    V[5,4] = q^1
    V[5,5] = q^3
    V[5,6] = q^2
    ###
    V[6,1] = q^4
    V[6,2] = q^2
    V[6,4] = q^2
    V[6,5] = q^2
    V[6,6] = q^3
    ###
    return V
end

end


