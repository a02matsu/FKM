################################
## Triangle-Square (TS)
module Graph
# 外から見える変数
export NAME, OMEGA, minl, NV, NE, RANK, e_free, S, T, poly, NL
# 外から見える関数
export Vmat, sof, qof

const NAME = "TS"
const OMEGA = 1.42405
const minl = 3
const NL = 1
const NV = 5
const NE = 6
const RANK = 2 # gauge fixingで残るUの数。
const poly = [3,4] # funcamental cycleの角数
const e_free = [1,3] # gauge fixingで残すedge

# cycle graphのsourceとtarget
S = [1,2,3,4,5,3]
T = [2,3,4,5,1,1]

# (q,u)からsへの変換
sof(q,u=0.0) = -log(q)/log((1-u)*(OMEGA-u))
# (s,u)からqへの変換
qof(s,u=0.0) = ( (1-u)*(OMEGA-u) )^(-s) 
############################################################################
# Path Ihara Zeta の matrix V
function Vmat(q)
    V = zeros(ComplexF64,4,4)
    ###
    V[1,1] = q^3
    V[1,2] = q
    V[1,4] = q^3
    ###
    V[2,1] = q^4
    V[2,2] = q^4
    V[2,3] = q^4
    ###
    V[3,2] = q^3
    V[3,3] = q^3
    V[3,4] = q^3
    ###
    V[4,1] = q^4
    V[4,3] = q^2
    V[4,4] = q^4
    return V
end

end


