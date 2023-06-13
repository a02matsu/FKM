#################################
## Triangle-Square (K4)
module Graph
# 外から見える変数
export NAME, OMEGA, minl, NV, NE, RANK, e_free, S, T, poly, NL
# 外から見える関数
export sof, qof, Vmat

const NAME = "K4"
const OMEGA = 2.0
const minl = 3
const NL = 4
const NV = 4
const NE = 6
const RANK = 3 # gauge fixingで残るUの数。
const poly = [3,3,3] # funcamental cycleの角数
const e_free = [4,5,6] # gauge fixingで残すedge

# sourceとtarget
S = [1,1,1,2,3,4]
T = [2,3,4,3,4,2]

# (q,u)からsへの変換
sof(q,u=0.0) = -log(q)/log((1-u)*(OMEGA-u))
# (s,u)からqへの変換
qof(s,u=0.0) = ( (1-u)*(OMEGA-u) )^(-s) 
############################################################################
# Path Ihara Zeta の matrix V
function Vmat(q)
    V = zeros(ComplexF64,6,6)
    V .= q^3
    for i in 1:3
      V[i,i+3] = 0.0
      V[i+3,i] = 0.0
    end
    ###
    V[1,2] = q
    ###
    V[2,3] = q
    ###
    V[3,1] = q
    ###
    V[4,6] = q
    ###
    V[5,4] = q
    ###
    V[6,5] = q
    return V
end

end


