## change log
# ver.01 
# ・グラフの情報をGraph.jlと分離した
# ・行列サイズが大きくなる時の不安定性を開所するため、actionを計算する際、LU分解をして、Uの対角成分のlogを足し算することにした。
# ver.02
# ・Graphデータを単純なincludeからmoduleに変えた
# ver.03 （2023/06/12）
# ・returnする量をUのconfigurationのみに変更。

module FKM
using Random
using LinearAlgebra
using SparseArrays
using DelimitedFiles
using LsqFit
using Plots
using Printf
using JLD

# 外から見える変数
export NV, NE, RANK, minl, OMEGA, NAME
# 外から見える関数
export 
Metropolis, 
HMC,
backup_file,
check_B,
check_gauge_inv, 
Phases, 
check_force_balance,
LFreverse, 
action, 
unitary,
sof,
qof

############################################################################
## 用いるグラフに合わせて Graph.jlを編集すること
include("Graph.jl")
using .Graph
############################################################################
# function that returns \theta^a T^a
# N : Integer
# theta :: N^2 個のパラメータ
function hermitian(N,theta)
    T=zeros(ComplexF64,N,N)
    a=0
    for i=1:N
        for j=i+1:N
            a+=1
            T[i,j]+=theta[a]/sqrt(2)
            T[j,i]+=theta[a]/sqrt(2)
            #######
            a+=1
            T[i,j] += -im*theta[a]/sqrt(2)
            T[j,i] += im*theta[a]/sqrt(2)
        end
    end
    for i=1:N-1
        a+=1
        for j=1:i
            T[j,j] += theta[a]/sqrt(i^2+i)
        end
        T[i+1,i+1] += -theta[a]*i/sqrt(i^2+i)
    end
    a+=1
    for i=1:N
        T[i,i] += theta[a]/sqrt(N)
    end
    return T
end

#################


############################################################################
# function that returns exp(i \theta^a T^a)
# N : Integer
# theta :: N^2 個のパラメータ
function unitary(N,theta)
    T=zeros(ComplexF64,N,N)
    a=0
    for i=1:N
        for j=i+1:N
            a+=1
            T[i,j]+=theta[a]/sqrt(2)
            T[j,i]+=theta[a]/sqrt(2)
            #######
            a+=1
            T[i,j] += -im*theta[a]/sqrt(2)
            T[j,i] += im*theta[a]/sqrt(2)
        end
    end
    for i=1:N-1
        a+=1
        for j=1:i
            T[j,j] += theta[a]/sqrt(i^2+i)
        end
        T[i+1,i+1] += -theta[a]*i/sqrt(i^2+i)
    end
    a+=1
    for i=1:N
        T[i,i] += theta[a]/sqrt(N)
    end
    return exp(im .* T)
end

############################################################################
# function that returns exp(i \theta^a T^a) with traceless Ta
# N : Integer
# theta :: N^2-1 個のパラメータ
function special_unitary(N,theta)
    T=zeros(ComplexF64,N,N)
    a=0
    for i=1:N
        for j=i+1:N
            a+=1
            T[i,j]+=theta[a]/sqrt(2)
            T[j,i]+=theta[a]/sqrt(2)
            #######
            a+=1
            T[i,j] += -im*theta[a]/sqrt(2)
            T[j,i] += im*theta[a]/sqrt(2)
        end
    end
    for i=1:N-1
        a+=1
        for j=1:i
            T[j,j] += theta[a]/sqrt(i^2+i)
        end
        T[i+1,i+1] += -theta[a]*i/sqrt(i^2+i)
    end
    #a+=1
    #for i=1:N
        #T[i,i] += theta[a]/sqrt(N)
    #end
    return exp(im .* T)
end


############################################################################
# function that returns matrix W(U)
# U : サイズNcのunitary matrixを RANK 個持つ配列
function Wmat(U,Nc)
    W = zeros(ComplexF64, 2*NE*Nc, 2*NE*Nc)
    # unitary 行列を準備
    Ufix = [Array{ComplexF64}(I, Nc, Nc) for i in 1:NE]
    for a in 1:RANK
        Ufix[e_free[a]] .= U[a]
    end
        
    for e=1:NE 
        # set (1,1)-block
        for ee in 1:NE
            if T[e] == S[ee]
                W[(e-1)*Nc+1:e*Nc, (ee-1)*Nc+1:ee*Nc] .= Ufix[e]
            end
        end
        # set (1,2)-block
        for ee in 1:NE
            if T[e] == T[ee] && e != ee
                W[(e-1)*Nc+1:e*Nc, (ee+NE-1)*Nc+1:(ee+NE)*Nc] .= Ufix[e]
            end
        end
    end
    for e in 1:NE
        # set (2,1)-block
        for ee in 1:NE
            if S[e] == S[ee] && e != ee
                W[(e+NE-1)*Nc+1:(e+NE)*Nc, (ee-1)*Nc+1:ee*Nc] .= adjoint(Ufix[e])
            end
        end
        # set (2,2)-block
        for ee in 1:NE
            if S[e] == T[ee]
                W[(e+NE-1)*Nc+1:(e+NE)*Nc, (ee+NE-1)*Nc+1:(ee+NE)*Nc] .= adjoint(Ufix[e])
            end
        end
    end
    return W
end

############################################################################
# function that returns J(U)
# U : サイズNcのunitary matrixを RANK 個持つ配列
function Jmat(U,Nc)
    J = zeros(ComplexF64, 2*NE*Nc, 2*NE*Nc)
    # unitary 行列を準備
    Ufix = [Array{ComplexF64}(I, Nc, Nc) for i in 1:NE]
    for a in 1:RANK
        Ufix[e_free[a]] .= U[a]
    end
    # set (1,2)-block
    for e in 1:NE 
        J[(e-1)*Nc+1:e*Nc, (e+NE-1)*Nc+1:(e+NE)*Nc] .= Ufix[e]
    end
    # set (2,1)-block
    for e in 1:NE
        J[(e+NE-1)*Nc+1:(e+NE)*Nc, (e-1)*Nc+1:e*Nc] .= adjoint(Ufix[e])
    end
    return J
end

##########################################################
# actionのexpがBartholdi zetaの逆数になっているはず
function check_B(Nc,q)
    u=0.0
    # Ihara zeta
    zeta_inv = (1.0-q^4)*(1.0+q^2-2.0*q^3)*(1.0-q^2-2.0*q^3)

    Uunit = [Array{ComplexF64}(I, Nc, Nc) for i in 1:RANK]
    #mat = Matrix{ComplexF64}(I, 2*NE*Nc, 2*NE*Nc) - q * Wmat(Uunit,Nc) - q*u * Jmat(Uunit,Nc)

    det_inv = exp(action(Uunit,Nc,q,u)/Nc)

    ## matの全成分を表示
    #for i in axes(mat[1])
    #    for j in axes(mat[2])
    #        print(real(mat[i,j]), "\t")
    #    end
    #    println()
    #end

    #det_inv = real(det(mat))

    println("ζ=",zeta_inv^Nc,"\t","e^{S/Nc}=",det_inv,"\t","diff=",zeta_inv^Nc-det_inv)

end

#check_B()

############################################################################
# function that returns  Nc * Tr( log(1-q(W+uJ)) ) 
# U : サイズNcのunitary matrixを RANK 個持つ配列
#function action(U, Nc::Int ,gamma::Float64, q::Float64, u::Float64)
function action(U, Nc::Int, q::Float64, u::Float64)
    A = Matrix{ComplexF64}(I, 2*NE*Nc, 2*NE*Nc) - q .* Wmat(U,Nc) - q*u .* Jmat(U,Nc)
    L, V, p = lu(A) # LU分解
    return Nc * sum(log.(abs.(diag(V)))) # log行列式
    #tmp = real( det(A)^2 )
    #println(isinf(tmp))
    #println(isnan(tmp))
    #if isinf(tmp) || isnan(tmp)
      #return Nc * real( sum(log.(eigvals(A))) )
    #else 
      #return Nc * 0.5 * log( tmp )
    #end
    #A = sparse(A) # 疎行列にする
    #return Nc * real( 0.5 * log( real( det(A)^2 ) ) ) 
end
############################################################################


############################################################################
# unitary行列拡張したVmat をVUとしたとき、
# 1 - VU
function ONE_minus_VUmat(U,Nc,q)
    #VU = zeros(ComplexF64, 2*RANK*Nc, 2*RANK*Nc)
    VU = Matrix{ComplexF64}(I, 2*RANK*Nc, 2*RANK*Nc)
    V = Vmat(q)
    for A in 1:2*RANK
        if A in 1:RANK
            Utmp = copy(U[A])
        elseif A in RANK+1:2*RANK
            Utmp = copy(adjoint(U[A-RANK]))
        end
        for B in 1:2*RANK
            VU[(A-1)*Nc+1:A*Nc, (B-1)*Nc+1:B*Nc] .+= -V[A,B] .* Utmp
        end
    end
    return VU
end 


############################################################################
# function that returns  Nc * Tr( log(1-VU)) ) 
# U : サイズNcのunitary matrixを RANK 個持つ配列
function action(U, Nc::Int, q::Float64)
    A = ONE_minus_VUmat(U,Nc,q)
    L, V, p = lu(A) 
    return Nc * sum(log.(abs.(diag(V))))
    #return Nc * real( 0.5 * log( real( det(ONE_minus_VUmat(U,Nc,q))^2 ) ) ) 
    #println(isinf(tmp))
    #println(isnan(tmp))
    #if isinf(tmp) || isnan(tmp)
      #return Nc * real( sum( log.( eigvals( ONE_minus_VUmat(U,Nc,q) ) ) ) )
    #else
      #return Nc * 0.5 * log( real( tmp ) ) 
    #end
    #return Nc * real( 0.5 * log( real( det(ONE_minus_VUmat(U,Nc,q))^2 ) ) ) 
end

############################################################################
## ゲージ不変性のチェック
## 今、ゲージ固定しているので、それぞれのUはひとつのgについて
##   U -> g U g^{-1}
## のように変換することに注意する。
function check_gauge_inv(Nc,q,u)

    U = [Array{Complex{Float64}}(I, Nc, Nc) for i in 1:RANK]
    #g = [Array{Complex{Float64}}(I, Nc, Nc) for i in 1:RANK]
    g = Array{Complex{Float64}}(I, Nc, Nc)

    for i=1:RANK
        theta = ( rand(Float64, Nc^2) .- 0.5 ) .* 2.0 .* pi
        U[i] .= unitary(Nc,theta)
    end

    action_org = action(U,Nc,q,u)

    theta = ( rand(Float64, Nc^2) .- 0.5 ) .* 2.0 .* pi
    g .= unitary(Nc,theta)

    for i=1:RANK
        U[i] .= g * U[i] * adjoint(g)
    end

    action_gU = action(U,Nc,q,u)

    println("org = ", action_org,"\t", "after = ", action_gU, "\t", "diff = ", action_gU - action_org)
end

############################################################################
##  Phases of the eigenvalues of U[1] and U[2]
function Phases(U)
    #eigenvalues1 = eigvals(U[1])
    #eigenvalues2 = eigvals(U[2])
    angles = []
    for a in 1:RANK
        eigenvalues = eigvals(U[a])
        push!(angles, angle(eigenvalues))
    end
    return angles
    #return map(angle,eigenvalues1), map(angle,eigenvalues2)
end

############################################################################
## Metropolis algorithm でシミュレーション
# NEW ::  0:初期状態から 1:読み込み
function Metropolis(NEW::Int, Nc::Int, gamma::Float64, q::Float64, u::Float64, niter::Int, step_size::Float64, configfile="config.txt")
    naccept = 1 # 初期化
    skip_step = 10 # このステップごとに測定値を格納する
    print_step = 500 # このステップごとに標準出力する
    theta = zeros(Complex{Float64}, Nc^2+Nc) 
    phase1 = []
    phase2 = []
    #M = []
    #M2 = []
    S = []
    #confU = []
    U = [Array{ComplexF64}(I, Nc, Nc) for i in 1:RANK]

    #################################################
    ## ファイルの読み込み
    ## 以前のconfigurationがあればそこから読み取る。
    ## ただし、読み込んだ後はファイル名を変更してバックアップにする
    #filename = "config.txt"
    # ファイルが存在する場合
    if isfile(configfile) && NEW != 0
        # ファイルから読み込む
        line = readdlm(configfile, ',')
        # 最初のエントリーは iteration number
        iter_init = line[1] + 1
        # 続いて U[r]
        for r=1:RANK
            for i=1:Nc
                for j=1:Nc
                    real = line[2*Nc*(i-1)+2*j]
                    imag = line[2*Nc*(i-1)+2*j+1]
                    U[r][i,j]= real + im * imag
                end
            end
        end
        # 以前のファイルをバックアップ
        backup_idx = 1
        while isfile("$(configfile).bak$(backup_idx)")
            backup_idx += 1
        end
        # 新しいファイル名を決定
        backup_filename = "$(configfile).bak$(backup_idx)"
        # ファイルをリネーム
        mv(configfile, backup_filename)
    # ファイルは存在するけどNEW==0の場合
    elseif isfile(configfile) && NEW == 0
        iter_init = 1
        # 以前のファイルをバックアップ
        backup_idx = 1
        while isfile("$(configfile).bak$(backup_idx)")
            backup_idx += 1
        end
        # 新しいファイル名を決定
        backup_filename = "$(configfile).bak$(backup_idx)"
        # ファイルをリネーム
        mv(configfile, backup_filename) 
    # ファイルが存在しない、または、NEW==0の場合
    else
        iter_init = 1
    end
    #################################################

    #################################################
    ## MCMC開始
    U_backup = copy(U)
    for iter = 1:niter 
        U_backup .= U
        # 初期action
        action_init = action(U,Nc,q,u)
        # Uを更新
        for r=1:RANK
            theta = ( rand(Float64, Nc^2) .- 0.5 ) .* 2.0 .* step_size
            dU = unitary(Nc,theta)
            #dU = special_unitary(Nc,theta)
            # 更新後action
            U[r] = dU * U[r]
        end
        action_fin = action(U,Nc,q,u)

        # メトロポリステスト
        if exp( gamma * (action_init - action_fin) )  > rand()
            naccept += 1
        else
            U .= U_backup
        end

        if iter % print_step == 0
            println()
            println(
                iter_init+iter-1,"\t",
                norm(U[1] * adjoint(U[1]) - Array{ComplexF64}(I, Nc, Nc), 2 ),"\t",
                action_fin - action_init,"\t", 
                naccept/iter)
        end
    
    
        if iter % skip_step == 0
            # 固有値を取得
            eigenvalues1 = eigvals(U[1])
            eigenvalues2 = eigvals(U[2])
            # 偏角を取得
            push!(phase1, map(angle,eigenvalues1) )
            push!(phase2, map(angle,eigenvalues2) )
            # tr(M)とtr(M^2)を取得
            #push!(M, TrM(U,Nc,q,u))
            #push!(M2, TrM2(U,Nc,q,u))
            # action と　action^2 を取得
            push!(S, action(U,Nc,q,u))
            # Uのconfigを保存
            #push!(confU, copy(U)) 
            # copyがないと、、confUにUの参照が複数回追加されているため、全ての要素が同じオブジェクトを参照する。
            # つまり、Uが更新されるたびに、すべての要素が最新のUの状態を参照するようになる。
            # そのため、confUにはUの最終的な状態しか含まれず、niter回の反復で同じ内容が続くことになる。
        end
    end
    #################################################



    #################################################
    # ファイルに書き出し
    # 書き出す
    n = iter_init+niter-1
    open(configfile, "w") do f
        write(f, "($n),   ")
        for r=1:RANK
            for i=1:Nc
                for j=1:Nc 
                    write(f, join(real(U[r][i,j]), ", "))
                    write(f, ",   ")
                    write(f, join(imag(U[r][i,j]), ", "))
                    write(f, ",   ")
                end
            end
        end
    end
    ##################################################

    return phase1, phase2, S
end


############################################################################
# ファイルをバックアップする関数
#function backup_file(filename)
#    if isfile(filename) 
#        backup_idx = 1
#        while isfile("$(filename).bak$(backup_idx)")
#            backup_idx += 1
#        end
#        # バックアップのファイル名を決定
#        backup_filename = "$(filename).bak$(backup_idx)"
#        # ファイルをリネーム
#        mv(filename, backup_filename) 
#    end
#end
############################################################################
# ファイルをバックアップする関数
function backup_file(body, opt)
    filename="$(body).$(opt)"
    if isfile(filename) 
        backup_idx = 1
        #while isfile("$(filename).bak$(backup_idx)")
        while isfile("$(body)_bak$(backup_idx).$(opt)")
            backup_idx += 1
        end
        # バックアップのファイル名を決定
        backup_filename="$(body)_bak$(backup_idx).$(opt)"
        # ファイルをリネーム
        cp(filename, backup_filename) 
    end
end

############################################################################
# configuration を読み込む関数
function read_config(Nc::Int,configfile)
    U = [Array{ComplexF64}(I, Nc, Nc) for i in 1:RANK]
    # ファイルから読み込む
    line = readdlm(configfile, ',')
    # 最初のエントリーは iteration number
    iter_init = line[1] + 1
    # 続いて U[r]
    for r=1:RANK
        for i=1:Nc
            for j=1:Nc
                real = line[2*Nc*(i-1)+2*j]
                imag = line[2*Nc*(i-1)+2*j+1]
                U[r][i,j]= real + im * imag
            end
        end
    end
    return iter_init, U
end

############################################################################
# configuration を書き出す関数
function write_config(configfile, iter_init, niter, U, Nc::Int)
    n = iter_init+niter-1
    open(configfile, "w") do f
        write(f, "$n,   ")
        for r=1:RANK
            for i=1:Nc
                for j=1:Nc 
                    write(f, join(real(U[r][i,j]), ", "))
                    write(f, ",   ")
                    write(f, join(imag(U[r][i,j]), ", "))
                    write(f, ",   ")
                end
            end
        end
    end
end

#################################################
## Hamiltonian
function hamiltonian(U,P,Nc::Int,gamma::Float64,q::Float64, u::Float64)
    H = 0.0
    for a=1:RANK
        H += 0.5 .* real(tr(P[a] * P[a]))
    end
    H +=  gamma * action(U,Nc,q,u)
    return H
end
#################################################
## Hamiltonian 
## u=0限定、path表示
function hamiltonian(U,P,Nc::Int,gamma::Float64,q::Float64)
    H = 0.0
    for a=1:RANK
        H += 0.5 .* real(tr(P[a] * P[a]))
    end
    H +=  gamma * action(U, Nc, q)
    return H
end

#################################################
## Force 
function force(U, Nc::Int,gamma::Float64, q::Float64, u::Float64)
    F = []
    Minv = inv( Matrix{ComplexF64}(I, 2*NE*Nc, 2*NE*Nc) - q .* Wmat(U,Nc) -  q*u .* Jmat(U,Nc) )
    MAT = zeros(Complex{Float64},Nc,Nc)

    for a in 1:RANK
        tmp = zeros(Complex{Float64}, Nc,Nc)
        for ee in 1:2*NE
            if ee <= NE 
                point = S[ee]
            else
                point = T[ee-NE]
            end
            if T[e_free[a]] == point
                MAT = copy(Minv[(ee-1)*Nc+1:ee*Nc, (e_free[a]-1)*Nc+1:e_free[a]*Nc])
                if ee == e_free[a]+NE 
                    tmp += u .* U[a] * MAT
                else
                    tmp += U[a] * MAT
                end
            end
            ##
            if S[e_free[a]] == point
                MAT = copy(Minv[(ee-1)*Nc+1:ee*Nc, (NE+e_free[a]-1)*Nc+1:(NE+e_free[a])*Nc])
                if ee == e_free[a]
                    tmp += -u .* MAT * adjoint(U[a])
                else
                    tmp += - MAT * adjoint(U[a])
                end
            end
        end
        # 係数の調整
        tmp .= im * gamma * q * Nc .* tmp 
        # forceがhermitianになるようにproject
        tmp = 0.5 .* (copy(tmp) + adjoint(copy(tmp)))
        # 格納
        push!( F, tmp )
    end
    return F
end
#################################################
## Force for path zeta function
function force(U, Nc::Int, gamma::Float64, q::Float64)
    F = []
    V = Vmat(q)
    Minv = inv( ONE_minus_VUmat(U,Nc,q) )

    for a in 1:RANK
        Fa = zeros(Complex{Float64}, Nc, Nc)
        tmp = zeros(Complex{Float64}, Nc ,Nc)
        for B in 1:2*RANK
            tmp += V[a,B] .* Minv[(B-1)*Nc+1:B*Nc,(a-1)*Nc+1:a*Nc]
        end
        Fa += U[a] * tmp
        tmp = zeros(Complex{Float64}, Nc, Nc)
        for B in 1:2*RANK
            tmp += V[a+RANK,B] .* Minv[(B-1)*Nc+1:B*Nc,(a+RANK-1)*Nc+1:(a+RANK)*Nc]
        end
        Fa .+= -tmp * adjoint(U[a])

        Fa .*= im * gamma * Nc 
        Fa = 0.5 .* ( copy(Fa) + adjoint(copy(Fa)) ) 

        push!( F, Fa )
    end
    return F
end

#################################################
## molecular dynamics
function molecular_evolution(U,P,Nc::Int,gamma::Float64,q::Float64,u::Float64,Ntau::Int,Dtau::Float64)
    # first step
    for a in 1:RANK
        U[a] = exp( im * 0.5 * Dtau .* P[a] ) * copy(U[a] )
    end
    # intermediate steps
    for i in 1:Ntau-1
        F = force(U,Nc,gamma,q,u)
        P += Dtau .* F
        for a=1:RANK
            U[a] = exp( im * Dtau .* P[a] ) * copy(U[a])
        end
    end
    # final step
    #println("force = ",norm(F[1]))
    F = force(U,Nc,gamma,q,u)
    P += Dtau .* F
    for a in 1:RANK
        U[a] = exp( im * 0.5 * Dtau .* P[a] ) * copy(U[a])
    end


    for a in 1:RANK
        if( check_unitary(U[a],Nc) > 1e-13 )
            V = copy(U[a])
            U[a] = closest_unitary_matrix(V)
        end
    end
    return U,P
end

#################################################
## molecular dynamics
## u=0 限定、伊原ゼータ関数のpath表示利用
function molecular_evolution(U,P,Nc::Int,gamma::Float64,q::Float64,Ntau::Int,Dtau::Float64)
    # first step
    for a in 1:RANK
        U[a] = exp( im * 0.5 * Dtau .* P[a] ) * copy(U[a])
    end
    # intermediate steps
    for i in 1:Ntau-1
        F = force(U,Nc,gamma,q)
        P += Dtau .* F
        for a=1:RANK
            U[a] = exp( im * Dtau .* P[a] ) * copy(U[a])
        end
    end
    # final step
    F = force(U,Nc,gamma,q)
    P += Dtau .* F
    for a in 1:RANK
        U[a] = exp( im * 0.5 * Dtau .* P[a] ) * copy(U[a])
    end


    for a in 1:RANK
        if( check_unitary(U[a],Nc) > 1e-13 )
            V = copy(U[a])
            U[a] = closest_unitary_matrix(V)
        end
    end
    return U,P
end

#################################################
## Leap Flogの反転チェック
function LFreverse(Nc::Int,gamma::Float64,q::Float64,u::Float64,Ntau::Int,Dtau::Float64)
    U = [unitary( Nc, 1.0.*(rand(Nc^2) .- 0.5).*pi ) for i in 1:RANK] 
    P = [hermitian( Nc, randn(Nc^2) ) for i in 1:RANK] 
    Uback = copy(U)
    Pback = copy(P)
    U2,P2 = molecular_evolution(copy(U),copy(P),Nc,gamma,q,u,Ntau,Dtau)

    U3,P3 = molecular_evolution(copy(U2),copy(P2).*(-1.0),Nc,gamma,q,u,Ntau,Dtau)

    for a in 1:RANK
       println("check = ",norm(U2[a] - U3[a]))
       println("dist = ",norm(U3[a] - Uback[a]))
    end
end


#################################################
## unitary の程度を判定
function check_unitary(U, Nc)
    return norm( Array{ComplexF64}(I, Nc, Nc) - copy(U) * inv(copy(U)) )
end

#################################################
## unitary行列にproject
function closest_unitary_matrix(U::Matrix{ComplexF64})
    # Step 1: Compute the SVD of U
    F = svd(U)
    # Step 2: Replace S with the identity matrix
    I = Diagonal(ones(size(U, 1)))
    # Step 3: Compute the new unitary matrix V
    V = F.U * I * F.Vt
    return V
end

######################################################
## HMC でのシミュレーション
function HMC(NEW::Int, Nc::Int, gamma::Float64, q::Float64, u::Float64, niter::Int, step_size::Float64, Ntau::Int,configbody="config")
    naccept = 1 # 初期化
    skip_step = 10 # このステップごとに測定値を格納する
    print_step = 10000 # このステップごとに標準出力する
    Dtau = step_size / Ntau
    Uconf = [] # Uのhistory
    #phases = [[] for _ in 1:RANK]
    #S = []
    # Hot start
    #U = [unitary(Nc, 2.0 .* (rand(Float64,Nc^2) .- 0.5) .* pi) for i in 1:RANK]
    # Cold start
    U = [Array{ComplexF64}(I, Nc, Nc) for i in 1:RANK]
    #U = [unitary(Nc, 2.0 .* (rand(Float64,Nc^2) .- 0.5) .* 10^(-6)) for i in 1:RANK]

    #################################################
    ## ファイルの読み込み
    ## 以前のconfigurationがあればそこから読み取る。
    ## ただし、読み込んだ後はファイル名を変更してバックアップにする
    #filename = "config.txt"
    # ファイルが存在する場合
    #configfile = "$(configbody).txt"
    configfile = "$(configbody).jld"
    if isfile(configfile) && NEW != 0
        iter_init = load(configfile, "iter")
        iter_init += 1
        U = load(configfile, "Uconf")
        backup_file(configbody,"jld")
        #iter_init, U = read_config(Nc, configfile)
        #backup_file(configbody,"txt")
    # ファイルは存在するけどNEW==0の場合
    elseif isfile(configfile) && NEW == 0
        iter_init = 1
        # 以前のファイルをバックアップ
        backup_file(configbody,"jld")
        #backup_file(configbody,"txt")
    # ファイルが存在しない、または、NEW==0の場合
    else
        iter_init = 1
    end
    #################################################

    #################################################
    ## MCMC開始
    U_backup = copy(U)
    for iter = 1:niter 
        # バックアップ
        U_backup .= U
        # momentumを乱数で指定
        P = [hermitian( Nc, randn(Nc^2) ) for i in 1:RANK]
        # 初期Hamiltonian
        H_init = hamiltonian(U,P,Nc,gamma,q,u)
        # UとPの更新
        U, P = molecular_evolution(copy(U),copy(P),Nc::Int,gamma::Float64,q::Float64,u::Float64,Ntau,Dtau)
        # 更新後のHamiltonian
        H_fin = hamiltonian(U,P,Nc,gamma,q,u)
        # メトロポリステスト
        if exp( H_init - H_fin )  > rand()
            naccept += 1
        else
            U .= U_backup
        end

        # 結果の書き出し
        if iter % print_step == 0
            println()
            println(
                "q:",@sprintf("%.3E",q),"\t",
                "Ntau:",Ntau,"\t",
                "iter:",iter_init+iter-1,"\t",
                "|U-1|",@sprintf("%.2E", norm(U[1] * adjoint(U[1]) - Array{ComplexF64}(I, Nc, Nc), 2 )) ,"\t",
                "dH:",@sprintf("%.2E",abs(H_fin - H_init)),"\t", 
                "acc:",@sprintf("%.3f",naccept/iter))
        end
        # 物理量の書き出し
        if iter % skip_step == 0
            # 偏角を取得
            #for a in 1:RANK
                #push!(phases[a], map(angle,eigvals(U[a])))
            #end
            push!(Uconf, copy(U))
            #push!(S, action(U,Nc,q,u))
            # Uのconfigを保存
            #push!(confU, copy(U)) 
            # copyがないと、、confUにUの参照が複数回追加されているため、全ての要素が同じオブジェクトを参照する。
            # つまり、Uが更新されるたびに、すべての要素が最新のUの状態を参照するようになる。
            # そのため、confUにはUの最終的な状態しか含まれず、niter回の反復で同じ内容が続くことになる。
        end
    end

    #################################################
    # 最終的なconfigurationをファイルに書き出し
    #write_config(configfile, iter_init, niter, U, Nc)
    jldopen(configfile, "w"; compress = true) do f 
      f["Uconf"] = U
      f["iter"] = iter_init+niter-1
      #f["large_array"] = zeros(10000)
      #save("$(config_file).jld", "Uconf", Uconf)
    end

    return naccept/niter, Uconf
end

######################################################
## HMC でのシミュレーション
##  u=0 に限定した、path表示を利用したHMC
function HMC(NEW::Int, Nc::Int, gamma::Float64, q::Float64, niter::Int, step_size::Float64, Ntau::Int,configbody="config")
    naccept = 1 # 初期化
    skip_step = 10 # このステップごとに測定値を格納する
    print_step = 10000 # このステップごとに標準出力する
    Dtau = step_size / Ntau
    #phases = [[] for _ in 1:RANK] # Uの固有値のhistory
    Uconf = [] # Uのhistory
    #S = []
    # Hot start
    #U = [unitary(Nc, 2.0 .* (rand(Float64,Nc^2) .- 0.5) .* pi) for i in 1:RANK]
    # Cold start
    U = [Array{ComplexF64}(I, Nc, Nc) for _ in 1:RANK]

    #################################################
    ## ファイルの読み込み
    ## 以前のconfigurationがあればそこから読み取る。
    ## ただし、読み込んだ後はファイル名を変更してバックアップにする
    #filename = "config.txt"
    # ファイルが存在する場合
    #configfile = "$(configbody).txt"
    configfile = "$(configbody).jld"
    if isfile(configfile) && NEW != 0
        iter_init = load(configfile, "iter")
        iter_init += 1
        U = load(configfile, "Uconf")
        backup_file(configbody,"jld")
        #iter_init, U = read_config(Nc, configfile)
        #backup_file(configbody,"txt")
    # ファイルは存在するけどNEW==0の場合
    elseif isfile(configfile) && NEW == 0
        iter_init = 1
        # 以前のファイルをバックアップ
        backup_file(configbody,"jld")
        #backup_file(configbody,"txt")
    # ファイルが存在しない、または、NEW==0の場合
    else
        iter_init = 1
    end
    #################################################

    #################################################
    ## MCMC開始
    U_backup = copy(U)
    for iter = 1:niter 
        # バックアップ
        U_backup .= U
        # momentumを乱数で指定
        P = [hermitian( Nc, randn(Nc^2) ) for i in 1:RANK]
        # 初期Hamiltonian
        H_init = hamiltonian(U,P,Nc,gamma,q)
        # UとPの更新
        U, P = molecular_evolution(copy(U),copy(P),Nc::Int,gamma::Float64,q::Float64,Ntau,Dtau)
        # 更新後のHamiltonian
        H_fin = hamiltonian(U,P,Nc,gamma,q)
        # メトロポリステスト
        if exp( H_init - H_fin )  > rand()
            naccept += 1
        else
            U .= U_backup
        end

        # 結果の書き出し
        if iter % print_step == 0
            println()
            println(
                "q:",@sprintf("%.3E",q),"\t",
                "Ntau:",Ntau,"\t",
                "iter:",iter_init+iter-1,"\t",
                "|U-1|",@sprintf("%.2E", norm(U[1] * adjoint(U[1]) - Array{ComplexF64}(I, Nc, Nc), 2 )) ,"\t",
                "dH:",@sprintf("%.2E",abs(H_fin - H_init)),"\t", 
                "acc:",@sprintf("%.3f",naccept/iter))
        end
        # 物理量の書き出し
        if iter % skip_step == 0
            # 偏角を取得
            # Uのhistoryを取得
            #for a in 1:RANK
                #push!(phases[a], map(angle,eigvals(U[a])))
            #end
            push!(Uconf, copy(U))
            # action と　action^2 を取得
            #push!(S, action(U,Nc,q))
            # Uのconfigを保存
            #push!(confU, copy(U)) 
            # copyがないと、、confUにUの参照が複数回追加されているため、全ての要素が同じオブジェクトを参照する。
            # つまり、Uが更新されるたびに、すべての要素が最新のUの状態を参照するようになる。
            # そのため、confUにはUの最終的な状態しか含まれず、niter回の反復で同じ内容が続くことになる。
        end
    end

    #################################################
    # 最終的なconfigurationをファイルに書き出し
    #write_config(configfile, iter_init, niter, U, Nc)
    jldopen(configfile, "w"; compress = true) do f 
      f["Uconf"] = U
      f["iter"] = iter_init+niter-1
      #f["large_array"] = zeros(10000)
      #save("$(config_file).jld", "Uconf", Uconf)
    end

    #return naccept/niter, phases, S
    return naccept/niter, Uconf
end

#################################################
# ハミルトニアンの収束をチェック
function check_force_balance(Nc::Int, gamma::Float64, q::Float64, u::Float64, tau::Float64, Nmax::Int)
    # 初期配位を乱数で決定
    U = [unitary(Nc, 2.0 .* (rand(Float64,Nc^2) .- 0.5) .* pi) for i in 1:RANK]
    P = [hermitian( Nc, randn(Nc^2) ) for i in 1:RANK]
    # バックアップ
    U_backup = copy(U)
    P_backup = copy(P)

    N = []
    dH = []
    #tau = 0.10
    Ntau = 10
    for m in 1:Nmax
        Ntau = Int(round(Ntau * 2))
        Dtau = tau / Ntau
        # 初期ハミルトニアン
        H_init = hamiltonian(U,P,Nc,gamma,q,u)
        # UとPの更新
        U, P = molecular_evolution(copy(U),copy(P),Nc::Int,gamma::Float64,q::Float64,u::Float64,Ntau,Dtau)
        # 更新後のHamiltonian
        H_fin = hamiltonian(U,P,Nc,gamma,q,u)
        # メトロポリステスト
        push!(N, Ntau)
        push!(dH, abs(H_init - H_fin) )
        #println(N,dH)
        println(Ntau,"\t",  abs(H_init - H_fin) )
        U .= U_backup
        P .= P_backup
    end 
    @. model(x, p) = p[1]*x^(-2.0)
    p0 = [0.001]
    fit = curve_fit(model, N, dH, p0)

    println("coef = ",coef(fit))
    println("err = ",stderror(fit))


    scatter(N, dH, xscale=:log10, yscale=:log10, label="\$\\Delta H\$")
    plot!(N, model(N,coef(fit)), label="fitting \$N^{-2}\$",markersize=10)
    xlabel!("Ntau")
    ylabel!("\$\\Delta H\$")
end

#################################################
# ハミルトニアンの収束をチェック
function check_force_balance(Nc::Int, gamma::Float64, q::Float64, tau::Float64, Nmax::Int)
    # 初期配位を乱数で決定
    U = [unitary(Nc, 2.0 .* (rand(Float64,Nc^2) .- 0.5) .* pi) for i in 1:RANK]
    P = [hermitian( Nc, randn(Nc^2) ) for i in 1:RANK]
    # バックアップ
    U_backup = copy(U)
    P_backup = copy(P)

    N = []
    dH = []
    #tau = 0.10
    Ntau = 10
    for m in 1:Nmax
        Ntau = Int(round(Ntau * 2))
        Dtau = tau / Ntau
        # 初期ハミルトニアン
        H_init = hamiltonian(U,P,Nc,gamma,q)
        # UとPの更新
        U, P = molecular_evolution(copy(U),copy(P),Nc::Int,gamma::Float64,q::Float64,Ntau,Dtau)
        # 更新後のHamiltonian
        H_fin = hamiltonian(U,P,Nc,gamma,q)
        # メトロポリステスト
        push!(N, Ntau)
        push!(dH, abs(H_init - H_fin) )
        #println(N,dH)
        println(Ntau,"\t",  abs(H_init - H_fin) )
        U .= U_backup
        P .= P_backup
    end 
    @. model(x, p) = p[1]*x^(-2.0)
    p0 = [0.001]
    fit = curve_fit(model, N, dH, p0)

    println("coef = ",coef(fit))
    println("err = ",stderror(fit))


    scatter(N, dH, xscale=:log10, yscale=:log10, label="\$\\Delta H\$")
    plot!(N, model(N,coef(fit)), label="fitting \$N^{-2}\$",markersize=10)
    xlabel!("Ntau")
    ylabel!("\$\\Delta H\$")
end

end


