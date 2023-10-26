include("FKM.jl")
using Random
using LinearAlgebra
using SparseArrays
using DelimitedFiles
using LsqFit
using Plots
using Printf
using .FKM

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
begin
  using .FKM
  using Random

  # 乱数を初期化する関数。毎回違うseedを使うためにtime()を使う。
  function init_rng()
    seed = round(Int, time()) % 6599 
    return MersenneTwister(seed)
  end
  # それぞれのプロセスで乱数を初期化
  init_rng()
  ## theory parameters
  Nc = 4
  u = 0e-1
  gamma = 1024.0
  # gammaを100倍して整数にし、ファイル名に利用
  gamma_int = Int(round(gamma*100))
  Nf = gamma * Nc 
  # シミュレーションのパラメータ
  niter = 50000
  step_size = 0.10
  SSint = Int(round(step_size*100)) 
end
configfile  = "config_test.txt"
configfile2  = "config_test2.txt"

U = [unitary(Nc, 2.0 .* (rand(Float64,Nc^2) .- 0.5) .* 10^(-6)) for i in 1:RANK]

iter_init, U = read_config(Nc, configfile)
println(action(U,Nc,0.5))
write_config(configfile2, iter_init, 10000, U, Nc)
iter_init, U = read_config(Nc, configfile2)
println(action(U,Nc,0.5))