
#JUST TESTING GIT!!!  -ANDI



using LinearAlgebra
using Plots
plotly()

function propagationM(n,D,k)
    P=[exp(im*n*D*k) 0;
       0 -exp(im*n*D*k)] #NOTE THE MINUS SIGN! I had to derive this from the start. Can easily overllok it.
    return P
end

function interfaceM(nᵢ₋₁,nᵢ)
    nr=nᵢ/nᵢ₋₁
    I=1/2*[1+nr   1-nr;
           1-nr   1+nr]
    return I
    # See some of the derivation for this 
    # matrix on paper (I almost derived it)
end

#@show(propagationM(1.5,1,1))
#@show(interfaceM(1,2))


# T can be represented as a product matrices, since each P, I pair 
# describes all propagation in one layer. Then the next 
# propagation is smiply yet another P, I pair. 
# The reason why the matrices are multiplied "backwards" is simpy 
# the definition of this ray/propagation ABCD matrices. 
# output=ABCD*input, where output is the next input for the next layer.


function transferM(n,D,k₀)
    if size(n,1)!=size(D,1)+2
        throw(DimensionMismatch("n must be 2 larger than D"))
    end
    T=LinearAlgebra.I
    for i in 1:size(D,1)
        I=interfaceM(n[i],n[i+1])
        # @show(I)
        P=propagationM(n[i+1],D[i],k₀)
        # @show(P)
        T=P*I*T
        # @show(T)
    end
    I=interfaceM(n[size(D,1)+1],n[size(D,1)+2])
    # @show(I)
    T=I*T
    # @show(T)
    return T
end

k₀_val=1

# N=10

# n=transpose([1 2 1])
# D=transpose([1])

n1=1
n2=2
n=transpose([1 n2 n1 n2 n1 n2 n1 n2 n1 n2 1])
D=transpose([1 1 1 1 1 1 1 1 1])


T=transferM(n,D,k₀_val)
display(T)

Trans=abs(T[1,1]-T[1,2]*T[2,1]/T[2,2])^2
Refl=abs(T[2,1]/T[2,2])^2


n1=1
n2=2
n=transpose([1 n2 n1 n2 n1 n2 n1 n2 n1 n2 1])
D=transpose([1 1 1 1 1 1 1 1 1])

N=1001
j=1
TransM=Matrix{Float64}(undef,N,1)
ReflM=Matrix{Float64}(undef,N,1)
k₀M=LinRange(0, 3, N)
for k₀ in k₀M
    T=transferM(n,D,k₀)
    # display(T)

    TransM[j]=abs(T[1,1]-T[1,2]*T[2,1]/T[2,2])^2
    ReflM[j]=abs(T[2,1]/T[2,2])^2
    j+=1
end


Plots.plot(k₀M,TransM, labels="T")
Plots.plot!(k₀M,ReflM, labels="R")

##
for i in 1:1:10
    display(i)
end
display(1:1:10)
##

# using Plots
# plotly()

# x=0:0.01:5
# y=x.^2
# # @show(y)
# Plots.plot(x,y)

# ## Lotka-Volterra coupled differential equation

# function LotkaVolterra(y::Vector,dt;α=1,β=0.1,γ=2,δ=0.05)
#     🐁=y[1]
#     🐈=y[2]
#     d🐁  = (α*🐁 - β*🐁*🐈)*dt
#     d🐈 = (-γ*🐈 + δ*🐁*🐈)*dt
#     return dy=[d🐁, d🐈]
# end



# dt=0.001
# N=Int.(10/dt) 
# t_tot=dt*N
# t=0:dt:t_tot-dt
# y=zeros(N,2)

# y_init=[50,15]

# y[1,:]=y_init
# for i in 1:(N-1)
#     y[i+1,:]=y[i,:]+LotkaVolterra(y[i,:],dt)
# end
# # @show(y)

# Plots.plot(t,y[:,1], lw = 5, c=:blue, labels="Mice")
# Plots.plot!(t,y[:,2], lw = 5, c=:red, labels="Cats", legend=:topleft)
# Plots.plot!(title = "Lotka Volterra diff.eq.", label = ["Mice" "Cats"])
# Plots.xlabel!("Time [e.g. years]")
# Plots.ylabel!("Population []")
# Plots.plot!(xlims=(0,t_tot),xticks=0:1:t_tot)
# Plots.plot!(ylims=(0,60))
