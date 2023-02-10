using LinearAlgebra # Just to get identity matrix
using Plots
plotly() # Interactive plots

#Subtask a
function propagationM(n,D,k)
    P=[exp(im*n*D*k) 0;
       0 exp(-im*n*D*k)] #NOTE THE MINUS SIGN! I had to derive this from the start. Can easily overlook it.
    return P
end

#Subtask b
function interfaceM(nᵢ₋₁,nᵢ)
    nr=nᵢ/nᵢ₋₁
    I=1/2*[1+nr   1-nr;
           1-nr   1+nr]
    return I
    # See some of the derivation for this 
    # matrix on paper (I almost derived it)
end
# display(propagationM(1.5,1,1))
# display(interfaceM(1,2))



#Subtask c

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

#Test transfer matrix:
# k₀_val=1
# T=transferM(n,D,k₀_val)
# display(T)
# Trans=abs(T[1,1]-T[1,2]*T[2,1]/T[2,2])^2
# Refl=abs(T[2,1]/T[2,2])^2

#Subtask d
function RT_Cal_Plot(n,D)
    N=1001
    j=1
    TransM=Array{Float64}(undef,N,1)
    ReflM=Array{Float64}(undef,N,1)
    lossM=Array{Float64}(undef,N,1)
    k₀M=LinRange(0, 3, N)
    for k₀ in k₀M
        T=transferM(n,D,k₀)
        # display(T)

        TransM[j]=abs(T[1,1]-T[1,2]*T[2,1]/T[2,2])^2
        ReflM[j]=abs(T[2,1]/T[2,2])^2
        lossM[j]=1-TransM[j]-ReflM[j]
        j+=1
    end


    Plots.plot(k₀M,lossM, labels="loss")
    Plots.plot!(title = "Relative energy loss of stack of dielectric layers")
    Plots.plot!(legend=:right)
    Plots.xlabel!("k₀ [m⁻¹]")
    loss_plot=Plots.ylabel!("Relative energy loss [ ]")
    display(loss_plot)

    Plots.plot(k₀M,TransM, labels="T")
    Plots.plot!(k₀M,ReflM, labels="R")
    Plots.plot!(title = "Reflectance and Transmittance of stack of dielectric layers")
    Plots.plot!(legend=:right)
    Plots.xlabel!("k₀ [m⁻¹]")
    TR_plot=Plots.ylabel!("Part of energy [ ]")
    display(TR_plot)
end

## 

n1=1
n2=2
# n=transpose([1 1.2 1.4 1 1.8 2 1.8 1 1.4 1.2 1]) # To play around
n=transpose([1 n2 n1 n2 n1 n2 n1 n2 n1 n2 1])
# D=transpose([1 2 3 0.5 1 0.5 3 2 1]) # To play around
D=transpose([1 1 1 1 1 1 1 1 1])
RT_Cal_Plot(n,D)

##

Fourpair=transpose([n2 n1 n2 n1 n2 n1 n2 n1])
n=vcat(1, Fourpair, Fourpair, n2, n1, n2 , 1)
D=ones(size(n,1)-2)
RT_Cal_Plot(n,D)

##

Fourpair=transpose([n2 n1 n2 n1 n2 n1 n2 n1])
n=vcat(1, Fourpair, Fourpair, Fourpair, Fourpair, n2, n1, n2, n1, n2, n1, n2, 1)
D=ones(size(n,1)-2)
RT_Cal_Plot(n,D)

##

#Subtask e
n1=1
n2=2+5*im
n=transpose([1 n2 n1 n2 n1 n2 n1 n2 n1 n2 1])
D=transpose([1 1 1 1 1 1 1 1 1])
RT_Cal_Plot(n,D)

##

n1=1
n2=2+10*im
n=transpose([1 n2 n1 n2 n1 n2 n1 n2 n1 n2 1])
D=transpose([1 1 1 1 1 1 1 1 1])
RT_Cal_Plot(n,D)

##

n1=1
n2=2+15*im
n=transpose([1 n2 n1 n2 n1 n2 n1 n2 n1 n2 1])
D=transpose([1 1 1 1 1 1 1 1 1])
RT_Cal_Plot(n,D)


