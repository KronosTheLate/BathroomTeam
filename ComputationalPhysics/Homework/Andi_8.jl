
using LinearAlgebra 
using FFTW
using Plots
plotly() # Interactive plots
# gr()

function RK4(f,x₀,t)
    x=copy(x₀)
    h=step(t)
    xtot = Matrix{Any}(undef,size(x,1),size(t,1)) 
    # xtot=zeros(size(x,1),size(t,1))
    xtot[:,1]=x
    for i in 2:size(t,1)
        k1=f(x)
        k2=f(x+h*k1/2)
        k3=f(x+h*k2/2)
        k4=f(x+h*k3)
        x+=1/6*(k1+2*k2+2*k3+k4)*h
        xtot[:,i]=x
    end
    
    return xtot
end
RK4(args)=RK4(args...)





xmin=-10
xmax=10
L=xmax-xmin

N=200

x=range(xmin,xmax,N+1)
x=x[1:end-1]
k0=2*pi/L
h=step(x)

t₀=0
t₁=4*pi
dt=1/N
Nt=Int(round(t₁/dt))
t=range(t₀,t₁,Nt)



κ=1/2
ω=1/2
ψ=@. exp(-κ*x^2)

function dψ(ψ)
    x=range(-10,10,N+1)
    xs=x[1:end-1]
    
    # Calculate second order derivative in space
    ψf=fft(ψ)
    A=im*k0*[range(0,N/2-1,Int(N/2)); range(-N/2,-1,Int(N/2))]
    A=A.^2 #The second order derivative
    ψf_ddx=A.*ψf
    ψ_ddx=ifft(ψf_ddx)

    ψ_dt=-im/2*(-ψ_ddx+xs.^2*ψ)
    return ψ_dt    
end

ψ_tot=RK4(dψ,ψ,t)


ψ_prop = @.(abs(ψ_tot)^2)


heatmap(ψ_prop)

# ψ_anal=@. exp(-κ*x^2-im*ω*t)

##

xmin=0
xmax=2*pi
L=xmax-xmin

N=100

x=range(xmin,xmax,N+1)
x=x[1:end-1]
k0=2*pi/L
h=step(x)

ψ=@. cos(2*x)


ψf=fft(ψ)
A=im*k0*[range(0,N/2-1,Int(N/2)); range(-N/2,-1,Int(N/2))]
A=A.^2 #The second order derivative
ψf_ddx=A.*ψf
ψ_ddx=ifft(ψf_ddx)

plot(x,ψ)
plot!(x,abs(ψ_ddx)
