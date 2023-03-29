
using LinearAlgebra 
using FFTW
using Plots
# using PlotlyJS
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


## Part b

# Initializing space and time dimension
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


#Initial condition
κ=1/2
ω=1/2
ψ=@. exp(-κ*x^2)
# ψ=@. exp(-1/4*x^2)
# ψ=@. exp(-κ*(x-4)^2)
# ψ =@. exp(-2(x-4)^2) - exp(-2(x + 4)^2)


function ℋψ(ψ) # Note that this function returns the Hamiltonian ACTING ON ψ
    # x=range(-10,10,N+1)
    # xs=x[1:end-1]
    xs=x # OBS, I take x-values from a globally defined variable. 
        #Not good practice, but made such that it works in the more general RK4

    # Calculate second order derivative in space
    ψf=fft(ψ)
    A=im*k0*[range(0,N/2-1,Int(N/2)); range(-N/2,-1,Int(N/2))]
    A=A.^2 #The second order derivative
    ψf_ddx=A.*ψf
    ψ_ddx=ifft(ψf_ddx)

    Hψ= 1/2*(-ψ_ddx+xs.^2 .*ψ)
    return Hψ    
end

function dψ(ψ)
    # # x=range(-10,10,N+1)
    # # xs=x[1:end-1]
    # xs=x # OBS, I take x-values from a globally defined variable. 
    #     #Not good practice, but made such that it works in the more general RK4

    # # Calculate second order derivative in space
    # ψf=fft(ψ)
    # A=im*k0*[range(0,N/2-1,Int(N/2)); range(-N/2,-1,Int(N/2))]
    # A=A.^2 #The second order derivative
    # ψf_ddx=A.*ψf
    # ψ_ddx=ifft(ψf_ddx)

    ψ_dt=-im*ℋψ(ψ)
    return ψ_dt    
end



# Calculating ψ over time with Runge Kutta 4
ψ_tot=RK4(dψ,ψ,t)



# Extra material. Calculating energy.
E_tot = Vector{ComplexF64}(undef,Nt)
ψ_current = Vector{ComplexF64}(undef,size(ψ,1))
for i in 1:Nt
    
    ψ_current[:]=ψ_tot[:,i]
    # ψ[:,1] = ψ_tot[:,1]

    E_tot[i]=h*ψ_current'*ℋψ(ψ_current) #OBS: times step size h to get the analytical value
end
E_tot
plot(t, real(E_tot), label="Real(E_tot)")
plot!(t, imag(E_tot), label="Imag(E_tot)")
xlabel!("Normalized time, t [ ]")
ylabel!("Normalized Energy [ ]")
plot!(title = "Energy of wavefunction")


##


ψ_prop = @. log10(abs(ψ_tot)^2 .+10^-7) #10^-7 to better see behaviour
# ψ_prop = @. (abs(ψ_tot)^2)

ψplot = heatmap(t, x, ψ_prop)
xlabel!("Normalized time, t [ ]")
ylabel!("Normalized position, x [ ]")
plot!(title = "Propapility of wavefunction |ψ|²")






##! Part c
κ=1/2
ω=1/2
ψ_anal = exp.(-κ*x.^2)*transpose(exp.(-im*ω*t))

# heatmap(real(ψ_anal))
# heatmap(real(ψ_tot))

ψ_anal_prop = @. log10(abs(ψ_tot)^2 .+10^-7)
# ψ_anal_prop = @. (abs(ψ_tot)^2)
ψplot_anal = heatmap(t, x, ψ_anal_prop)
xlabel!("Normalized time, t [ ]")
ylabel!("Normalized position, x [ ]")
plot!(title = "Analytical propapility of wavefunction |ψ|²")

residual = ψ_tot-ψ_anal
log_abs_residual = @. log10(abs(residual).+10^(-14)) #+10^-14 to not get -∞ in plot

ψplot_residual = heatmap(t, x, log_abs_residual)
xlabel!("Normalized time, t [ ]")
ylabel!("Normalized position, x [ ]")
plot!(title = "Residual for wavefunction |ψ|²")

total_plot=plot(ψplot, ψplot_anal, ψplot_residual, layout=(3,1))
display(total_plot)

function error_cal1(residual,h)
    Nt=size(residual,2)
    err_N = Vector{Float64}(undef,Nt)
    for i in 1:Nt
        err_N[i]=√(sum(h.*(abs.(residual[:,i])).^2)) # L2 norm in space

        # err_N[i]= abs(residual[Int(round(size(residual,1)/2)),i]) #Just a single point
    end
    return err_N
end

function error_cal2(residual,h)
    err_N=√(sum(h.*(abs.(residual)).^2)) # L2 norm in time
    return err_N
end

residual_vs_t = error_cal1(residual, step(x))
residual_vs_t_plot=plot(t, residual_vs_t, label="residual")
xlabel!("Normalized time, t [ ]")
ylabel!("Residual [ ]")
plot!(title = "Total spatial residual vs. time")
display(residual_vs_t_plot)


error_tot=error_cal2(residual_vs_t, step(t))











##

# Ns=[20, 50, 100, 200, 300, 500, 800]
Ns=20:2:200

error_tot_vec = Vector{Float64}(undef,size(Ns,1))
time_step_vec = Vector{Float64}(undef,size(Ns,1))


for (jjj,N) in enumerate(Ns)
# Initializing space and time dimension
xmin=-10
xmax=10
L=xmax-xmin

# N=200

x=range(xmin,xmax,N+1)
x=x[1:end-1]
k0=2*pi/L
h=step(x)

t₀=0
t₁=4*pi
dt=1/N
Nt=Int(round(t₁/dt))
t=range(t₀,t₁,Nt)

time_step_vec[jjj]=step(t)

#Initial condition
κ=1/2
ω=1/2
ψ=@. exp(-κ*x^2)

# NEED this here to update xs in it.
function ℋψ(ψ) # Note that this function returns the Hamiltonian ACTING ON ψ
    # x=range(-10,10,N+1)
    # xs=x[1:end-1]
    xs=x # OBS, I take x-values from a globally defined variable. 
        #Not good practice, but made such that it works in the more general RK4

    # Calculate second order derivative in space
    ψf=fft(ψ)
    A=im*k0*[range(0,N/2-1,Int(N/2)); range(-N/2,-1,Int(N/2))]
    A=A.^2 #The second order derivative
    ψf_ddx=A.*ψf
    ψ_ddx=ifft(ψf_ddx)

    Hψ= 1/2*(-ψ_ddx+xs.^2 .*ψ)
    return Hψ    
end


# Calculating ψ over time with Runge Kutta 4
ψ_tot=RK4(dψ,ψ,t)

# Analytical result
κ=1/2
ω=1/2
ψ_anal = exp.(-κ*x.^2)*transpose(exp.(-im*ω*t))

# Residual and errors
residual = ψ_tot-ψ_anal

residual_vs_t = error_cal1(residual, step(x))

# display(residual_vs_t)
error_tot_vec[jjj]=error_cal2(residual_vs_t, step(t))

end

# display(error_tot_vec)

# plot(time_step_vec,error_tot_vec)
# plot(time_step_vec,error_tot_vec, xaxis=:log, yaxis=:log)
# plot(time_step_vec,error_tot_vec, yaxis=:log)
plot(Ns,error_tot_vec, yaxis=:log)
plot(Ns,error_tot_vec, xaxis=:log, yaxis=:log)
xlabel!("N")
ylabel!("Error")
plot!(title="Convergence plot")



