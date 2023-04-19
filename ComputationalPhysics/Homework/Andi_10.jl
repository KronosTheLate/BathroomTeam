
using LinearAlgebra 
using Plots
# using PlotlyJS
# plotly() # Interactive plots
gr()
using Trapz


# Initializing space dimension
xmin=0
xmax=1
L=xmax-xmin

# N=100

function error_cal2(residual,h) #Can maybe use trapz() instead
    err_N=√(sum(h.*(abs.(residual)).^2)) # L2 norm in time, and where error accounts for stepsize
    # err_N=√(h)*norm(residual) # More elegant way to calculate error 
    return err_N
end

# Setting up different parameters
Ns=[10, 20, 50, 100, 200, 500, 1000]
error=zeros(size(Ns,1))
ω=zeros(size(Ns,1))
ω_anal=π

# For loop over several amount of discretization points: Ns
for (jj,N) in enumerate(Ns)

    #Make 1D mesh
    mesh=range(xmin,xmax,N+2)
    h=step(mesh)

    # Initialize and calculate mass and stiffness matrix. See formula on slides
    ℳ = zeros(N,N) # Could optimize by not initilializing and store full matrices, when only tri-diagonal
    𝒮 = zeros(N,N) # Could optimize by not initilializing and store full matrices, when only tri-diagonal
    v=ones(N+1) # NOTE that v is N+1 long
    # N+1 is the number of elements. So one v-value for each element

    for i in 1:N
        ii=Int(i+1) #Mesh offset compared to for-loop index
        ℳ[i,i]=2*(mesh[ii+1]-mesh[ii-1])
        if i!=1
            ℳ[i,i-1]=mesh[ii]-mesh[ii-1]
        end
        if i!=N
            ℳ[i,i+1]=mesh[ii+1]-mesh[ii]
        end
        

        # Explanation of which v to choose (I was not sure of this in the beginning):
            # Note that size(mesh,1)=N+2, and size(v,1)=N+1. So the element between mesh[ii]
            #   and mesh[ii+1] is the v[ii] element.
        𝒮[i,i]=(v[ii]^2/(mesh[ii+1]-mesh[ii]) + v[ii-1]^2/(mesh[ii]-mesh[ii-1]))
        if i!=1
            𝒮[i,i-1]=-v[ii-1]^2/(mesh[ii]-mesh[ii-1])
        end
        if i!=N
            𝒮[i,i+1]=-v[ii]^2/(mesh[ii+1]-mesh[ii])
        end

    end
    ℳ*=1/6

    # display(ℳ)
    # display(𝒮)

    #Solve general eigenvalue problem:
    F = eigen(𝒮, ℳ)
    # display(F.values) # Eigenvalues
    # F.vectors # Eigenvectors
    # plot(F.vectors[:,1]) #Not displayed. But can be used to show eigen-modes
    ω[jj]=sqrt(F.values[1]) # Eigenfrequencies, ω = √(λ)
    

    # Get and normalize numerical calculated eigenmode/eigenshape
    mode=F.vectors[:,1]
    maximum(mode)
    mode=mode/mode[findmax(abs,mode)[2]]
    maximum(mode)

    # Analytical eigenmode (already normalized to max ampltidue being 1)
    mode_anal=sin.(ω_anal*mesh[2:end-1]) #Note this is normalized to 1 (max()=1)
    # Residual
    residual=mode-mode_anal

    # Plotting
    mode_plot=plot(mesh[2:end-1], mode, label="num")
    plot!(mesh[2:end-1], mode_anal, label="anal")
    xlabel!("x in mesh []")
    ylabel!("Normalized mode shape []")
    plot!(title="First mode - shape, N = $N")
    display(mode_plot)

    residual_plot=plot(mesh[2:end-1], residual, label="")
    xlabel!("x in mesh []")
    ylabel!("Residual []")
    plot!(title="Residual of first mode, N = $N")
    display(residual_plot)

    # Write the error for the current N to a vector
    # error[jj]=error_cal2(residual,h)
    error[jj]=√(trapz(mesh[2:end-1],residual.^2))

end


# Calculate error on eigenvalue, and plot together with error on eigenmode calculated in for-loop
eigenvalue_error=ω.-ω_anal

hs=L./Ns
error_plot=plot(hs,error, xaxis=:log, yaxis=:log, label="Mode shape error", markershape=:cross)
plot!(hs,eigenvalue_error, xaxis=:log, yaxis=:log, label="Eigenvalue error", markershape=:cross)
xlabel!("Step size []")
ylabel!("Error []")
plot!(title="Convergence plot")
display(error_plot)







## Subtask d

N=20

# mesh=range(xmin,xmax,N+2)
# h=step(mesh)

mesh1=collect(range(Float64(xmin),√(0.5),Int(N/2)+1))
mesh2=collect(range(√(0.5),Float64(xmax),Int(N/2)+2))[2:end]
mesh=[mesh1; mesh2]
# plot(mesh, markershape=:cross)

v=sqrt.([1*ones(Int(N/2)); 0.1*ones(Int(N/2)+1)]) # NOTE that v is N+1 long
    # N+1 is the number of elements. So one v-value for each element
    # NOTE: not equal amount of elements on each side of the discontinuity.
    #   N=20 leads to 10 elements to the left and 11 elements to the right.
# @show v

ℳ = zeros(N,N)
𝒮 = zeros(N,N)


for i in 1:N
    ii=Int(i+1) #Mesh offset compared to for-loop index
    ℳ[i,i]=2*(mesh[ii+1]-mesh[ii-1])
    if i!=1
        ℳ[i,i-1]=mesh[ii]-mesh[ii-1]
    end
    if i!=N
        ℳ[i,i+1]=mesh[ii+1]-mesh[ii]
    end
    

    # Explanation of which v to choose (I was not sure of this in the beginning):
        # Note that size(mesh,1)=N+2, and size(v,1)=N+1. So the element between mesh[ii]
        #   and mesh[ii+1] is the v[ii] element.
    𝒮[i,i]=(v[ii]^2/(mesh[ii+1]-mesh[ii]) + v[ii-1]^2/(mesh[ii]-mesh[ii-1]))
    if i!=1
        𝒮[i,i-1]=-v[ii-1]^2/(mesh[ii]-mesh[ii-1])
    end
    if i!=N
        𝒮[i,i+1]=-v[ii]^2/(mesh[ii+1]-mesh[ii])
    end

end
ℳ*=1/6

# display(ℳ)
# display(𝒮)

F = eigen(𝒮, ℳ)

ω=√(F.values[1])
@show ω #Print the angular frequency to compare to problem set value: ω≈2.06

mode=F.vectors[:,1]
maximum(mode)
mode=mode/mode[findmax(abs,mode)[2]]
maximum(mode)


mode_plot=plot(mesh[2:end-1], mode, label="num", markershape=:cross)
xlabel!("x in mesh []")
ylabel!("Normalized mode shape []")
plot!(title="First mode - shape, N = $N")
display(mode_plot)



