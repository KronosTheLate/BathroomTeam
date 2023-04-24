
using LinearAlgebra 
using Plots
# plotly() # Interactive plots
gr()

## Part a

# Defining the Green functions, electric field functions, and electric potential functions.

function 𝒢_Φ(x,x′,y,y′)
    return @. -1/(4π)*log((x-x′)^2+(y-y′)^2)
end

function 𝒢_Ex(x,x′,y,y′)
    return @. 1/(4π)*1/((x-x′)^2+(y-y′)^2) * 2*(x-x′)
end

function 𝒢_Ey(x,x′,y,y′)
    return @. 1/(4π)*1/((x-x′)^2+(y-y′)^2) * 2*(y-y′)
end

function Eabsᵢₙₜ(x,x′,y,y′)
    return @. √(𝒢_Ex(x,x′,y,y′)^2+𝒢_Ey(x,x′,y,y′)^2)
end

# Note the usage of a and b in several of the below functions.
a1=1 #a1=1 and b1=0, leads to homogeneous field with E_x = -1 and E_y = 0. 
b1=0
function Φₑₓₜ(x,y)
    return @. a1*x+b1*y
end

function Exₑₓₜ(x,y) # E = -∇Φ
    return @. -a1
end

function Eyₑₓₜ(x,y) # E = -∇Φ
    return @. -b1
end

## Make an NNxNN grid space in specified range; and calculate and plot the different quantities.

NN=200

xs = range(-2,2,NN)
ys = range(-2,2,NN)

x′=1
y′=0
𝒢s = [𝒢_Φ(x, x′, y, y′) for x in xs, y in ys]
𝒢s_plot = heatmap(xs, ys, reverse(rotl90(𝒢s), dims = 1))
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="Internal Electric potential")
display(𝒢s_plot)

Exs = [𝒢_Ex(x,x′,y,y′) for x in xs, y in ys]
Exs_plot = heatmap(xs, ys, reverse(rotl90(Exs), dims = 1), clim=(-5,5)) #clim=(-5,5)
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="Internal E_x, Electric field in x-direction")
display(Exs_plot)

Eys = [𝒢_Ey(x,x′,y,y′) for x in xs, y in ys]
Eys_plot = heatmap(xs, ys, reverse(rotl90(Eys), dims = 1), clim=(-5,5)) #clim=(-5,5)
# scatter!(Eys_plot, [1], [0.5], markersize=2, markercolor=:white)
        # reverse(rotl90(Eys), dims = 1)
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="Internal E_y, Electric field in y-direction")
display(Eys_plot)

Eabss = [Eabsᵢₙₜ(x,x′,y,y′) for x in xs, y in ys]
Eabss_plot = heatmap(xs, ys, reverse(rotl90(Eabss), dims = 1), clim=(0,5)) #, clim=(0,5)
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="Internal abs(E_tot), Absoulte of electric field")
display(Eabss_plot)

Φs = [Φₑₓₜ(x,y) for x in xs, y in ys]
Φs_plot = heatmap(xs, ys, reverse(rotl90(Φs), dims = 1))
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="External electric potential")
display(Φs_plot)


# heatmap([1 2 3;4 5 6;7 8 9]) # To see how to rotate/orient matrix

## Part b

# Create N point charges on the surface of a circle.
N=10
θ=range(0,2π,N+1)
θ=θ[1:end-1] # Remove the double point in the start/end.

rxs = cos.(θ)
rys = sin.(θ)
rs = hcat(rxs,rys)

#Plot charge positions
charge_pos_plot=scatter(rs[:,1],rs[:,2],aspect_ratio = :equal, markershape=:circle, label="Charge positions")
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="Charge locations, N=$N")
display(charge_pos_plot)

## Part c
# Calculate the charge of the point cahrges. MAIN PART OF THE PROBLEM.
b = -Φₑₓₜ(rxs,rys)

M =  [i==j ? 0.0 : 𝒢_Φ(rxs[i], rxs[j], rys[i], rys[j]) for i in eachindex(rxs), j in eachindex(rys)] 
# M =  [i==j ? 0.0 : 1*i+0*j for i in eachindex(rxs), j in eachindex(rys)] 
    # Note the indexes. Just have to loop over 
# display(M)

M = issymmetric(M) ? Symmetric(M) : M

ρ = M\b
charge_pos_plot=scatter(rs[:,1],rs[:,2],aspect_ratio = :equal, markershape=:circle, label="Charges", zcolor=ρ)
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="Charge locations and relative charge, N=$N")
display(charge_pos_plot)


# Calculate the total electric potential
Φₜₒₜ = Φs
for i in 1:size(rxs,1)
    # print(i)
    x′ = rxs[i]
    y′ = rys[i]
    Φₜₒₜ += ρ[i].*[𝒢_Φ(x, x′, y, y′) for x in xs, y in ys]
end

# NOTE the circle marker color shows the relative charge of the point charge.
Φₜₒₜ_plot = heatmap(xs, ys, reverse(rotl90(Φₜₒₜ), dims = 1), aspect_ratio=:equal)
scatter!(Φₜₒₜ_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, label="", zcolor=ρ.*(N/15)) #, markercolor = :white
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="Φₜₒₜ, Total electric potential, N=$N")
display(Φₜₒₜ_plot)


##

# Calculate the total E-field in x- and y-direction as well as the total abs(Eₜₒₜ)
Exs_tot = [Exₑₓₜ(x, y) for x in xs, y in ys]
for i in 1:size(rxs,1)
    # print(i)
    x′ = rxs[i]
    y′ = rys[i]
    Exs_tot += ρ[i].*[𝒢_Ex(x, x′, y, y′) for x in xs, y in ys]
end

# Note the clim, to not have exploding values
Exs_tot_plot = heatmap(xs, ys, reverse(rotl90(Exs_tot), dims = 1), aspect_ratio=:equal, clim=(-3,3)) #clim=(-3,3)
scatter!(Exs_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="E_x, Electric field in x-direction, N=$N")
display(Exs_tot_plot)


Eys_tot = [Eyₑₓₜ(x, y) for x in xs, y in ys]
for i in 1:size(rxs,1)
    # print(i)
    x′ = rxs[i]
    y′ = rys[i]
    Eys_tot += ρ[i].*[𝒢_Ey(x, x′, y, y′) for x in xs, y in ys]
end

# Note the clim, since one point explode in value. 
# Most likely due to the spatial resolution having a point very close to a point charge wheres very close to a pole/(±)infinity value
Eys_tot_plot = heatmap(xs, ys, reverse(rotl90(Eys_tot), dims = 1), aspect_ratio=:equal, clim=(-3,3)) #clim=(-1.5,1.5)
scatter!(Eys_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="E_y, Electric field in y-direction, N=$N")
display(Eys_tot_plot)


E_tot = @. sqrt(Exs_tot^2+Eys_tot^2)

# Note the clim, to not have exploding values
E_tot_plot = heatmap(xs, ys, reverse(rotl90(E_tot), dims = 1), aspect_ratio=:equal,clim=(0,3)) #clim=(-1.5,1.5)
scatter!(E_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
xlabel!("x [ ]")
ylabel!("y [ ]")
plot!(title="abs(E_tot), Absoulte of total electric field, N=$N")
display(E_tot_plot)



## Part d
# The electric field should be zero inside a metal tube/block.
# Therefore I set the analytical result as having |E_tot|=0 for the center point.

# The below code is the above code just modified to only account for the electric 
#  field in the center of the metal bar. 
#  It is also made more efficient by only calculating the E-field in the center. 
#  So I think almost only computation time spend on sovling matrix problem  ρ = M\b


Ns=[10, 20, 50, 100, 200, 500, 1000, 2000]
# Ns=[10]
E_error=zeros(size(Ns,1))

center=Int(round(NN/2)) # Note this only applies for a square space. Else take both dimensions into account.

for (jj,N) in enumerate(Ns)
    θ=range(0,2π,N+1)
    θ=θ[1:end-1]
    
    rxs = cos.(θ)
    rys = sin.(θ)
    rs = hcat(rxs,rys)
    
    # charge_pos_plot=scatter(rs[:,1],rs[:,2],aspect_ratio = :equal, markershape=:circle, label="Charge positions")
    # display(charge_pos_plot)
    
    
    b = -Φₑₓₜ(rxs,rys)
    
    M =  [i==j ? 0.0 : 𝒢_Φ(rxs[i], rxs[j], rys[i], rys[j]) for i in eachindex(rxs), j in eachindex(rys)] 
    # M =  [i==j ? 0.0 : 1*i+0*j for i in eachindex(rxs), j in eachindex(rys)] 
        # Note the indexes. Just have to loop over 
    # display(M)
    
    M = issymmetric(M) ? Symmetric(M) : M
    
    ρ = M\b
    # charge_pos_plot=scatter(rs[:,1],rs[:,2],aspect_ratio = :equal, markershape=:circle, label="Charges", zcolor=ρ)
    # display(charge_pos_plot)
    
    
    # Φs = Φₑₓₜ(center,center) # Note only evaluated at center
    # Φₜₒₜ = Φs
    # for i in 1:size(rxs,1)
    #     # print(i)
    #     x′ = rxs[i]
    #     y′ = rys[i]
    #     # Φₜₒₜ += ρ[i].*[𝒢_Φ(x, x′, y, y′) for x in xs, y in ys]
    #     Φₜₒₜ += ρ[i].*𝒢_Φ(center, x′, center, y′)# Note only evaluated at center
    # end
    
    # NOTE the circle marker color.
    # Φₜₒₜ_plot = heatmap(xs, ys, reverse(rotl90(Φₜₒₜ), dims = 1), aspect_ratio=:equal)
    # scatter!(Φₜₒₜ_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, label="", zcolor=ρ.*(N/15)) #, markercolor = :white
    # display(Φₜₒₜ_plot)
    
    ##
    
    # Exs_tot = [Exₑₓₜ(x, y) for x in xs, y in ys] 
    Exs_tot = Exₑₓₜ(xs[center], ys[center]) # Note only evaluated at center
    for i in 1:size(rxs,1)
        # print(i)
        x′ = rxs[i]
        y′ = rys[i]
        # Exs_tot += ρ[i].*[𝒢_Ex(x, x′, y, y′) for x in xs, y in ys]
        Exs_tot += ρ[i].*𝒢_Ex(xs[center], x′, ys[center], y′) # Note only evaluated at center
    end
    
    # Exs_tot_plot = heatmap(xs, ys, reverse(rotl90(Exs_tot), dims = 1), aspect_ratio=:equal)
    # scatter!(Exs_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
    # display(Exs_tot_plot)
    
    
    # Eys_tot = [Eyₑₓₜ(x, y) for x in xs, y in ys]
    Eys_tot = Eyₑₓₜ(xs[center], ys[center]) # Note only evaluated at center
    for i in 1:size(rxs,1)
        # print(i)
        x′ = rxs[i]
        y′ = rys[i]
        # Eys_tot += ρ[i].*[𝒢_Ey(x, x′, y, y′) for x in xs, y in ys]
        Eys_tot += ρ[i].*𝒢_Ey(xs[center], x′, ys[center], y′) # Note only evaluated at center
    end
    
    # Note the clim, since one point explode in value
    # Eys_tot_plot = heatmap(xs, ys, reverse(rotl90(Eys_tot), dims = 1), aspect_ratio=:equal) #clim=(-1.5,1.5)
    # scatter!(Eys_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
    # display(Eys_tot_plot)
    
    
    E_tot = @. sqrt(Exs_tot^2+Eys_tot^2)
    
    # E_tot_plot = heatmap(xs, ys, reverse(rotl90(E_tot), dims = 1), aspect_ratio=:equal,clim=(0,3)) #clim=(-1.5,1.5)
    # scatter!(E_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
    # display(E_tot_plot)

    # E_error[jj]=E_tot[center,center]
    E_error[jj]=E_tot #Note that analytical result is 0, whereas E_error = E_tot - 0 = E_tot
    # print(jj," : ", Exs_tot, "\n")
    # print(jj," : ", Eys_tot, "\n")
    # print(jj," : ", E_tot, "\n")
end


initial_error_plot=plot(Ns, E_error, label="")
xlabel!("Number of points [ ]")
ylabel!("Error [ ]")
plot!(title="Error plot")
display(initial_error_plot)

error_plot=plot(Ns,E_error, xaxis=:log, yaxis=:log, label="Center E-field error", markershape=:cross)
plot!(Ns, 10*Ns.^(-1), xaxis=:log, yaxis=:log, label="1st order convergence")
xlabel!("Number of points [ ]")
ylabel!("Error [ ]")
plot!(title="Convergence plot")
display(error_plot)











# -------------------
# -------------------
# -------------------

## Part e

# The below code is similar to above
# Used after defining the new metal surface and place the point charges. 
# See new shapes further below.

# NOTE. Not good practice to have function use global variables. 
#  But this is just to get more readable code.
function calculate_and_plot_all()
    charge_pos_plot=scatter(rs[:,1],rs[:,2],aspect_ratio = :equal, markershape=:circle, label="Charge positions")
    xlabel!("x [ ]")
    ylabel!("y [ ]")
    plot!(title="Charge locations, N=$N")
    display(charge_pos_plot)

    b = -Φₑₓₜ(rxs,rys)

    M =  [i==j ? 0.0 : 𝒢_Φ(rxs[i], rxs[j], rys[i], rys[j]) for i in eachindex(rxs), j in eachindex(rys)] 
    # M =  [i==j ? 0.0 : 1*i+0*j for i in eachindex(rxs), j in eachindex(rys)] 
        # Note the indexes. Just have to loop over 
    # display(M)

    M = issymmetric(M) ? Symmetric(M) : M

    ρ = M\b
    charge_pos_plot=scatter(rs[:,1],rs[:,2],aspect_ratio = :equal, markershape=:circle, label="Charges", zcolor=ρ)
    xlabel!("x [ ]")
    ylabel!("y [ ]")
    plot!(title="Charge locations and relative charge, N=$N")
    display(charge_pos_plot)


    Φs = [Φₑₓₜ(x,y) for x in xs, y in ys]
    Φs_plot = heatmap(xs, ys, reverse(rotl90(Φs), dims = 1))
    xlabel!("x [ ]")
    ylabel!("y [ ]")
    plot!(title="Φₑₓₜ, External electric potential, N=$N")
    display(Φs_plot)

    Φₜₒₜ = Φs
    for i in 1:size(rxs,1)
        # print(i)
        x′ = rxs[i]
        y′ = rys[i]
        Φₜₒₜ += ρ[i].*[𝒢_Φ(x, x′, y, y′) for x in xs, y in ys]
    end

    # NOTE the circle marker color.
    Φₜₒₜ_plot = heatmap(xs, ys, reverse(rotl90(Φₜₒₜ), dims = 1), aspect_ratio=:equal)
    scatter!(Φₜₒₜ_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, label="") #, markercolor = :white
    xlabel!("x [ ]")
    ylabel!("y [ ]")
    plot!(title="Φₜₒₜ, Total electric potential, N=$N")
    display(Φₜₒₜ_plot)


    ##

    Exs_tot = [Exₑₓₜ(x, y) for x in xs, y in ys]
    for i in 1:size(rxs,1)
        # print(i)
        x′ = rxs[i]
        y′ = rys[i]
        Exs_tot += ρ[i].*[𝒢_Ex(x, x′, y, y′) for x in xs, y in ys]
    end

    # Note the clim, to not have exploding values
    Exs_tot_plot = heatmap(xs, ys, reverse(rotl90(Exs_tot), dims = 1), aspect_ratio=:equal, clim=(-3,3)) # clim=(-3,3)
    scatter!(Exs_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
    xlabel!("x [ ]")
    ylabel!("y [ ]")
    plot!(title="E_x, Electric field in x-direction, N=$N")
    display(Exs_tot_plot)


    Eys_tot = [Eyₑₓₜ(x, y) for x in xs, y in ys]
    for i in 1:size(rxs,1)
        # print(i)
        x′ = rxs[i]
        y′ = rys[i]
        Eys_tot += ρ[i].*[𝒢_Ey(x, x′, y, y′) for x in xs, y in ys]
    end

    # Note the clim, to not have exploding values
    Eys_tot_plot = heatmap(xs, ys, reverse(rotl90(Eys_tot), dims = 1), aspect_ratio=:equal, clim=(-3,3)) #clim=(-1.5,1.5)
    scatter!(Eys_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
    xlabel!("x [ ]")
    ylabel!("y [ ]")
    plot!(title="E_y, Electric field in y-direction, N=$N")
    display(Eys_tot_plot)


    E_tot = @. sqrt(Exs_tot^2+Eys_tot^2)

    # Note the clim, to not have exploding values
    E_tot_plot = heatmap(xs, ys, reverse(rotl90(E_tot), dims = 1), aspect_ratio=:equal, clim=(0,3)) # clim=(0,3)
    scatter!(E_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
    xlabel!("x [ ]")
    ylabel!("y [ ]")
    plot!(title="abs(E_tot), Absoulte of total electric field, N=$N")
    display(E_tot_plot)
end



## Defining square shape metal

N=80
N=Int(round(N/4)*4)

N_fourth = Int(N/4)+1 # +1 here since 1 point is removed later, due to corners overlap

xlow=-1
xhigh=1
ylow=-1
yhigh=1

side_xs = range(xlow,xhigh,N_fourth)
side_ys = range(ylow,yhigh,N_fourth)
const_xlow_val = range(xlow,xlow,N_fourth)
const_xhigh_val = range(xhigh,xhigh,N_fourth)
const_ylow_val = range(ylow,ylow,N_fourth)
const_yhigh_val = range(yhigh,yhigh,N_fourth)

# side_xs = side_xs[1:end-1]
# side_ys = side_ys[1:end-1]
# const_xlow_val = const_xlow_val[1:end-1]
# const_xhigh_val = const_xhigh_val[1:end-1]
# const_ylow_val = const_ylow_val[1:end-1]
# const_yhigh_val = const_yhigh_val[1:end-1]

rxs = [side_xs[1:end-1]; const_xhigh_val[1:end-1]; reverse(side_xs)[1:end-1]; const_xlow_val[1:end-1]] # Note automatic deletation of multiple end points
rys = [const_yhigh_val[1:end-1]; reverse(side_ys)[1:end-1]; const_ylow_val[1:end-1]; side_ys[1:end-1]] # Note automatic deletation of multiple end points
rs = hcat(rxs,rys)

calculate_and_plot_all()

## 

# Make flowershape

N=80
N=Int(round(N/4)*4)

N_fourth = Int(N/4)

θ1=range(-π/2,π,N_fourth+1)
θ1=θ1[1:end-1] # Remove the double point in the start/end.

rxs1 = cos.(θ1)*0.5 .+ 0.5
rys1 = sin.(θ1)*0.5 .+ 0.5
rs1 = hcat(rxs1,rys1)


θ2=range(0,3/2*π,N_fourth+1)
θ2=θ2[1:end-1] # Remove the double point in the start/end.

rxs2 = cos.(θ2)*0.5 .- 0.5
rys2 = sin.(θ2)*0.5 .+ 0.5
rs2 = hcat(rxs2,rys2)


θ3=range(π/2,2*π,N_fourth+1)
θ3=θ3[1:end-1] # Remove the double point in the start/end.

rxs3 = cos.(θ3)*0.5 .- 0.5
rys3 = sin.(θ3)*0.5 .- 0.5
rs3 = hcat(rxs3,rys3)


θ4=range(π,2.5*π,N_fourth+1)
θ4=θ4[1:end-1] # Remove the double point in the start/end.

rxs4 = cos.(θ4)*0.5 .+ 0.5
rys4 = sin.(θ4)*0.5 .- 0.5
rs4 = hcat(rxs4,rys4)

rs = vcat(rs1, rs2, rs3, rs4)
rxs = rs[:,1]
rys = rs[:,2]

calculate_and_plot_all()


##

# Make pincushion distortion shape

N=80
N=Int(round(N/4)*4)

N_fourth = Int(N/4)

displacement_in_center=cos(π/4)*2

θ1=range(-π/4,-3/4*π,N_fourth+1)
θ1=θ1[1:end-1] # Remove the double point in the start/end.

rxs1 = cos.(θ1)
rys1 = sin.(θ1) .+ displacement_in_center
rs1 = hcat(rxs1,rys1)


θ2=range(π/4,-π/4,N_fourth+1)
θ2=θ2[1:end-1] # Remove the double point in the start/end.

rxs2 = cos.(θ2) .- displacement_in_center
rys2 = sin.(θ2)
rs2 = hcat(rxs2,rys2)


θ3=range(3/4*π,π/4,N_fourth+1)
θ3=θ3[1:end-1] # Remove the double point in the start/end.

rxs3 = cos.(θ3)
rys3 = sin.(θ3) .- displacement_in_center
rs3 = hcat(rxs3,rys3)


θ4=range(-3/4*π,-5/4*π,N_fourth+1)
θ4=θ4[1:end-1] # Remove the double point in the start/end.

rxs4 = cos.(θ4) .+ displacement_in_center
rys4 = sin.(θ4)
rs4 = hcat(rxs4,rys4)

rs = vcat(rs1, rs2, rs3, rs4)
rxs = rs[:,1]
rys = rs[:,2]

calculate_and_plot_all()


##

# Make sinusoid shape (open shape)

N=80

rxs = range(-1,1,N)

θ=range(0,2π,N)
rys = -0.5*sin.(θ)

rs = hcat(rxs,rys)

calculate_and_plot_all()



## Still part e

# Point charges in external field 
# Run this an afterwards some of the other sections above 
#  (Just to avoid more copy pasting and having more readable code)

function Φₑₓₜ(x,y)
    return @. -𝒢_Φ(x,-1.3,y,0) + 𝒢_Φ(x,1.3,y,0)
end

function Exₑₓₜ(x,y) # E = -∇Φ
    return @. -𝒢_Ex(x,-1.3,y,0) + 𝒢_Ex(x,1.3,y,0)
end

function Eyₑₓₜ(x,y) # E = -∇Φ
    return @. -𝒢_Ey(x,-1.3,y,0) + 𝒢_Ey(x,1.3,y,0)
end


##

# Make sharp edges

N=80
N=Int(round(N/2)*2)

N_half = Int(N/2)

rxs1 = range(-1,1,N_half+1)
rxs1 = rxs1[1:end-1]
rxs2 = range(1,-1,N_half+1)
rxs2 = rxs2[1:end-1]

θ=range(0,π,N_half+1)
θ=θ[1:end-1] 
rys1 = 0.25*sin.(θ)
rys2 = -0.25*sin.(θ)

rs1 = hcat(rxs1,rys1)
rs2 = hcat(rxs2,rys2)

rs = vcat(rs1,rs2)
rxs = rs[:,1]
rys = rs[:,2]

calculate_and_plot_all()
