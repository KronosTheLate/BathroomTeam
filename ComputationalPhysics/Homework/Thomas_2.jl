#using GLMakie; Makie.inline!(true) #This is for plotting. 'true' argument denotes inline plots in vscode
using LinearAlgebra #For Linear Algebra math
using Plots
plotly() #For interactive plots

##---------------A--------------##

#Construct a function handle that returns P matrix

#n: refrace i
#D: thickness of layer
#k: vacuum wave number

function P_matrix(n,D,k)
    P = [exp(im*n*D*k) 0; 0 exp(-im*n*D*k)]
    return P
end

#P is a vector that translates the amplitudes from one side of the medium to another

##--------------B---------------##

#nᵢ: refractive i to the *right* of the interface
#nᵢ₋₁: refractive i to the *left* of the interface


function I_matrix(nᵢ,nᵢ₋₁)
    I = 0.5 * [     1+nᵢ/nᵢ₋₁ 1-nᵢ/nᵢ₋₁;
                    1-nᵢ/nᵢ₋₁ 1+nᵢ/nᵢ₋₁]
    return I
end

#I is a matrix that relates amplitudes to the left and right of a material interface
#due to reflection and transmissions

##-------------C----------------##

#Create a total transfer matrix T

# ~~--!!n_vect = D_vect + 2 must hold!!--~~

function T_matrix(n_vect,D_vect,k)
    T = LinearAlgebra.I
    for i in 1:length(n_vect)-1
        
        I=I_matrix(n_vect[i+1],n_vect[i])
        
        if i <= length(D_vect)
            P=P_matrix(n_vect[i+1],D_vect[i],k)
        else
            P=LinearAlgebra.I
        end

        T=P*I*T #Add another propagation set 
    end
    return T
end

#Test T function

n_vect = [1 50 20 4 1]
D_vect = [1 1 1]
k = 1
t_mat = T_matrix(n_vect,D_vect,k)
println("test")
@show t_mat

##-----------------D----------------##



#Reflectance function

function Reflectance(Mat)
    R = abs2(Mat[2,1]/Mat[2,2])
    return R
end

#Transmittance function

function Transmittance(Mat)
    T = abs2(Mat[1,1] - (Mat[1,2]*Mat[2,1])/Mat[2,2])
    return T
end


#Define plotting function
function superPlot(n_vect,D_vect,N,bool)

    reflectance_vect =      Array{Float64}(undef, N)
    transmittance_vect =    Array{Float64}(undef, N)
    loss_vect =             Array{Float64}(undef, N)
    k_range =               LinRange(0,3,N)

    for i in 1:N
        transfer_matrix = T_matrix(n_vect,D_vect,k_range[i])
        reflectance_vect[i] = Reflectance(transfer_matrix)
        transmittance_vect[i] = Transmittance(transfer_matrix)
        loss_vect[i] = 1 - transmittance_vect[i] - reflectance_vect[i]
    end

    plot1 = Plots.plot(k_range,transmittance_vect, labels = "Transmittance")
    Plots.plot!(k_range,reflectance_vect,labels = "Reflectance")
    Plots.plot!(title = "R & T")
    Plots.plot!(legend=:right)
    if bool == true
        Plots.plot!(yaxis=:log)
    end
    display(plot1)

    plot2 = Plots.plot(k_range,loss_vect, labels = "Error (lost energy)")
    Plots.plot!(title = "Error for each k, ideally zero")
    Plots.plot!(legend=:right)
    if bool == true
        Plots.plot!(yaxis=:log)
    end
    display(plot2)
end



#5 pairs
n2 = 2
n1 = 1

n_vect = [1 n2 n1 n2 n1 n2 n1 n2 n1 n2 1]
D_vect = [1 1 1 1 1 1 1 1 1]

#superPlot(n_vect,D_vect,1000,false)

#19 pairs

n_vect = [1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 1 ]
D_vect = ones(length(n_vect)-2)

#superPlot(n_vect,D_vect,1000,false)

#39 pairs

n_vect = [1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 1 ]
D_vect = ones(length(n_vect)-2)

#superPlot(n_vect,D_vect,1000,false)

#----------E--------------#

#trying n = 2 + 5i
n_vect = [1 n2 n1 n2 n1 (2+5*im) n1 n2 n1 n2 1]
D_vect = [1 1 1 1 1 1 1 1 1]

superPlot(n_vect,D_vect,1000,false)

#trying n = 2 + 10i
n_vect = [1 n2 n1 n2 n1 (2+10*im) n1 n2 n1 n2 1]
D_vect = [1 1 1 1 1 1 1 1 1]

superPlot(n_vect,D_vect,1000,false)

#trying n = 2 + 15i
n_vect = [1 n2 n1 n2 n1 (2+15*im) n1 n2 n1 n2 1]
D_vect = [1 1 1 1 1 1 1 1 1]

superPlot(n_vect,D_vect,1000,true)


