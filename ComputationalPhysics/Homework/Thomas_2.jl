using GLMakie; Makie.inline!(true) #This is for plotting. 'true' argument denotes inline plots in vscode
using LinearAlgebra #For Linear Algebra math

##---------------A--------------##

#Construct a function handle that returns P matrix

#n: refrace index
#D: thickness of layer
#k: vacuum wave number

function P_matrix(n,D,k)
    P = [exp(im*n*D*k) 0; 0 exp(-im*n*D*k)]
    return P
end

#P is a vector that translates the amplitudes from one side of the medium to another

##--------------B---------------##

#nᵢ: refractive index to the *right* of the interface
#nᵢ₋₁: refractive index to the *left* of the interface


function I_matrix(nᵢ,nᵢ₋₁)
    I = [   1+nᵢ/nᵢ₋₁ 1-nᵢ/nᵢ₋₁;
            1-nᵢ/nᵢ₋₁ 1+nᵢ/nᵢ₋₁]
    return I
end

#I is a matrix that relates amplitudes to the left and right of a material interface
#due to reflection and transmissions

##-------------C----------------##

#Create a total transfer matrix T

# ~~--!!n_vect = D_vect + 2 must hold!!--~~

function T(n_vect,D_vect,k)
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

n_vect = [1 1 1 1]
D_vect = [1 1]
k = 1

T(n_vect,D_vect,k)
