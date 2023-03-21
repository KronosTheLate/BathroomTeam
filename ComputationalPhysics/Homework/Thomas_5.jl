#Andis includes
using LinearAlgebra
using Random
using Distributions
using Plots

#PROBLEM SET 4

include("Thomas_5_leapfrog.jl")
include("Thomas_5_kepler.jl")

#1A 

#No of particles
N = 50 

#Initial conditions
init_cond = rand(-1:0.1:1,6,N)

#Center of box
R = [0; 0; 0]

#Leapfrog integration

numeric = leapfrog_integration(kepler_eom,R,init_cond,0,20,0.1)