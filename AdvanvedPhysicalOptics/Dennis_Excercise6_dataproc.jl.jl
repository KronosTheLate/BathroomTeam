if occursin("dennishb", homedir())
    using Pkg
    Pkg.activate("AdvancedPhysicalOptics", shared=true)
    cd("/home/dennishb/Uni/Semester/8. Sem/Advanced Physical Optics/Lab excercise 6")
end
using XLSX
using DataFrames
using GLMakie; Makie.inline!(true)
using Statistics  # For `mean`
using Unitful

#* Find which angle is in middle of first 2 minimi
#* offset all angles by that amount to get middle of 
#* 2 first minima at 0
##

sheet = XLSX.readxlsx("data.xlsx")

##¤ Linbo
begin
    sheet_linbo = sheet["Linbo"]
    angle_linbo, current_linbo = eachcol(sheet_linbo["A2:B84"].|>float)#|>vscodedisplay
    sortmask = sortperm(angle_linbo)
    original_angle_linbo = angle_linbo[sortmask]
    current_linbo = current_linbo[sortmask]
end
scatterlines(original_angle_linbo, current_linbo, axis=(title="Linbo", xlabel="Angle [°]", ylabel="Current [mA]"))
DataInspector(); display(current_figure())
original_first_minima_angles_linbo = (-8.0, 11.5)
corrected_angle_linbo = original_angle_linbo .- mean(original_first_minima_angles_linbo)
scatterlines(corrected_angle_linbo, current_linbo, axis=(title="Linbo", xlabel="Angle [°]", ylabel="Current [mA]"))

##* We have ordinary+ordinary -> extraordinary
#*  we know that linbo is negative uniaxial, meaning nₒ > nₑ
thickness_linbo = 500e-6
linbo_nₒ_532nm = 2.3232
linbo_nₒ_1064nm = 2.232
linbo_nₑ_532nm = 2.2342
linbo_nₑ_1064nm = 2.156

#¤ In words: ordinary 1064 nm light enters from air (n₁=1), refracts (n₂ = linbo_nₒ_1064nm), 
#¤ and then propegates with `n=linbo_nₑ_532nm`
#* OPL = thick/(cos(asin(n₁/n₂ * sin(θᵢ)))) * n_prop
OPL_linbo = @. thickness_linbo/(cos(asin(1/linbo_nₒ_1064nm * sind(corrected_angle_linbo)))) * linbo_nₑ_532nm
scatterlines(OPL_linbo.*10^3, current_linbo, axis=(title="Linbo", xlabel="Optical path length [mm]", ylabel="Current [mA]"))
minima_linbo = [1.12037, 1.1278, 1.1353, 1.1428]
vlines!(minima_linbo, color=:red, label="Minima")
maxima_linbo = minima_linbo .+ 0.47*diff(minima_linbo[1:2])
vlines!(maxima_linbo, color=:green, label="Maxima")
axislegend("Manually determined", position=:rt)
minima_linbo|>diff
display(current_figure())

##¤ KTP
begin
    sheet_KTP = sheet["KTP"]
    angle_KTP, current_KTP = eachcol(sheet_KTP["A2:B84"].|>float)#|>vscodedisplay
    sortmask = sortperm(angle_KTP)
    original_angle_KTP = angle_KTP[sortmask]
    current_KTP = current_KTP[sortmask]
end
scatterlines(original_angle_KTP, current_KTP, axis=(title="KTP", xlabel="Angle [°]", ylabel="Current [mA]"))
DataInspector(); display(current_figure())
original_first_minima_angles_KTP = (-13, 15)
corrected_angle_KTP = original_angle_KTP .- mean(original_first_minima_angles_KTP)
scatterlines(corrected_angle_KTP, current_KTP, axis=(title="KTP", xlabel="Angle [°]", ylabel="Current [mA]"))

#! OPL = thickness/(cos(asin(n₁/nₜ * sin(θ_incident)))) * n₂