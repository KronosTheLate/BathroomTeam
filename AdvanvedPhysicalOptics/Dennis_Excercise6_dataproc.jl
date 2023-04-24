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

sheet = XLSX.readxlsx("data.xlsx")

##¤ Linbo
begin
    sheet_linbo = sheet["Linbo"]
    angle_linbo, current_linbo = eachcol(sheet_linbo["A2:B84"].|>float)#|>vscodedisplay
    sortmask = sortperm(angle_linbo)
    original_angle_linbo = angle_linbo[sortmask]
    current_linbo = current_linbo[sortmask]
end


begin
    scatterlines(original_angle_linbo, current_linbo, axis=(title=L"\text{LiNbO}_3", xlabel="Angle [°]", ylabel="Current [mA]"))
    display(current_figure())
    save(joinpath(homedir(), "Pictures", "Lab6_rawdata_linbo.png"), current_figure())
end
# DataInspector(); display(current_figure())
original_first_minima_angles_linbo = (-8.0, 11.5)
corrected_angle_linbo = original_angle_linbo .- mean(original_first_minima_angles_linbo)
scatterlines(corrected_angle_linbo, current_linbo, axis=(title=L"\text{LiNbO}_3", xlabel="Corrected angle [°]", ylabel="Current [mA]"))

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

begin
    OPL_linbo = @. thickness_linbo/(cos(asin(1/linbo_nₒ_1064nm * sind(corrected_angle_linbo)))) * linbo_nₑ_532nm
    minima_OPL_linbo = [1.12037, 1.1278, 1.1353, 1.1428] ./ 10^3
    maxima_OPL_linbo = minima_OPL_linbo .+ 0.5*diff(minima_OPL_linbo[1:2])
    scatterlines(OPL_linbo .* 10^3, current_linbo, axis=(title=L"\text{LiNbO}_3", xlabel="Optical path length inside crystal [mm]", ylabel="Current [mA]", xtickformat=vals->string.(round.(vals*1000, digits=2), "μm")))
    vlines!(minima_OPL_linbo .* 10^3, color=:red, label="Minima")
    vlines!(maxima_OPL_linbo .* 10^3, color=:green, label="Maxima")
    axislegend("Manually determined", position=:rt)
    display(current_figure())
    save(joinpath(homedir(), "Pictures", "Lab6_optical_path_length_linbo.png"), current_figure())
end
begin
    GPL_linbo = OPL_linbo ./ linbo_nₑ_532nm
    minima_GPL_linbo = minima_OPL_linbo ./ linbo_nₑ_532nm
    maxima_GPL_linbo = maxima_OPL_linbo ./ linbo_nₑ_532nm
    scatterlines(GPL_linbo .* 10^3, current_linbo, axis=(title=L"\text{LiNbO}_3", xlabel="Geometrical path length inside crystal", ylabel="Current [mA]", xtickformat=vals->string.(round.(vals*1000, digits=2), "μm")))
    vlines!(minima_GPL_linbo .* 10^3, color=:red, label="Minima")
    vlines!(maxima_GPL_linbo .* 10^3, color=:green, label="Maxima")
    axislegend("Manually determined", position=:rt)
    display(current_figure())
    save(joinpath(homedir(), "Pictures", "Lab6_geometrical_path_length_linbo.png"), current_figure())
end
coherence_length_linbo = mean(diff(sort([minima_GPL_linbo; maxima_GPL_linbo])))
coherence_length_linbo * 1e6  # in μm
theoretical_coherence_length_linbo = 121e-6
rel_dev_linbo = (coherence_length_linbo - theoretical_coherence_length_linbo)/theoretical_coherence_length_linbo


##¤ KTP
##TODO - take n inside KTP to be (n_y_1064nm+n_z_1064nm)/2
begin
    sheet_KTP = sheet["KTP"]
    angle_KTP, current_KTP = eachcol(sheet_KTP["A2:B84"].|>float)#|>vscodedisplay
    sortmask = sortperm(angle_KTP)
    original_angle_KTP = angle_KTP[sortmask]
    current_KTP = current_KTP[sortmask]
end

begin
    scatterlines(original_angle_KTP, current_KTP, axis=(title="KTP", xlabel="Angle [°]", ylabel="Current [mA]"))
    display(current_figure())
    save(joinpath(homedir(), "Pictures", "Lab6_rawdata_KTP.png"), current_figure())
end

thickness_KTP = 1000e-6
KTP_nx_532nm = 1.77803
KTP_ny_532nm = 1.7886
KTP_nz_532nm = 1.8887
KTP_nx_1064nm = 1.73773
KTP_ny_1064nm = 1.7453
KTP_nz_1064nm = 1.8297

begin
    original_first_minima_angles_KTP = (-13.3, 15) # Manually tuned to match end
    corrected_angle_KTP = original_angle_KTP .- mean(original_first_minima_angles_KTP)
    GPL_KTP = thickness_KTP ./ cos.(asin.(1/mean([KTP_ny_1064nm, KTP_nz_1064nm]) * sind.(corrected_angle_KTP)))
    minima_GPL_KTP = [1.009, 1.023, 1.037, 1.051] .* 1e-3#[1.12037, 1.1278, 1.1353, 1.1428]
    maxima_GPL_KTP = minima_GPL_KTP .+ 0.5*diff(minima_GPL_KTP[1:2])
    scatterlines(GPL_KTP, current_KTP, axis=(title="KTP", xlabel="Geometrical path length inside crystal", ylabel="Current [mA]", xtickformat=vals->string.(round.(vals*1000, sigdigits=3), "mm")))
    xlims!(nothing, 1.085e-3)
    vlines!(minima_GPL_KTP, color=:red, label="Minima")
    vlines!(maxima_GPL_KTP, color=:green, label="Maxima")
    axislegend("Manually determined", position=:rt)
    display(current_figure())
    save(joinpath(homedir(), "Pictures", "Lab6_geometrical_path_length_KTP.png"), current_figure())
end

coherence_length_KTP = mean(diff(sort([minima_GPL_KTP; maxima_GPL_KTP])))
coherence_length_KTP * 1e6  # in μm
theoretical_coherence_length_KTP = 242e-6
rel_dev_KTP = (coherence_length_KTP - theoretical_coherence_length_KTP)/theoretical_coherence_length_KTP

# DataInspector(); display(current_figure())
scatterlines(corrected_angle_KTP, current_KTP, axis=(title="KTP", xlabel="Angle [°]", ylabel="Current [mA]"))
map(vals->string.(vals*10^3, "m"), 1:10)
#! OPL = thickness/(cos(asin(n₁/nₜ * sin(θ_incident)))) * n₂