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

##! OPL = thickness/(cos(asin(n₁/nₜ * sin(θ_incident)))) * n₂

##¤ From Vladimir's script, for LN
no2 = 2.3232; ne2 = 2.2342; # at 532 nm
no1 = 2.232; ne1 = 2.156; # at 1064 nm
L = 500e-6; # crystal length in um

n_eff(nₑ_532nm, nₒ_532nm, θ) = 1 ./ sqrt(cos(θ).^2/nₑ_532nm ^2+sin(θ).^2/nₒ_532nm^2);#refractive index of extraordinary way

λ_i = 1.064e-6; #fundamental harmonic in um
# L_c_linbo(θ, λ_i=1064e-9) = λ_i ./ (4*(linbo_nₒ_1064nm - n_eff(linbo_nₑ_532nm+0.012, linbo_nₒ_532nm, θ))); # coherence length
L_c_linbo(θ, λ_i=1064e-9) = λ_i ./ (4*(linbo_nₒ_1064nm - n_eff(linbo_nₑ_532nm, linbo_nₒ_532nm, θ))); # coherence length
L_c_linbo(0) * 10^6
# ϕ_linbo(θ) = 2π/λ_i*(thickness_linbo ./cos(θ)) * 2 .* (linbo_nₒ_1064nm - n_eff(linbo_nₑ_532nm+0.012,linbo_nₒ_532nm,θ))/2; #phase mismatch
ϕ_linbo(θ) = 2π/λ_i*(thickness_linbo ./cos(θ)) * 2 .* (linbo_nₒ_1064nm - n_eff(linbo_nₑ_532nm,linbo_nₒ_532nm,θ))/2; #phase mismatch

I2w_linbo(θ) = ((thickness_linbo./cos(θ)).*sin(ϕ_linbo(θ))./(ϕ_linbo(θ))).^2; #power of second harmonic

θ_i = range(0, π/2, 500)     #th_deg = 0:0.1:90; # incident angle
θ_f = @. asin(sin(θ_i)/linbo_nₒ_1064nm) # th = angle in crystal - by Snells law
actual_length = thickness_linbo ./ cos.(θ_f)
using LinearAlgebra
let 
    fig = Figure()
    ax = Axis(fig[1, 1], xtickformat = vals -> string.(round.(vals.*10^6, digits=3), 'μ', 'm'), 
    xlabel="Geometrical path length inside crystal based on incident angle", ylabel="Normalized intensity")
    lines!(actual_length .- 1.46e-6, normalize(I2w_linbo.(θ_f), Inf), label="Analytical (offset)")
    GPL_linbo = OPL_linbo ./ linbo_nₑ_532nm
    minima_GPL_linbo = minima_OPL_linbo ./ linbo_nₑ_532nm
    maxima_GPL_linbo = maxima_OPL_linbo ./ linbo_nₑ_532nm
    lines!(GPL_linbo, normalize(current_linbo, Inf), label="Experimental")
    minima_actual_length = [504.7, 507.89, 511.05, 514.25] ./ 10^6  .- 3.3e-6
    vlines!(minima_actual_length, color=Cycled(1), linestyle=:dash, label="Minima analytical")
    vlines!(minima_GPL_linbo, color=Cycled(2), linestyle=:dash, label="Minima experimental")
    # vlines!(maxima_GPL_linbo, color=:green, label="Maxima")
    # axislegend("Manually determined", position=:rt)
    xlims!(498e-6, 525e-6)
    axislegend()
    save(joinpath(homedir(), "Pictures", "Lab6_geometrical_path_length_w_analytical_linbo.png"), fig)
    display(fig)
    μ1 = mean(diff(minima_actual_length)) / 2
    μ2 = mean(diff(minima_GPL_linbo)) / 2

    @show μ1 * 10^6
    @show μ2 * 10^6
    @show (μ2 - μ1) / μ1 * 100
end
##
let 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Angle inside crystal [rad]", ylabel="Refractive index")
    scatterlines!(θ_f, n_eff.(linbo_nₑ_532nm, linbo_nₒ_532nm, θ_f), label="n_eff at 532 nm")
    scatterlines!(θ_f, ϕ_linbo.(θ_f) .|> deg2rad, label="Phase mismatch at 532 nm")
    scatterlines!(θ_f, θ_i, label="Indicent angle")
    axislegend()
    display(fig)
end

##
1
#=
figure; %birefrigence
yyaxis left;
plot(th_deg,ne_th(ne2,no2,th_deg*pi/180),'b','DisplayName','n_e(\theta) @532 nm');
title('Effect of birefrigence');
hold on;
yline(ne2, 'b--','DisplayName','n_e(0) @532 nm');
yline(no2, 'b:','DisplayName','n_o @532 nm');
xlabel('Angle inside crystal [deg]');
ylabel('Refractive index');
legend('Location','best');
yyaxis right;
plot(asin(sin(th_deg*pi/180)/no1)*180/pi,th_deg,'DisplayName','\theta_i');
ylabel('Incident angle, \theta_i [deg]');

figure('Position',[700,200,500,700]);
tiledlayout('flow','TileSpacing','none','Padding','compact');
nexttile;
plot(th_deg,ne_th(ne2,no2,th),'b','DisplayName','n_e(\theta) @ 532 nm');
title('Maker firnges for lithium niobate');
hold on;
yline(ne2, 'b--','DisplayName','n_e(0) @532 nm');
yline(no2, 'b:','DisplayName','n_o @532 nm');
yline(ne1, 'r--','DisplayName','n_e @1064 nm');
yline(no1, 'r:','DisplayName','n_o @1064 nm');
% xlabel('Incident angle [deg]');
ylabel('Refractive index');
legend('Location','best');
set(gca,'XTick',[]);

nexttile;
yyaxis left;
plot(th_deg,abs(lc(th)));
% xlabel('Incident angle [deg]');
ylabel('Coherence length, \pi/\Deltak');
set(gca,'XTick',[]);
yyaxis right;
plot(th_deg,L./cos(th));
ylabel('Optical path, L/cos(\theta)');

nexttile;
yyaxis left;
plot(th_deg,I2w(th));
ylabel('SHG intensity');
xlabel('Incident angle, \theta_i [deg]');

yyaxis right;
plot(th_deg,(abs(phi(th)))/pi);
ylabel('Phase mismatch in \pi, \DeltakL/2\pi');
hold on;
yline(2, ':');
yline(3, ':');
yline(4, ':');
yline(5, ':');
=#
