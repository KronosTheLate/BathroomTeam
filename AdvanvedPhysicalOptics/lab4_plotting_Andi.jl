using LinearAlgebra 
using Plots
gr()
using QuadGK

## Task 1b

ω₁=96
ω₂=100

Δω=1/2*(ω₂-ω₁)
ω_avg=1/2*(ω₂+ω₁)

I(t)=1/4*(1+cos(2*ω_avg*t))*(1+cos(2*Δω*t))

function integral_of_beat(N,startp,endp)
    # startp=0
    # endp=4
    # N=60
    t=range(startp,endp,N)
    A=Vector{Float64}(undef,N-1)
    for (i,ti) in enumerate(t[1:end-1])
        integral=quadgk(I, ti, ti+step(t), rtol=1e-5)
        A[i]=integral[1]
    end
    tnew=collect(t)
    tnew=tnew.+step(t)/2
    tnew=tnew[1:end-1]
    return tnew, A
end

startp=0
endp=1

plot()

plot!(I, linewidth=3, label="Real intensity")
# plot!(title = "Intensity")
# plot!(legend=:none)
# xlabel!("t [ ]")
# # ylabel!("Intensity [ ]")
# Int_plot=ylabel!("Intensity [ ]")
# display(Int_plot)


for (i,N) in enumerate(reverse([10, 100, 500]))
    tnew, A =integral_of_beat(N,startp,endp)
    # t_int=round(tnew[2]-tnew[1],digits=4)
    f_int=round(1/(tnew[2]-tnew[1]),digits=4)
    plot!(tnew, A./maximum(A), linewidth=4, label="1/T_avg = $f_int",
    linestyle=:solid, linealpha=1, xlims=(startp,endp), dpi=300)
    T_fast=round(1/ω_avg,digits=4)
    plot!(title = "Intensity, ω_avg=$ω_avg, Δω=$Δω")
    # plot!(legend=)
    xlabel!("Normalized t [ ]")
    Int_plot=ylabel!("Normalized intensity [ ]")
    display(Int_plot)
end
savefig("Int_pattern_meas")


endp=25
plot()
I2(t)=1/2*(1+cos(2*Δω*t))
plot!(I2, linewidth=3, label="Beat envelope",xlims=(startp,endp))

for (i,N) in enumerate(reverse([7, 50, 300]))
    tnew, A =integral_of_beat(N,startp,endp)
    f_int=round(1/(tnew[2]-tnew[1]),digits=4)
    plot!(tnew, A./maximum(A), linewidth=3, label="1/T_avg = $f_int",
    linestyle=:solid, linealpha=1, xlims=(startp,endp), dpi=300)
    plot!(title = "Intensity, Δω=$Δω")
    plot!(legend=:bottomright)
    xlabel!("Normalized t [ ]")
    Int_plot=ylabel!("Normalized intensity [ ]")
    display(Int_plot)
end
savefig("BeatEnvelope_pattern_meas")



## Task 1f


ω_avg=1

Δω_array=ω_avg*[0; 0.01; 0.05; 0.1]

# plot(layout = (4, 1))
# plot()
# Int_plot=Array{Any}(nothing, 4)
# Int_plot=repeat([plot(1)], 4)


for (i,Δω) in enumerate(Δω_array)
    ω₁=ω_avg+Δω
    ω₂=ω_avg-Δω

    I(h)=2*(1+cos(ω_avg*2*h/1)*cos(Δω*2*h/1))

    plot(I, linewidth=4, label="Intensity", xlim=(0, 100))
    plot!(title = "Intensity")
    plot!(legend=:none)

    plot!(title = "Intensity, ω_avg=$ω_avg, Δω=$Δω")
    # plot!(legend=)
    xlabel!("Normalized h [ ]")
    Int_plot=ylabel!("Normalized intensity [ ]")
    display(Int_plot)
    savefig(string("Beat_pattern_Delta_omega",string(i)))
end

# Tot_Int_plot=plot(Int_plot[1],Int_plot[2],Int_plot[3],Int_plot[4], layout = (2,2))
# print(Int_plot)



## Task 1f. Virker ikke som det skal... Jeg kan ikke få lavet subplots
if false
ω_avg=100

Δω_array=ω_avg*[0; 0.01; 0.05; 0.1]

# plot(layout = (4, 1))
# Int_plot=plot()
# Int_plot=Array{Any}(nothing, 4)
# Int_plot=repeat([plot(1)], 4)

I_tot=Array{Any}(nothing, 4)
for (i,Δω) in enumerate(Δω_array)
    ω₁=ω_avg+Δω
    ω₂=ω_avg-Δω

    I(τ)=2*(1+cos(ω_avg*τ)*cos(Δω*τ))
    I_tot[i]=I


end




Int_plot1=plot(I_tot[1],linewidth=2, label="Intensity", xlim=(0, 2))
Int_plot2=plot(I_tot[2],linewidth=2, label="Intensity", xlim=(0, 2))
Int_plot3=plot(I_tot[3],linewidth=2, label="Intensity", xlim=(0, 2))
Int_plot4=plot(I_tot[4],linewidth=2, label="Intensity", xlim=(0, 2))
Int_tot=plot(Int_plot1,Int_plot2,Int_plot3,Int_plot4, layout = (2,2))
# plot!(I, linewidth=4, label="Intensity", xlim=(0, 2))
plot!(title = "Intensity")
plot!(legend=:none)
plot!(title = "Intensity, ω_avg=$ω_avg, Δω=$Δω")
# plot!(legend=)
xlabel!("Normalized t [ ]")
ylabel!("Normalized intensity [ ]")
display(Int_tot)
end

# Tot_Int_plot=plot(Int_plot[1],Int_plot[2],Int_plot[3],Int_plot[4], layout = (2,2))
# print(Int_plot)


## Task 2b

# NOTE \nu / ν and not v [vee] is used!!!
Δν=0.2
ν₀=10
ν_next=1
I₀=1
τ=0:0.001:15
ν=5:0.002:35

Iarr = @. I₀*(1+exp(-(Δν)^2*τ^2)*cos(2*π*ν₀*τ))
spectral_distribution=@. I₀*exp(-( (ν-ν₀) / (2*Δν) ))
Int_plot=plot()
spectral_plot=plot()
# spectral_distribution=zeros(size(spectral_distribution,1))
Iarr=zeros(size(Iarr,1))

for (i,N) in enumerate([5])#, 2, 5, 20])
spectral_distribution=zeros(size(spectral_distribution,1))
    for j in 1:N
        Iarr+=@. I₀*(1+exp(-(Δν)^2 *τ^2)*cos(2*π*(ν₀+(j-1)*ν_next)*τ))
        spectral_distribution+=@. I₀*exp(-( (ν-(ν₀+(j-1)*ν_next)) / (2*Δν) )^2)
        
    end
linewidths4plot=[8, 4, 2, 1]
plot!(spectral_plot, ν, spectral_distribution, linewidth=linewidths4plot[i],
label="Number of gauss peaks: N=$N", dpi=300)
Iarr=@. Iarr/Iarr[1]
# I(τ)=I₀*(3+exp(-(Δν)^2 *τ^2)*cos(2*π*ν₀*τ)+exp(-(Δν)^2 *τ^2)*cos(2*π*(ν₀+0.5)*τ)+exp(-(Δν)^2 *τ^2)*cos(2*π*(ν₀+1)*τ)+exp(-(Δν)^2 *τ^2)*cos(2*π*(ν₀+1.5)*τ)+exp(-(Δν)^2 *τ^2)*cos(2*π*(ν₀+2)*τ)+exp(-(Δν)^2 *τ^2)*cos(2*π*(ν₀+2.5)*τ)+exp(-(Δν)^2 *τ^2)*cos(2*π*(ν₀+3)*τ))

# N=200
startp=0
endp=8
# tau=range(startp,endp,N)


plot!(Int_plot,τ, Iarr, linewidth=3, label="Number of gauss peaks: N=$N",
linestyle=:solid, linealpha=1, xlims=(startp,endp), dpi=300)


end
plot!(Int_plot, title = "Intensity, ν₀=$ν₀, Δν=$Δν, next peak +$ν_next")
plot!(Int_plot, legend=:bottomright)
xlabel!(Int_plot, "Normalized τ [ ]")
ylabel!(Int_plot, "Normalized intensity [ ]")
display(Int_plot)
savefig("Intensity_for_Gaussian_density")

plot!(spectral_plot, title = "Spectral density, ν₀=$ν₀, Δν=$Δν, next peak +$ν_next")
plot!(spectral_plot, legend=:bottomright)
xlabel!(spectral_plot, "Normalized ν [ ]")
ylabel!(spectral_plot, "Normalized intensity density [ ]")
display(spectral_plot)
savefig("Intensity_spectral_density")



#-----------------------------------------
## Test to see if one can find coherence time from cosine transform.

a=3*10^8 / (1000*10^-9)
b=3*10^8 / (600*10^-9)
c=3*10^8 / (200*10^-9)

(a+c)/2


frequ(λv)=3*10^8 / (λv)

plot(frequ,500*10^-9,800*10^-9)


#-----------------------------------------
## Coherence length calculated from cosine transfrom 
# for data obtained in experiment
spacing_in_spectrometer_data_in_nm=0.35
transfrom_wavelength=633
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=1/bdiff # delay spacing after cosine transform

index_in_cos_trans=220 # 10 as an example

lc=3*10^8*τ1*(index_in_cos_trans)



#-----------------------------------------
## HeNe laser

spacing_in_spectrometer_data_in_nm=0.33
transform_wavelength=633
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=1/bdiff # delay spacing after cosine transform

index_in_cos_trans=136 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)


##For FWHM
spacing_in_spectrometer_data_in_nm=5
transform_wavelength=633-spacing_in_spectrometer_data_in_nm/2
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=3.33/bdiff # delay spacing after cosine transform

index_in_cos_trans=1 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)




#-----------------------------------------
## LASER_diode_larger

spacing_in_spectrometer_data_in_nm=0.33
transform_wavelength=655
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=1/bdiff # delay spacing after cosine transform

index_in_cos_trans=132 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)


##For FWHM
spacing_in_spectrometer_data_in_nm=7.5
transform_wavelength=655-spacing_in_spectrometer_data_in_nm/2
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=3.33/bdiff # delay spacing after cosine transform

index_in_cos_trans=1 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)


#-----------------------------------------
## LASER_diode_smaller

spacing_in_spectrometer_data_in_nm=0.33
transform_wavelength=673
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=1/bdiff # delay spacing after cosine transform

index_in_cos_trans=133 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)


##For FWHM
spacing_in_spectrometer_data_in_nm=8
transform_wavelength=673-spacing_in_spectrometer_data_in_nm/2
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=3.33/bdiff # delay spacing after cosine transform

index_in_cos_trans=1 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)



#-----------------------------------------
## LED_red

spacing_in_spectrometer_data_in_nm=0.33
transform_wavelength=638
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=1/bdiff # delay spacing after cosine transform

index_in_cos_trans=46 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)


##For FWHM
spacing_in_spectrometer_data_in_nm=17.5
transform_wavelength=638-spacing_in_spectrometer_data_in_nm/2
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=3.33/bdiff # delay spacing after cosine transform

index_in_cos_trans=1 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)


#-----------------------------------------
## LED_yellow

spacing_in_spectrometer_data_in_nm=0.33
transform_wavelength=585
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=1/bdiff # delay spacing after cosine transform

index_in_cos_trans=25 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)


##For FWHM
spacing_in_spectrometer_data_in_nm=33
transform_wavelength=585-spacing_in_spectrometer_data_in_nm/2
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=3.33/bdiff # delay spacing after cosine transform

index_in_cos_trans=1 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)




#-----------------------------------------
## Coherence length calculated from cosine transfrom 
# for data delivered; "Laserdiode HR"

spacing_in_spectrometer_data_in_nm=0.02
transform_wavelength=670
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
τ1=1/bdiff # delay spacing after cosine transform

#Inner envelope
index_in_cos_trans=14 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)

## Outer envelope

index_in_cos_trans=188 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)

##For FWHM
spacing_in_spectrometer_data_in_nm=0.04
transform_wavelength=670
b1=3*10^8 / (transform_wavelength*10^-9)
b2=3*10^8 / ((transform_wavelength+spacing_in_spectrometer_data_in_nm)*10^-9)

bdiff=b1-b2
bdiff/3.33
τ1=3.33/bdiff # delay spacing after cosine transform

index_in_cos_trans=1 # 10 as an example
lc=3*10^8*τ1*(index_in_cos_trans)




