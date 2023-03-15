using XLSX, DataFrames, CSV
using GLMakie, FFTW; Makie.inline!(true)
pic_dir = joinpath(homedir(), "Pictures")  # Override in below if-statement if you want pictures elsewhere
if occursin("dennishb", homedir())  # My hacky way of checking which computer I am on
    my_dir = "/home/dennishb/Uni/Semester/8. Sem/Advanced Physical Optics/Lab excercise 4"
elseif occursin("your_username", homedir())
    # enter path to files ("Sweep_mittor_data.xlsx" and "Laserdiode_HR_spectrum.csv") on your computer
else
    throw("Don't know what computer I am on, and no file path is defined. Files will not be found, and the script will error. Please provide path.")
end
sweep_mirror_datapath = joinpath(my_dir, "Sweep_mirror_data.xlsx")
laserdiode_spectrum_datapath = joinpath(my_dir, "Laserdiode_HR_spectrum.csv")

##? Defining data and titles
all_dataframes = [XLSX.readtable(sweep_mirror_datapath, "Sheet$i", "A:B", first_row=2)|>DataFrame.|>identity for i in 1:3]
titles = [XLSX.readxlsx(sweep_mirror_datapath)[i]["A1"] for i in 1:3]

#? Styling plot
#? In current workflow, I open a standalone window, and reuse that (empty!), never closing it.
# fig = current_figure()
update_theme!(fontsize=20, markersize=10)

##? Plotting
begin
    fig = Figure()
    for i in eachindex(all_dataframes)
        dataframe = all_dataframes[i]
        title = titles[i]
        ax, _ = scatterlines(fig[i, 1], eachcol(dataframe)..., linestyle="-")
        ax.xlabel, ax.ylabel = propertynames(dataframe) .|> string
        ax.title = title
        # Makie.inline!(false)
        # DataInspector()
        fig|>display
        # Makie.inline!(true)
    end
end
dataframe = all_dataframes[1][[begin+3, end-2], :]
d1, d2 = 20.19, 21.35
(d2-d1)/2
d3, d4 = 18.15, (19.2+19.25)/2
(d4-d3)/2

##
df3 = all_dataframes[3]
xs_df3 = df3[:, 1]
ys_df3 = df3[:, 2]
using Optim
using Statistics
normal(x, μ, σ, A) = A/(σ*√(2π)) * exp(-1/2 * ((x-μ)/σ)^2)

function lossfunc(params, xs=xs_df3, ys=ys_df3)
    μ, σ, A = params
    
    predictions = normal.(xs, μ, σ, A)
    residuals = predictions - ys
    return sum(abs2, residuals)
end
sol = optimize(lossfunc, [mean(ys_df3), std(ys_df3), maximum(ys_df3)])
sol.minimizer
fitted_normal(x) = normal(x, sol.minimizer...)
begin
    fig, ax, _ = scatterlines(xs_df3, fitted_normal.(xs_df3), label="Fitted")
    scatterlines!(xs_df3, ys_df3, label="Data")
    ax.xlabel="Position [mm]"
    ax.ylabel="Intensity"
    axislegend(position=(1, 0))
    fig
end


##? Reading in spectra
data_itsl_spectrum = CSV.read(laserdiode_spectrum_datapath, DataFrame)
begin
    fig = Figure()
    Label(fig[0, 1], "Spectrum (top) and DCT of spectrum (bottom)\nfor \"Laserdiode HR\"", tellwidth=false)
    ax, _ = scatterlines(fig[1, 1], data_itsl_spectrum.lam, data_itsl_spectrum[:, " I"])
    ax.xlabel = "λ [nm]"
    ax.ylabel = "Intensity"
    ax, _ = scatterlines(fig[2, 1], dct(data_itsl_spectrum[:, " I"]))
    ax.xlabel="Index"
    ax.ylabel="DCT of intensity"
    fig|>display
    # save(joinpath(pic_dir, "Spectrum_DCT_Laserdiode_HR.png"), fig)
end
##
spectrum_filepaths = readdir(joinpath(my_dir, "Spectrum data"), join=true)
for i in (1:5)
    spectrum_filepath = spectrum_filepaths[i]
    df = CSV.read(spectrum_filepath, DataFrame, decimal=',', delim='\t', header=["lam", "I"])
    empty!(fig)
    Label(fig[0, 1], "Spectrum (top) and DCT of spectrum (bottom)\nfor \"$(splitext(splitpath(spectrum_filepath)[end])[1])\"", tellwidth=false)
    ax, _ = scatterlines(fig[1, 1], df.lam, df.I)
    ax.xlabel="λ [nm]"
    ax.ylabel="Intensity"
    ax, _ = scatterlines(fig[2, 1], dct(df.I))
    ax.xlabel="Index"
    ax.ylabel="DCT of intensity"
    # Makie.inline!(false)
    # DataInspector(fig)
    fig|>display
    # Makie.inline!(true)
    # save(joinpath(pic_dir, "Spectrum_DCT_$(splitext(splitpath(spectrum_filepath)[end])[1]).png"), fig)
end