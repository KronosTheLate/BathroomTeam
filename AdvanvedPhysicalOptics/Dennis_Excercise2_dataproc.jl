using GLMakie; GLMakie.activate!(); Makie.inline!(true)
using Images
Axis = Makie.Axis  # Images also exports Axis
using Statistics
using DataFrames
imdir = "/home/dennishb/Uni/Semester/8. Sem/Advanced Physical Optics/Lab excercise 2/Pictures"  # Assumed "Pictures" is inside pwd, and contains the thorcam pictures.

let #! Loading pictures
    images_paths_all = readdir(imdir, join=true)
    image_paths_bmp = filter(p->occursin(".bmp", p), images_paths_all)
    global image_names = [splitext(splitpath(image_path)[end])[1] for image_path in image_paths_bmp]

    global images = Images.load.(image_paths_bmp)
    images = [image .|> Gray for image in images] # Convert to grayscale
end


##! Task 3 - Linearity of "amplitude response function of SLM"
begin #? Masaging data
    task3_mask = occursin.("task3", image_names)  # Logical mask, used for indexing to get only pictures with "task3" in filename
    task3_uniform_grayscale_vals = image_names[task3_mask] .|> x->split(x, '_') |> x->x[2] .|> x->parse(Int64, x)
    sorted_inds = sortperm(task3_uniform_grayscale_vals)
    task3_uniform_grayscale_vals = task3_uniform_grayscale_vals[sorted_inds]
    task3_images = images[task3_mask][sorted_inds]
    task3_imagemats = task3_images .|> x->gray.(x) .|> float
end
task3_images[5]
task3_images[5] .|> Gray .|> gray .|> Float64 |> heatmap |> display
task3_images[5] .|> Gray .|> gray .|> Float64 |> rotr90 |> heatmap |> display
task3_images[5] .|> Gray .|> gray .|> Float64 |> x->permutedims(x, (2, 1)) |> heatmap |> display
task3_images[5] .|> Gray .|> gray .|> Float64 |> transpose |> heatmap |> display


with_theme(resolution=(1920÷2, 1080÷2), fontsize = 30) do#? Plotting
    means = mean.(task3_imagemats) .* 255
    medians = median.(task3_imagemats) .* 255
    # bla = extrema.(task3_imagemats)
    # mins, maxs = first.(bla)[sorted_inds], last.(bla)
    # stddevs = std.(task3_imagemats)
 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Uniform grayscal value of image", ylabel="Grayscale value measured")
    
    scatterlines!(task3_uniform_grayscale_vals, means, label="Mean")
    scatterlines!(task3_uniform_grayscale_vals, medians, label="Median")
    # scatterlines!(task3_uniform_grayscale_vals, mins, label="Minium")
    # scatterlines!(task3_uniform_grayscale_vals, maxs, label="Maximum")
    # scatterlines!(task3_uniform_grayscale_vals, stddevs, label="Std dev")
    # Legend(fig[1, 2], ax)
    axislegend(position=(1, 0.5))
    
    save(joinpath(homedir(), "Pictures", "SLM_linrange.png"), fig)
    fig |> display
end

##! Task 2 - diffraction pattern of blank SLM screen
camera_pixel_pitch = 6.784e-3 / 1280
SLM_pixel_pitch = 32e-6
λ = 632.8e-9
dist_to_period(dist, f) = f*λ/dist

begin #? Masaging data
    task2_mask = occursin.("task2", image_names)  # Logical mask, used for indexing to get only pictures with "task3" in filename
    task2_image = Gray.(only(images[task2_mask]))

    task2_imagemat = task2_image .|> gray .|> float
    heatmap(task2_imagemat)
    DataInspector()
    current_figure()
    task2_spotcenters_h =[[97, 646], [541, 654], [892, 663]]  # by manual determination
    Δx, Δy = task2_spotcenters_h[3] .- task2_spotcenters_h[1]
    dist_h = hypot(Δx, Δy)/2 * camera_pixel_pitch
    task2_spotcenters_v =[[522, 1064], [541, 654], [531, 254]]  # by manual determination
    Δx, Δy = task2_spotcenters_v[3] .- task2_spotcenters_v[1]
    dist_v = hypot(Δx, Δy)/2 * camera_pixel_pitch
    period_h = dist_to_period(dist_h)
    period_v = dist_to_period(dist_v)
    true_period = 32e-6
    period_h, period_v, true_period
end


##! Task 4 - diffraction pattern of blank SLM screen
function lam_from_filename(filename)
    i_start = findall(==('m'), filename)
    if length(i_start) == 0
        return NaN
    elseif length(i_start) > 2
        @warn "Found more than 1 'm' in filename. Returning NaN"
        return NaN
    end
    i_start = only(i_start)
    i_end = findnext(==('_'), filename, i_start)
    i_end === nothing  && (i_end=length(filename)+1)  # Handle case when lamX is final content in string
    return parse(Int64, filename[i_start+1:i_end-1])
end

function pit_from_filename(filename)
    i_start = findall(==('h'), filename)
    if length(i_start) == 0
        return NaN
    elseif length(i_start) > 2
        @warn "Found more than 1 'h' in filename. Returning NaN"
        return NaN
    end
    i_start = only(i_start)
    i_end = findnext(==('_'), filename, i_start)
    i_end === nothing  && (i_end=length(filename)+1)  # Handle case when lamX is final content in string
    return parse(Int64, filename[i_start+1:i_end-1])
end

function exp_from_filename(filename)
    i_start = findall(==('x'), filename) |> only
    i_start += 1  # we find x, not p
    i_end = findnext(==('_'), filename, i_start)
    i_end === nothing  && (i_end=length(filename)+1)  # Handle case when lamX is final content in string
    return filename[i_start+1:i_end-1]
end

let i=1 #? Locating centers
    task4_mask = occursin.("task4", image_names)  # Logical mask, used for indexing to get only pictures with "task3" in filename
    global task4_images = images[task4_mask] .|> x->Gray.(x)
    global task4_imagemats = task4_images  .|> x->gray.(x) .|> x->float.(x)
    # task4_imagemats = rotr90.(task4_imagemats)  # Correct, but messes up rest of code

    n_images = length(task4_images)
    @show i
    image_name = image_names[task4_mask][i]
    @show image_name
    lam = lam_from_filename(image_name)
    pit = pit_from_filename(image_name)
    exp_ = exp_from_filename(image_name)
    if exp_[1]=='0'
        exp_ = string(exp_[1], '.', exp_[2:end])
        exp_ = parse(Float64, exp_)
    else
        exp_ = parse(Int64, exp_)
    end
    rot = occursin("rot", image_name)
    data = (lam=lam, pit=pit, rot=rot, exp=exp_, image_name = image_name, centers=Vector{Float64}[])
    heatmap(task4_imagemats[i].|>log, axis=(title=image_name,))

    di = DataInspector()
    on(events(current_figure()).mousebutton) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                mp = mouseposition(current_axis()) |> Vector .|> Float64 .|> x->round(x, sigdigits=5)# events(current_figure()).mouseposition[]
                push!(data.centers, mp)
                println()
                println(data)
            end
        end
    end
    #! Note - log of intensity makes small peaks much clearer.
    current_figure()
end
# task4_imagemats[1] |> rotr90
task4_imagecenters_old = [  # manually determined centers, left to right
    (lam = 10, pit = NaN, rot = false, exp = 0.12, image_name = "task4_exp012_lam10",   centers = [[471.96, 715.82],[473.86, 610.25], [476.63, 500.07]]),
    (lam = 20, pit = NaN, rot = false, exp = 0.12, image_name = "task4_exp012_lam20",   centers = [[472.04, 664.1], [473.55, 611.38], [475.37, 551.56]]),
    (lam = 5, pit = NaN, rot = false, exp = 0.12, image_name = "task4_exp012_lam5",     centers = [[471.23, 762.8], [473.72, 610.34], [476.37, 458.22]]),
    (lam = 5, pit = NaN, rot = true, exp = 0.12, image_name = "task4_exp012_lam5_rot",  centers = [[112.1,  604.76],[291.6,  608.43], [473.71, 611.04], [655.67, 613.53], [838.61, 614.59]]),
    (lam = NaN, pit = 10, rot = false, exp = 0.12, image_name = "task4_exp012_pitch10", centers = [[472.08, 762.95],[471.4,  714.83], [473.84, 610.52], [476.73, 504.1], [476.61, 456.23]]),
    (lam = NaN, pit = 20, rot = false, exp = 0.12, image_name = "task4_exp012_pitch20", centers = [[473.3,  698.78],[472.87, 663.81], [474.08, 610.96], [476.05, 551.51], [475.69, 524.74]]),
    (lam = NaN, pit = 5, rot = false, exp = 0.12, image_name = "task4_exp012_pitch5",   centers = [[475.02, 785.41],[477.4,  634.43], [479.53, 481.17]]),
    (lam = NaN, pit = 5, rot = true, exp = 0.12, image_name = "task4_exp012_pitch5_rot",centers = [[109.39, 604.47],[291.83, 607.76], [474.37, 609.72], [655.5, 612.81], [836.67, 614.75]]),
    (lam = 5, pit = NaN, rot = false, exp = 10, image_name = "task4_exp10_lam5",        centers = [[471.51, 763.93],[473.96, 611.49], [476.44, 459.41]]),
    (lam = 5, pit = NaN, rot = true, exp = 10, image_name = "task4_exp10_lam5_rot",     centers = [[111.28, 605.66],[292.65, 608.73], [475.05, 612.44], [656.49, 613.54], [838.0, 616.48]]),
    (lam = NaN, pit = 20, rot = false, exp = 10, image_name = "task4_exp10_pitch20",    centers = [[473.02, 665.2], [474.66, 611.07], [474.84, 554.1]]),
    (lam = NaN, pit = 5, rot = false, exp = 10, image_name = "task4_exp10_pitch5",      centers = [[470.81, 761.87],[473.41, 610.34], [476.53, 457.55]]),
    (lam = NaN, pit = 5, rot = true, exp = 10, image_name = "task4_exp10_pitch5_rot",   centers = [[110.45, 606.49],[292.24, 609.06], [473.79, 611.63], [655.65, 614.05], [837.27, 615.51]])
]
task4_imagecenters = [  # manually determined centers, center first, pairwise from there
    (lam = 10, pit = NaN, rot = false, exp = 0.12, image_name = "task4_exp012_lam10",   centers = [[473.86, 610.25], [471.96, 715.82], [476.63, 500.07]]),
    (lam = 20, pit = NaN, rot = false, exp = 0.12, image_name = "task4_exp012_lam20",   centers = [[473.55, 611.38], [472.04, 664.1],  [475.37, 551.56]]),
    (lam = 5, pit = NaN, rot = false, exp = 0.12, image_name = "task4_exp012_lam5",     centers = [[473.72, 610.34], [471.23, 762.8],  [476.37, 458.22]]),
    (lam = 5, pit = NaN, rot = true, exp = 0.12, image_name = "task4_exp012_lam5_rot",  centers = [[473.71, 611.04], [291.6,  608.43], [655.67, 613.53], [112.1,  604.76], [838.61, 614.59]]),
    (lam = NaN, pit = 10, rot = false, exp = 0.12, image_name = "task4_exp012_pitch10", centers = [[473.84, 610.52], [471.4,  714.83], [476.73, 504.1],  [472.08, 762.95], [476.61, 456.23]]),
    (lam = NaN, pit = 20, rot = false, exp = 0.12, image_name = "task4_exp012_pitch20", centers = [[474.08, 610.96], [472.87, 663.81], [476.05, 551.51], [473.3,  698.78], [475.69, 524.74]]),
    (lam = NaN, pit = 5, rot = false, exp = 0.12, image_name = "task4_exp012_pitch5",   centers = [[477.4,  634.43], [475.02, 785.41], [479.53, 481.17]]),
    (lam = NaN, pit = 5, rot = true, exp = 0.12, image_name = "task4_exp012_pitch5_rot",centers = [[474.37, 609.72], [291.83, 607.76], [655.5, 612.81],  [109.39, 604.47], [836.67, 614.75]]),
    (lam = 5, pit = NaN, rot = false, exp = 10, image_name = "task4_exp10_lam5",        centers = [[473.96, 611.49], [471.51, 763.93], [476.44, 459.41]]),
    (lam = 5, pit = NaN, rot = true, exp = 10, image_name = "task4_exp10_lam5_rot",     centers = [[475.05, 612.44], [292.65, 608.73], [656.49, 613.54], [111.28, 605.66], [838.0, 616.48]]),
    (lam = NaN, pit = 20, rot = false, exp = 10, image_name = "task4_exp10_pitch20",    centers = [[474.66, 611.07], [473.02, 665.2],  [474.84, 554.1]]),
    (lam = NaN, pit = 5, rot = false, exp = 10, image_name = "task4_exp10_pitch5",      centers = [[473.41, 610.34], [470.81, 761.87], [476.53, 457.55]]),
    (lam = NaN, pit = 5, rot = true, exp = 10, image_name = "task4_exp10_pitch5_rot",   centers = [[473.79, 611.63], [292.24, 609.06], [655.65, 614.05], [110.45, 606.49], [837.27, 615.51]])
]
begin #! Making processed_centerdf
    processed_centerdf = DataFrame()
    processed_centerdf.centers_rel = map(task4_imagecenters) do nt
        vec2imag(v) = v[1] + v[2]*im
        centerloc = copy(nt.centers[1]) |> vec2imag
        return [vec2imag(center) .- centerloc for center in nt.centers[begin+1:end]]
    end
    processed_centerdf.centerdiffs = map(processed_centerdf.centers_rel) do center_vector
        [center_vector[i+1] - center_vector[i] for i in 1:2:length(center_vector)] ./ 2
    end
    processed_centerdf.distances = processed_centerdf.centerdiffs .|> x->abs.(x)
    processed_centerdf.angles = processed_centerdf.centerdiffs .|> x->rad2deg.(angle.(x))
    processed_centerdf = processed_centerdf .|> x->round.(x, digits=3)
end

processed_centerdf.centerdiffs
processed_centerdf.distances
processed_centerdf.angles

let  #! Making processed_data_px
    header = [:vec, :mag, :angle, :lam, :pit, :rot, :exp, :filename]
    mat = Matrix{Any}(undef, 0, length(header))
    for i in eachindex(task4_imagecenters)
        for j in eachindex(processed_centerdf.centerdiffs[i])
            datavec = hcat(
                processed_centerdf.centerdiffs[i][j],
                processed_centerdf.distances[i][j],
                processed_centerdf.angles[i][j],
                tuple(task4_imagecenters[i]...)[1:5]..., 
            )
            mat = vcat(mat, datavec)
        end
    end
    global processed_data_px = DataFrame(mat, header) .|> identity
end
processed_data_px |> vscodedisplay

processed_data_px_1order = filter(processed_data_px) do row
    count(==(row.filename), processed_data_px.filename) ≤ 1 ? true :
    row.mag == minimum(filter(x->x.filename==row.filename, processed_data_px).mag) ? true : false
end  # same as unique(processed_data_px, :filename), but more explicit to be sure
processed_data_px_1order |> vscodedisplay

relative_error(calc_val, true_val) = (calc_val - true_val)/true_val
let input_df = processed_data_px
    calculated_spatial_period_100 = dist_to_period.(input_df.mag .* camera_pixel_pitch, 100e-3)
    calculated_spatial_period_300 = dist_to_period.(input_df.mag .* camera_pixel_pitch, 300e-3)
    true_spatial_periods = [isnan(a) ? b : a for (a, b) in zip(
        input_df.pit .* SLM_pixel_pitch,
        input_df.lam .* SLM_pixel_pitch
    )] .|> x->round(x, sigdigits=5)

    global relerr100_all = relative_error.(calculated_spatial_period_100, true_spatial_periods) .* 100
    global relerr300_all = relative_error.(calculated_spatial_period_300, true_spatial_periods) .* 100
    df = DataFrame((; true_spatial_periods, calculated_spatial_period_100, calculated_spatial_period_300, relerr100_all, relerr300_all)) 
    df |> x->round.(x, sigdigits=3) |> vscodedisplay
end

let input_df = processed_data_px_1order
    calculated_spatial_period_100 = dist_to_period.(input_df.mag .* camera_pixel_pitch, 100e-3)
    calculated_spatial_period_300 = dist_to_period.(input_df.mag .* camera_pixel_pitch, 300e-3)
    true_spatial_periods = [isnan(a) ? b : a for (a, b) in zip(
        input_df.pit .* SLM_pixel_pitch,
        input_df.lam .* SLM_pixel_pitch
    )] .|> x->round(x, sigdigits=5)

    global relerr100_1st = relative_error.(calculated_spatial_period_100, true_spatial_periods) .* 100
    global relerr300_1st = relative_error.(calculated_spatial_period_300, true_spatial_periods) .* 100
    df = DataFrame((; true_spatial_periods, calculated_spatial_period_100, calculated_spatial_period_300, relerr100_1st, relerr300_1st)) 
    df |> x->round.(x, sigdigits=3) |> vscodedisplay
end

mean(abs, relerr100_all), mean(abs, relerr300_all)
mean(relerr100_all),           mean(relerr300_all)
median(relerr100_all),       median(relerr300_all)

mean(abs, relerr100_1st), mean(abs, relerr300_1st)
mean(relerr100_1st),           mean(relerr300_1st)
median(relerr100_1st),       median(relerr300_1st)


##¤ Diffraction efficiency
processed_data_px_1order

with_theme(fontsize=40) do
    overexposed_mask = [i ∈ [3, 4, 7, 8] ? true : false for i in 1:13] #map(x -> x.exp==0.12, task4_imagecenters)
    not_overexposed_centers = task4_imagecenters[overexposed_mask]
    not_overexposed_imagemats = task4_imagemats[overexposed_mask]
    not_overexposed_filenames = getproperty.(task4_imagecenters, :image_name)[overexposed_mask]

    all_centers = getproperty.(not_overexposed_centers, :centers)
    relevant_centers = [centervecs[1:3] for centervecs in all_centers]  # 0st order and two 1st orders
    radius = 20
    for i in eachindex(not_overexposed_centers)
        println(not_overexposed_filenames[i])
        imagemat = not_overexposed_imagemats[i]
        fig, ax1, _ = heatmap(imagemat)
        ax1.title="Linear scale"
        ax2, _ = heatmap(fig[1, 2], imagemat.|>log10)
        ax2.title = "Logarithmic scale"
        ax2.yaxisposition = :right
        # heatmap(imagemat.|>log10)
        sums = zeros(3)
        xmin = 0
        xmax = 0
        ymin = 0
        ymax = 0
        for (j, center) in enumerate(relevant_centers[i])
            # center = reverse(center)  # Rotated pictures
            box_inds = [(i, j) for i in -radius:radius, j in -radius:radius]
            circ_inds = filter(box_inds) do ind
                dist_from_center = hypot(ind...)
                dist_from_center ≤ radius ? true : false
            end
            sum_around_peak = sum(circ_inds) do ind
                imagemat[(tuple(round.(Int64, center)...) .+ ind)...]
            end
            sums[j] = sum_around_peak
            θs = range(0, 2π, 100)
            xs = radius .* cos.(θs) .+ center[1]
            ys = radius .* sin.(θs) .+ center[2]
            # xs, ys = ys, xs  # Rotated pictures
            xmin = min(xmin, minimum(xs))
            ymin = min(ymin, minimum(ys))
            xmax = max(xmax, maximum(xs))
            ymax = max(ymax, maximum(ys))
            lines!(ax1, xs, ys, color=:white, linewidth=0.8)
            lines!(ax2, xs, ys, color=:white, linewidth=0.8)
            # lines!([center[1]-radius, center[1]+radius], [center[2], center[2]], color=Cycled(j))
            # lines!([center[1], center[1]], [center[2]-radius, center[2]+radius], color=Cycled(j))
            println("Sum around $center: $sum_around_peak")
        end
        plotsidelength = max(xmax-xmin, ymax-ymin)
        current_cen = relevant_centers[i][1]
        blaval = 3
        xlims!(ax1, current_cen[1]-plotsidelength/blaval, current_cen[1]+plotsidelength/blaval)
        ylims!(ax1, current_cen[2]-plotsidelength/blaval, current_cen[2]+plotsidelength/blaval)
        xlims!(ax2, current_cen[1]-plotsidelength/blaval, current_cen[1]+plotsidelength/blaval)
        ylims!(ax2, current_cen[2]-plotsidelength/blaval, current_cen[2]+plotsidelength/blaval)
        # hypot((relevant_centers[1] .- relevant_centers[2])...)
        avg_1st_order = mean(sums[2:3])
        ratio = round(avg_1st_order/sums[1], sigdigits=5)
        @show ratio
        println()
        hidedecorations!.((ax1, ax2))
        colsize!(fig.layout, 1, Aspect(1, 1.))
        colsize!(fig.layout, 2, Aspect(1, 1.))
        filename = not_overexposed_filenames[i]
        Label(fig[0, :], (isnan(pit_from_filename(filename)) ? "Sinusoidal" : "Binary") * " grating" * ", Ratio = $ratio")
        resize_to_layout!(fig)
        Makie.save(joinpath(homedir(), "Pictures", "DiffractionEfficiency_$(filename).png"), fig)
        fig |> display
    end
end
##
1
#=
all_centers = getproperty.(task4_imagecenters, :centers)

distance_to_all_orders_px = [begin 
    centervec_pairs = []
    let centervec = all_centers[i]
        @assert length(centervec) |> isodd
        for j in 0:length(centervec)÷2-1
            push!(centervec_pairs, (centervec[begin+j], centervec[end-j]))
        end
        centervec_pairs .|> identity
        pair_distances = [begin
            Δx = pair[1][1] - pair[2][1]
            Δy = pair[1][2] - pair[2][2]
            hypot(Δx, Δy)
        end for pair in centervec_pairs]
        avg_distances = pair_distances ./ 2
        output = reverse(avg_distances)  # reverse to get inner (lower) orders first
        output
    end
end for i in eachindex(all_centers)]

distance_to_all_orders = distance_to_all_orders_px .* camera_pixel_pitch

function distance_order(n)
    [try
        getindex(distances, n)
    catch e
        NaN
    end for distances in distance_to_all_orders]
end

distance_order(1)
getproperty.(task4_imagecenters, :image_name)
=#