using GLMakie; GLMakie.activate!(); Makie.inline!(false)
using Images
Axis = Makie.Axis  # Images also exports Axis
using Statistics
imdir = joinpath(pwd(), "Pictures")  # Assumed "Pictures" is inside pwd, and contains the thorcam pictures.


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

let #? Plotting
    means = mean.(task3_imagemats) .* 255
    medians = median.(task3_imagemats) .* 255
    # bla = extrema.(task3_imagemats)
    # mins, maxs = first.(bla)[sorted_inds], last.(bla)
    # stddevs = std.(task3_imagemats)

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Uniform grayscal value", ylabel="Grayscale value")
    
    scatterlines!(task3_uniform_grayscale_vals, means, label="Mean")
    scatterlines!(task3_uniform_grayscale_vals, medians, label="Median")
    # scatterlines!(task3_uniform_grayscale_vals, mins, label="Minium")
    # scatterlines!(task3_uniform_grayscale_vals, maxs, label="Maximum")
    # scatterlines!(task3_uniform_grayscale_vals, stddevs, label="Std dev")
    Legend(fig[1, 2], ax)
    
    fig |> display
end

##! Task 2 - diffraction pattern of blank SLM screen
camera_pixel_pitch = 6.784e-3 / 1280
SLM_pixel_pitch = 32e-6
λ = 632.8e-9
f = 100e-3
dist_to_period(dist) = f*λ/dist
begin #? Masaging data
    task2_mask = occursin.("task2", image_names)  # Logical mask, used for indexing to get only pictures with "task3" in filename
    task2_image = Gray.(only(images[task2_mask]))
    
    #? Filtering approch did not work
    # filtered_image = task2_image
    # for _ in 1:100
    #     filtered_image = imfilter(filtered_image, ImageFiltering.Kernel.gaussian(2), NA())
    # end
    # findlocalmaxima(filtered_image)

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
    task4_images = images[task4_mask] .|> x->Gray.(x)
    task4_imagemats = task4_images  .|> x->gray.(x) .|> x->float.(x)

    n_images = length(task4_images)
    # task4_spotcenters = Vector{Any}(undef, n_images)

    # task4_spotcenters = [
    #     (lam=10, rot=false,
    #     centers=[[472, 716]
    #              [473, 610]
    #              [476, 500]]),
    #     (lam=10, rot=false,
    #     centers=[[472, 716]
    #              [473, 610]
    #              [476, 500]])
    # ]

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

task4_imagecenters = [  # manually determined centers
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

distance_order_1 = distance_order(1)
distance_order(2)
getproperty.(task4_imagecenters, :image_name)