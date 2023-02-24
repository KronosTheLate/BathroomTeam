using GLMakie; GLMakie.activate!(); Makie.inline!(false)
using Images
Axis = Makie.Axis  # Images also exports Axis
imdir = "/home/dennishb/Uni/Semester/8. Sem/Advanced Physical Optics/Lab excercise 3/Images 3 lab"  # Assumed "Pictures" is inside pwd, and contains the thorcam pictures.


let #! Loading pictures
    images_paths_all = readdir(imdir, join=true)
    image_paths_bmp = filter(p->occursin(".bmp", p), images_paths_all)
    global image_names = [splitext(splitpath(image_path)[end])[1] for image_path in image_paths_bmp]

    global images = Images.load.(image_paths_bmp)
    images = [image .|> Gray for image in images] # Convert to grayscale
end

camera_pixel_pitch = 6.784e-3 / 1280
SLM_pixel_pitch = 32e-6
λ = 632.8e-9
f = 300e-3
dist_to_period2(dist) = f*λ/dist

##! Task 1 - Linearity of "amplitude response function of SLM"
task1_mask = occursin.("grating2D_", image_names)
task1_images = images[task1_mask]
task1_imagemats = task1_images .|> x->Gray.(x) .|> x->gray.(x) .|> x->Float64.(x)
task1_imagenames = image_names[task1_mask]
task1_imagenames = replace(task1_imagenames) do name
    if occursin("rot", name)
        return name
    else
        return name * "_rot0"
    end
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

function rot_from_filename(filename)
    i_start = findall(==('o'), filename) |> only
    i_start += 1  # we find x, not p
    i_end = findnext(==('_'), filename, i_start)
    i_end === nothing  && (i_end=length(filename)+1)  # Handle case when lamX is final content in string
    return parse(Int64, filename[i_start+1:i_end-1])
end

pit_from_filename.(task1_imagenames)
rot_from_filename.(task1_imagenames)

let i=5 #? Locating centers. Max i is 10
    @show i
    image_name = task1_imagenames[i]
    @show image_name
    pit = pit_from_filename(image_name)
    rot = rot_from_filename(image_name)
    data = (pit=pit, rot=rot, image_name = image_name, centers=Vector{Float64}[])
    heatmap(task1_imagemats[i].|>log, axis=(title=image_name,))

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
    DataInspector()
    current_figure()
end


task1_imagecenters = [  # manually determined centers. First point is center, pairwise from there are opposite
    (pit = 7, rot = 0, image_name = "grating2D_graymax100_pitch7_rot0", centers = [[520.47, 652.54], [392.06, 655.34], [649.03, 651.24]]),
    (pit = 11, rot = 0, image_name = "grating2D_graymax255_pitch11_rot0", centers = [[520.8, 653.73], [602.75, 652.49], [438.7, 655.03], [359.77, 655.21], [683.87, 650.84], [766.35, 649.82], [275.42, 656.89], [195.81, 656.84], [846.24, 648.28], [932.24, 646.91], [108.91, 659.67]]),
    (pit = 11, rot = 30, image_name = "grating2D_graymax255_pitch11_rot30", centers = [[518.8, 628.54], [448.72, 573.12], [591.1, 685.27], [407.63, 672.31], [629.8, 584.11], [700.99, 640.51], [335.96, 614.77], [399.9, 820.27], [638.9, 434.11]]),
    (pit = 11, rot = 45, image_name = "grating2D_graymax255_pitch11_rot45", centers = [[519.23, 629.11], [403.6, 631.48], [635.97, 627.11], [464.76, 711.27], [577.36, 546.22], [460.3, 548.04], [579.91, 711.27], [521.62, 739.69], [516.84, 515.62], [287.86, 631.31], [751.76, 624.83]]),
    (pit = 15, rot = 0, image_name = "grating2D_graymax255_pitch15_rot0", centers = [[520.66, 652.45], [461.16, 653.05], [580.92, 651.74], [523.44, 735.49], [521.43, 552.11], [400.46, 653.31], [642.7, 650.07], [700.5, 649.24], [341.08, 654.52], [279.86, 654.9], [764.26, 647.66], [821.44, 647.66], [219.02, 656.76], [465.87, 734.82], [578.06, 506.67], [456.13, 507.17], [580.18, 736.32], [517.26, 422.26], [524.25, 797.09], [457.62, 423.49], [583.02, 798.02]]),#, [462.79, 799.05]
    (pit = 15, rot = 30, image_name = "grating2D_graymax255_pitch15_rot30", centers = [[519.19, 628.91], [466.11, 586.93], [571.89, 670.89], [437.79, 662.14], [601.97, 594.67], [651.91, 634.9], [386.67, 619.91]]),
    (pit = 15, rot = 45, image_name = "grating2D_graymax255_pitch15_rot45", centers = [[519.33, 629.53], [434.93, 630.82], [605.25, 626.84], [560.92, 568.1], [477.43, 690.84], [477.2, 570.15], [562.52, 688.8]]),
    (pit = 7, rot = 0, image_name = "grating2D_graymax255_pitch7_rot0", centers = [[520.92, 652.51], [392.16, 655.43], [649.56, 651.24], [779.13, 649.05], [263.57, 657.03], [136.28, 658.27], [907.87, 647.44]]),
    (pit = 7, rot = 30, image_name = "grating2D_graymax255_pitch7_rot30", centers = [[519.43, 628.44], [470.93, 490.45], [568.58, 765.05], [694.2, 558.19], [344.0, 698.13], [409.98, 538.76], [631.72, 717.01], [405.61, 592.9], [630.41, 656.76]]),
    (pit = 7, rot = 45, image_name = "grating2D_graymax255_pitch7_rot45", centers = [[519.97, 629.28], [337.99, 633.6], [702.39, 626.2], [517.13, 481.82], [521.5, 775.82], [339.09, 778.6], [700.43, 478.73], [336.24, 483.67], [704.14, 776.13]])
]

all_centers2 = getproperty.(task1_imagecenters, :centers)

distance_to_all_orders_px2 = [begin 
    centervec_pairs = []
    let centervec = all_centers2[i]
        @assert length(centervec) |> isodd
        for j in 2:2:length(centervec)
            push!(centervec_pairs, (centervec[j], centervec[j+1]))
        end
        centervec_pairs .|> identity
        pair_distances = [begin
            Δx = pair[1][1] - pair[2][1]
            Δy = pair[1][2] - pair[2][2]
            hypot(Δx, Δy)
        end for pair in centervec_pairs]
        avg_distances = pair_distances ./ 2
        output = sort(avg_distances)  # reverse to get inner (lower) orders first
        # output = first(output)  # only get inner order
        output
    end
end for i in eachindex(all_centers2)]

function distance_order2(n)
    [try
        getindex(distances, n)
    catch e
        NaN
    end for distances in distance_to_all_orders2]
end

begin
    estimated_periods_based_on_nearest_peaks = dist_to_period2.(distance_order2(1))
    true_periods = getproperty.(task1_imagecenters, :pit) .* SLM_pixel_pitch
    println("filenames")
    display(getproperty.(task1_imagecenters, :image_name))
    println("estimated_periods_based_on_nearest_peaks")
    display(estimated_periods_based_on_nearest_peaks)
    println("true_periods")
    display(true_periods)
end
##? Finding spatial period of source