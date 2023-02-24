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
f = 100e-3
dist_to_period(dist) = f*λ/dist

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

task4_imagecenters_old = [  # manually determined centers. If rot!=0, points are pairwise opposite, and 0 order is last center
    (pit = 7, rot = 0, image_name = "grating2D_graymax100_pitch7_rot0", centers = [[392.31, 655.0], [520.9, 653.23], [649.01, 651.1]]),
    (pit = 11, rot = 0, image_name = "grating2D_graymax255_pitch11_rot0", centers = [[109.17, 658.95], [196.46, 657.66], [275.13, 656.94], [359.22, 655.95], [438.58, 654.52], [520.29, 653.06], [602.5, 652.09], [683.48, 651.4], [766.57, 649.23], [846.39, 647.59], [932.17, 646.91]]),
    (pit = 11, rot = 30, image_name = "grating2D_graymax255_pitch11_rot30", centers = [[399.55, 821.18], [639.47, 435.52], [336.55, 615.59], [702.07, 642.16], [630.42, 584.68], [407.12, 672.19], [448.47, 572.68], [591.46, 685.52], [518.49, 628.04]]),
    (pit = 11, rot = 45, image_name = "grating2D_graymax255_pitch11_rot45", centers = [[403.65, 631.39], [635.71, 628.42], [517.61, 519.02], [522.31, 740.61], [464.2, 711.73], [576.79, 547.06], [460.72, 548.9], [580.14, 711.16], [518.82, 629.61]]),
    (pit = 15, rot = 0, image_name = "grating2D_graymax255_pitch15_rot0", centers = [[159.3, 657.95], [220.01, 658.31], [279.69, 656.68], [341.55, 655.6], [400.94, 654.34], [461.4, 653.42], [521.22, 653.11], [580.74, 652.48], [642.3, 651.68], [701.13, 649.93], [764.27, 649.03], [822.34, 649.15], [883.66, 648.38], [516.62, 422.37], [518.59, 507.86], [521.28, 553.54], [523.96, 736.32], [524.02, 797.4], [524.34, 880.81], [466.17, 736.26], [580.35, 542.44], [577.89, 507.59], [465.83, 735.85], [582.27, 797.28], [456.66, 508.22]]),
    (pit = 15, rot = 30, image_name = "grating2D_graymax255_pitch15_rot30", centers = [[635.34, 443.0], [402.53, 814.16], [386.13, 620.9], [653.29, 636.97], [601.91, 594.52], [438.49, 660.41], [466.59, 586.14], [572.45, 670.84], [519.24, 628.89]]),
    (pit = 15, rot = 45, image_name = "grating2D_graymax255_pitch15_rot45", centers = [[434.16, 631.29], [604.94, 628.76], [563.16, 689.14], [476.32, 570.42], [561.04, 568.78], [478.02, 691.6], [518.66, 629.44]]),
    (pit = 7, rot = 0, image_name = "grating2D_graymax255_pitch7_rot0", centers = [[136.39, 658.77], [264.06, 657.86], [391.96, 655.88], [521.51, 653.53], [649.4, 651.55], [779.26, 650.18], [908.4, 647.56]]),
    (pit = 7, rot = 30, image_name = "grating2D_graymax255_pitch7_rot30", centers = [[394.86, 836.2], [645.28, 422.06], [470.58, 491.28], [569.04, 766.26], [695.28, 559.32], [344.31, 699.21], [519.09, 628.72]]),
    (pit = 7, rot = 45, image_name = "grating2D_graymax255_pitch7_rot45", centers = [[154.07, 635.64], [886.02, 624.98], [700.78, 479.55], [339.63, 778.88], [522.23, 777.48], [517.27, 482.1], [334.61, 485.29], [705.62, 774.22], [703.07, 627.08], [337.34, 632.64], [519.33, 631.14]])
]
task4_imagecenters = [  # manually determined centers. First point is center, pairwise from there are opposite
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

all_centers2 = getproperty.(task4_imagecenters, :centers)

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

distance_to_all_orders2 = distance_to_all_orders_px2 .* camera_pixel_pitch

function distance_order2(n)
    [try
        getindex(distances, n)
    catch e
        NaN
    end for distances in distance_to_all_orders2]
end
getproperty.(task4_imagecenters, :image)
dist_to_period.(distance_order2(1))
true_periods = getproperty.(task4_imagecenters, :pit) .* SLM_pixel_pitch
##? Finding spatial period of source