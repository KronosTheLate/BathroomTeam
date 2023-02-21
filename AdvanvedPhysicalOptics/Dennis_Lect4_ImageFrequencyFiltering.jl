using DelimitedFiles
using Images
using FFTW
using GLMakie

cd(joinpath(homedir(), "Uni", "Semester", "8. Sem", "Advanced Physical Optics", "BathroomTeamGithubFolder"))
begin
    # image_matrix = readdlm("lena.dat") |> x->reverse(x, dims=1)
    image_matrix = Images.load("Vilde's pus, nå med farger!.png") .|> Gray .|> Float64
    image_matrix = image_matrix ./ maximum(image_matrix)
    image_matrix .|> Gray
end

##¤ Coupld of functions
mynorm(A) = A ./ maximum(A)
viz(A) = display(mynorm(A) .|> Gray)

##¤ High pass filtering
for n_freqs_pass = range(5, step=25, stop=200)
    input_image = image_matrix
    image_dft = fft(input_image)|>fftshift
    zerofreq_inds = size(image_dft) .÷ 2

    output_dft = ComplexF64[(abs(i - zerofreq_inds[1]) < n_freqs_pass && abs(j - zerofreq_inds[2]) < n_freqs_pass ? 0 : image_dft[i, j]) for i in axes(image_dft, 1), j in axes(image_dft, 2)]

    show1 = input_image |> mynorm
    show2 = image_dft .|>abs .|> log10 |> mynorm
    show3 = output_dft .|> abs .|>log10 |> mynorm
    show4 = ifft(output_dft|>fftshift) .|>abs |> mynorm
    [show1 show2; show3 show4] |> viz
end

##¤ Low pass filtering
for n_freqs_pass = range(5, step=25, stop=200) |> reverse
    input_image = image_matrix
    image_dft = fft(input_image)|>fftshift
    zerofreq_inds = size(image_dft) .÷ 2

    output_dft = ComplexF64[(abs(i - zerofreq_inds[1]) > n_freqs_pass || abs(j - zerofreq_inds[2]) > n_freqs_pass ? 0 : image_dft[i, j]) for i in axes(image_dft, 1), j in axes(image_dft, 2)]

    show1 = input_image |> mynorm
    show2 = image_dft .|>abs .|> log10 |> mynorm
    show3 = output_dft .|> abs .|>log10 |> mynorm
    show4 = ifft(output_dft|>fftshift) .|>abs |> mynorm
    [show1 show2; show3 show4] |> viz
end