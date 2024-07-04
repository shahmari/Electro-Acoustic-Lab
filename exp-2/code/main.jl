using FileIO: load
using Plots, FFTW, StatsBase, LibSndFile
using CSV, DataFrames

DataDir = "../data/" # Directory of the data
FiguresDir = "../figure/" # Directory of the figures
cd(dirname(@__FILE__)) # Changing the current directory to the directory of the file

############################################## Exprience 1 Part 1 ##############################################
# Loading the data
y, fs = load(DataDir * "processed/exno1.wav")
fs = Int(fs)
y = mean(y, dims=2)
y = y[y.!=0]

# Calculating the Fast Fourier Transform of the signal and the magnitude in dB
fft_result = fft(y)
freqaxis = fftfreq(length(y), fs)
magnitude_db = 20 * log10.(abs.(fft_result))

# Calculating the peak magnitude for each angle
angles = round.(Int, LinRange(0, 360, 20))
times = LinRange(0, div(length(y), fs), 21)
peaks = zeros(20)
for i ∈ 1:20
    angle = angles[i]
    i_start = round(Int, times[i] * fs)
    i_end = round(Int, times[i+1] * fs)
    ansemble = magnitude_db[1+i_start:i_end]
    histfit = fit(Histogram, ansemble, nbins=1000)
    peak = histfit.edges[1][argmax(histfit.weights)]
    peaks[i] = peak
end

# Saving the data and the plot
CSV.write(DataDir * "processed/AngularMagnitude.csv", DataFrame(Angle=angles, Magnitude=peaks))

plot(angles, peaks, label="Peak Magnitude", title="Peak Magnitude vs. Angle",
    xlabel="Angle (Degrees)", ylabel="Magnitude (dB)", c=:black, lw=2,
    xticks=0:20:360, yticks=-20:5:15, legend=false, frame=:box, size=(800, 600))
savefig(FiguresDir * "AngularMagnitude.png")

# Creating a gif of the histogram of the magnitude
# lensec = round(Int, length(y) / fs) - 30
# @gif for i ∈ 1:10:lensec
#     histogram(magnitude_db[i*fs:i*fs+30*fs],
#         xticks=(-100:10:100), normed=true, xlims=(-40, 70), ylims=(-0.005, 0.1), c=:black, fill=true,
#         label="Time: $i s", xlabel="Magnitude (dB)", ylabel="Density")
# end

############################################## Exprience 1 Part 2 ##############################################
# Loading the data
y, fs = load(DataDir * "processed/exno2.wav")
fs = Int(fs)
y = mean(y, dims=2)
y = y[y.!=0]
# Calculating the Fast Fourier Transform of the signal and the magnitude in dB
fft_result = fft(y)
freqaxis = fftfreq(length(y), fs)
magnitude_db = 20 * log10.(abs.(fft_result))

# Calculating the peak magnitude for each frequency
outfreq = 100 * (2 .^ collect(0:7))
times = LinRange(0, div(length(y), fs), 9)
peaks = zeros(8)
for i ∈ 1:8
    i_start = round(Int, times[i] * fs)
    i_end = round(Int, times[i+1] * fs)
    ansemble = magnitude_db[1+i_start:i_end]
    histfit = fit(Histogram, ansemble, nbins=1000)
    peak = histfit.edges[1][argmax(histfit.weights)]
    peaks[i] = peak
end

# Saving the data and the plot
CSV.write(DataDir * "processed/FrequentialMagnitude.csv", DataFrame(Frequency=outfreq, Magnitude=peaks))
plot(outfreq, peaks, xaxis=:log2, xticks=(outfreq, string.(outfreq)),
    label="Peak Magnitude", title="Peak Magnitude vs. Frequency", xlabel="Frequency (Hz)",
    ylabel="Magnitude (dB)", c=:black, lw=2, legend=false, frame=:box, size=(800, 600))
savefig(FiguresDir * "FrequentialMagnitude.png")

############################################## Exprience 1 Part 3 ##############################################
# Loading the data and Calculating the Fast Fourier Transform of the signal and the magnitude in dB
y, fs = load(DataDir * "processed/exno3-out.wav")
fs_out = Int(fs)
y = mean(y, dims=2)
y = y[y.!=0]
fft_result = fft(y)
freqaxis = fftfreq(length(y), fs)
magnitude_db_out = 20 * log10.(abs.(fft_result))

y, fs = load(DataDir * "processed/exno3-in.wav")
fs_in = Int(fs)
y = mean(y, dims=2)
y = y[y.!=0]
fft_result = fft(y)
freqaxis = fftfreq(length(y), fs)
magnitude_db_in = 20 * log10.(abs.(fft_result))

# Calculating the peak magnitude for each trial
times_in = LinRange(0, div(length(magnitude_db_in), fs_in), 11)
times_out = LinRange(0, div(length(magnitude_db_out), fs_out), 11)
peaks_in = zeros(10)
peaks_out = zeros(10)
for i ∈ 1:10
    i_start = round(Int, times_in[i] * fs_in)
    i_end = round(Int, times_in[i+1] * fs_in)
    ansemble = magnitude_db_in[1+i_start:i_end]
    histfit = fit(Histogram, ansemble, nbins=1000)
    peak = histfit.edges[1][argmax(histfit.weights)]
    peaks_in[i] = deepcopy(peak)

    i_start = round(Int, times_out[i] * fs_out)
    i_end = round(Int, times_out[i+1] * fs_out)
    ansemble = magnitude_db_out[1+i_start:i_end]
    histfit = fit(Histogram, ansemble, nbins=1000)
    peak = histfit.edges[1][argmax(histfit.weights)]
    peaks_out[i] = deepcopy(peak)
end

# Saving the data and the plot
CSV.write(DataDir * "processed/InboundMagnitude.csv", DataFrame(Input=peaks_in, Output=peaks_out))

plot(peaks_in, xlims=(0.5, 10.5), label="Input audio", c=:black, lw=2)
plot!(peaks_out, xlims=(0.5, 10.5), xticks=1:10, label="Output audio", title="Peak Magnitude for different trials", xlabel="Trial",
    ylabel="Magnitude (dB)", c=:blue, lw=2, legend=200, frame=:box, size=(800, 600))
savefig(FiguresDir * "AmplitudeLinearity.png")