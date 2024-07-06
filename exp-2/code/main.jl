using FileIO: load
using Plots, FFTW, StatsBase, LibSndFile, Plots.Measures
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

# Calculating the peak magnitude for each angle
angles = round.(Int, LinRange(0, 360, 19))
times = LinRange(0, div(length(y), fs), 20)
peaks = zeros(19)
for i ∈ 1:19
    angle = angles[i]
    i_start = round(Int, times[i] * fs)
    i_end = round(Int, times[i+1] * fs)
    Y = y[1+i_start:i_end]

    # Calculating the Fast Fourier Transform of the signal and the magnitude in dB
    fft_result = fft(Y)
    freqaxis = fftfreq(length(Y), fs)
    N_half = div(length(Y), 2)
    magnitude_db = 20 * log10.(abs.(fft_result[2:N_half]))
    v_, i_ = findmax(magnitude_db)
    peaks[i] = v_
end

# Saving the data and the plot
CSV.write(DataDir * "processed/AngularMagnitude.csv", DataFrame(Angle=angles, Magnitude=peaks))

plot(angles, peaks, label="Peak Magnitude", title="Peak Magnitude vs. Angle",
    xlabel="Angle (Degrees)", ylabel="Magnitude (dB)", c=:black, lw=2,
    xticks=0:20:360, yticks=-20:5:15, legend=false, frame=:box, size=(800, 600), leftmargin=5mm)
savefig(FiguresDir * "AngularMagnitude.png")

############################################## Exprience 1 Part 2 ##############################################
# Loading the data
y, fs = load(DataDir * "processed/exno2.wav")
fs = Int(fs)
y = mean(y, dims=2)
y = y[y.!=0]

# Calculating the peak magnitude for each frequency
outfreq = 100 * (2 .^ collect(0:7))
times = LinRange(0, div(length(y), fs), 9)
peaks = zeros(8)
for i ∈ 1:8
    i_start = round(Int, times[i] * fs)
    i_end = round(Int, times[i+1] * fs)
    Y = y[1+i_start:i_end]

    # Calculating the Fast Fourier Transform of the signal and the magnitude in dB
    fft_result = fft(Y)
    freqaxis = fftfreq(length(Y), fs)
    N_half = div(length(Y), 2)
    magnitude_db = 20 * log10.(abs.(fft_result[2:N_half]))
    v_, i_ = findmax(magnitude_db)
    peaks[i] = v_
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
y_out = mean(y, dims=2)
y_out = y_out[y_out.!=0]

y, fs = load(DataDir * "processed/exno3-in.wav")
fs_in = Int(fs)
y_in = mean(y, dims=2)
y_in = y_in[y_in.!=0]

# Calculating the peak magnitude for each trial
times_in = LinRange(0, div(length(y_in), fs_in), 11)
times_out = LinRange(0, div(length(y_out), fs_out), 11)
peaks_in = zeros(10)
peaks_out = zeros(10)
for i ∈ 1:10
    i_start = round(Int, times_in[i] * fs_in)
    i_end = round(Int, times_in[i+1] * fs_in)
    Y = y_in[1+i_start:i_end]
    fft_result = fft(Y)
    freqaxis = fftfreq(length(Y), fs)
    N_half = div(length(Y), 2)
    magnitude_db = 20 * log10.(abs.(fft_result[2:N_half]))
    v_, i_ = findmax(magnitude_db)
    peaks_in[i] = deepcopy(v_)

    i_start = round(Int, times_out[i] * fs_out)
    i_end = round(Int, times_out[i+1] * fs_out)
    Y = y_out[1+i_start:i_end]
    fft_result = fft(Y)
    freqaxis = fftfreq(length(Y), fs)
    N_half = div(length(Y), 2)
    magnitude_db = 20 * log10.(abs.(fft_result[2:N_half]))
    v_, i_ = findmax(magnitude_db)
    peaks_out[i] = deepcopy(v_)
end

# Saving the data and the plot
CSV.write(DataDir * "processed/InboundMagnitude.csv", DataFrame(Input=peaks_in, Output=peaks_out))

plot(peaks_in, xlims=(0.5, 10.5), label="Input audio", c=:black, lw=2)
plot!(peaks_out, xlims=(0.5, 10.5), xticks=1:10, label="Output audio", title="Peak Magnitude for different trials", xlabel="Trial",
    ylabel="Magnitude (dB)", c=:blue, lw=2, legend=200, frame=:box, size=(800, 600))
savefig(FiguresDir * "AmplitudeLinearity.png")