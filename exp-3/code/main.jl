using FileIO: load
using Plots, FFTW, StatsBase, LibSndFile
using CSV, DataFrames, LsqFit

DataDir = "../data/" # Directory of the data
FiguresDir = "../figure/" # Directory of the figures
cd(dirname(@__FILE__)) # Changing the current directory to the directory of the file

heights = [10, 9, 8, 7, 6, 5]
trialtimes = [0, 5, 10]

################################################################ Exp 3 - Part 1 ################################################################
function GetFrequency(Y, fs; window_size=100)
    fft_result = fft(Y)
    freqs = fftfreq(length(Y), fs)
    N_half = div(length(Y), 2)
    magnitude_db = 20 * log10.(abs.(fft_result))
    MA_freqs = [mean(freqs[i:i+window_size]) for i in 1:length(freqs[2:N_half])-window_size]
    MA_mags = [mean(magnitude_db[i:i+window_size]) for i in 1:length(magnitude_db[2:N_half])-window_size]
    # finding the first harmonics:
    i_1 = argmax(MA_mags)
    i_2 = argmax(MA_mags[2*i_1:4*i_1]) + 2 * i_1
    i_3 = argmax(MA_mags[4*i_1:6*i_1]) + 4 * i_1
    i_4 = argmax(MA_mags[6*i_1:8*i_1]) + 6 * i_1
    return [MA_freqs[i_1], MA_freqs[i_2], MA_freqs[i_3], MA_freqs[i_4]]
end

Freqs = zeros(6, 4, 3)

for i in 1:6
    y, fs = load(DataDir * "processed/Exp3-P1-T$i.wav")
    fs = Int(fs)
    y = y .+ eps()
    Y1 = y[1+fs*trialtimes[1]:fs*trialtimes[2]]
    Y2 = y[1+fs*trialtimes[2]:fs*trialtimes[3]]
    Y3 = y[1+fs*trialtimes[3]:end]

    Freqs[i, :, 1] = GetFrequency(Y1, fs)
    Freqs[i, :, 2] = GetFrequency(Y2, fs)
    Freqs[i, :, 3] = GetFrequency(Y3, fs)
end

errs = std(Freqs, dims=3)[:, :, 1]
harmonics = mean(Freqs, dims=3)[:, :, 1]

# Fitting the data to a logarithmic model
model(x, p) = p[1] .+ p[2] * log.(x)
X = heights
Y_1 = harmonics[:, 1]
Y_2 = harmonics[:, 2]
Y_3 = harmonics[:, 3]
Y_4 = harmonics[:, 4]

fit_1 = curve_fit(model, X, Y_1, [1.0, 1.0])
fit_2 = curve_fit(model, X, Y_2, [1.0, 1.0])
fit_3 = curve_fit(model, X, Y_3, [1.0, 1.0])
fit_4 = curve_fit(model, X, Y_4, [1.0, 1.0])

plot(heights, model(heights, fit_1.param), c=:darkblue, ribbon=stderror(fit_1),
    label="1st: y = a + b×log(x), a = $(round(fit_1.param[1], digits=2)), b = $(round(fit_1.param[2], digits=2))")
plot!(heights, harmonics[:, 1], label=nothing, ls=:dash, c=:black, yerr=errs[:, 1])
plot!(heights, model(heights, fit_2.param), c=:darkred, ribbon=stderror(fit_2),
    label="2nd: y = a + b×log(x), a = $(round(fit_2.param[1], digits=2)), b = $(round(fit_2.param[2], digits=2))")
plot!(heights, harmonics[:, 2], label=nothing, ls=:dash, c=:black, yerr=errs[:, 2])
plot!(heights, model(heights, fit_3.param), c=:darkgreen, ribbon=stderror(fit_3),
    label="3rd: y = a + b×log(x), a = $(round(fit_3.param[1], digits=2)), b = $(round(fit_3.param[2], digits=2))")
plot!(heights, harmonics[:, 3], label=nothing, ls=:dash, c=:black, yerr=errs[:, 3])
plot!(heights, model(heights, fit_4.param), c=:darkorange, ribbon=stderror(fit_4),
    label="4th: y = a + b×log(x), a = $(round(fit_4.param[1], digits=2)), b = $(round(fit_4.param[2], digits=2))")
plot!(heights, harmonics[:, 4], label=nothing, ls=:dash, c=:black, yerr=errs[:, 4])
plot!(xlabel="Height (cm)", ylabel="Frequency (Hz)", title="Frequency vs Height - Straw", legend=:topright,
    yscale=:log10, ylims=(500, 25000), yticks=exp10.(2:0.5:4), frame=:box)
savefig(FiguresDir * "FrequencyVsHeight.png")

# Saving the data
df = DataFrame(Height=heights, Harmonic1=harmonics[:, 1], Harmonic2=harmonics[:, 2],
    Harmonic3=harmonics[:, 3], Harmonic4=harmonics[:, 4], Error1=errs[:, 1], Error2=errs[:, 2],
    Error3=errs[:, 3], Error4=errs[:, 4])
CSV.write(DataDir * "processed/Exp3-P1.csv", df)

################################################################ Exp 3 - Part 2 ################################################################
volumes = [50, 100, 150, 200, 250]
trialtimes = [0, 5, 10]

function GetFrequency(Y, fs; window_size=100)
    fft_result = fft(Y)
    freqs = fftfreq(length(Y), fs)
    N_half = div(length(Y), 2)
    magnitude_db = 20 * log10.(abs.(fft_result))
    MA_freqs = [mean(freqs[i:i+window_size]) for i in 1:length(freqs[2:N_half])-window_size]
    MA_mags = [mean(magnitude_db[i:i+window_size]) for i in 1:length(magnitude_db[2:N_half])-window_size]
    # finding the first harmonics:
    i_1 = argmax(MA_mags)
    i_2 = argmax(MA_mags[5*i_1:8*i_1]) + 5 * i_1
    i_3 = argmax(MA_mags[8*i_1:12*i_1]) + 8 * i_1
    i_4 = argmax(MA_mags[12*i_1:16*i_1]) + 12 * i_1
    return [MA_freqs[i_1], MA_freqs[i_2], MA_freqs[i_3], MA_freqs[i_4]]
end

Freqs = zeros(5, 4, 3)

for i in 1:5
    y, fs = load(DataDir * "processed/Exp3-P2-T$i.wav")
    fs = Int(fs)
    y = y .+ eps()
    Y1 = y[1+fs*trialtimes[1]:fs*trialtimes[2]]
    Y2 = y[1+fs*trialtimes[2]:fs*trialtimes[3]]
    Y3 = y[1+fs*trialtimes[3]:end]

    Freqs[i, :, 1] = GetFrequency(Y1, fs)
    Freqs[i, :, 2] = GetFrequency(Y2, fs)
    Freqs[i, :, 3] = GetFrequency(Y3, fs)
end

errs = std(Freqs, dims=3)[:, :, 1]
harmonics = mean(Freqs, dims=3)[:, :, 1]

# Fitting the data to a logarithmic model

plot(volumes, harmonics[:, 1], label=nothing, ls=:dash, c=:blue, yerr=errs[:, 1])
plot!(volumes, harmonics[:, 2], label=nothing, ls=:dash, c=:red, yerr=errs[:, 2])
plot!(volumes, harmonics[:, 3], label=nothing, ls=:dash, c=:green, yerr=errs[:, 3])
plot!(volumes, harmonics[:, 4], label=nothing, ls=:dash, c=:orange, yerr=errs[:, 4])
plot!(xlabel="Volume (cc)", ylabel="Frequency (Hz)", title="Frequency vs Volume - Bottle", legend=:topright, frame=:box)
savefig(FiguresDir * "FrequencyVsVolume.png")

# Saving the data
df = DataFrame(Volume=volumes, Harmonic1=harmonics[:, 1], Harmonic2=harmonics[:, 2],
    Harmonic3=harmonics[:, 3], Harmonic4=harmonics[:, 4], Error1=errs[:, 1], Error2=errs[:, 2],
    Error3=errs[:, 3], Error4=errs[:, 4])
CSV.write(DataDir * "processed/Exp3-P2.csv", df)