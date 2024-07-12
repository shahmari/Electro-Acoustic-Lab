using FileIO: load
using Plots, Plots.Measures, FFTW, StatsBase, LibSndFile
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
model(x, p) = p[1] .+ p[2] * x
Y = 1 ./ heights
X_1 = harmonics[:, 1]
X_2 = harmonics[:, 2]
X_3 = harmonics[:, 3]
X_4 = harmonics[:, 4]

fit_1 = curve_fit(model, X_1, Y, [1.0, 1.0])
fit_2 = curve_fit(model, X_2, Y, [1.0, 1.0])
fit_3 = curve_fit(model, X_3, Y, [1.0, 1.0])
fit_4 = curve_fit(model, X_4, Y, [1.0, 1.0])

plot(X_1, model(X_1, fit_1.param), ribbon=stderror(fit_1)[1], c=:darkblue, label="1st: y = a + b.x, a = $(round(fit_1.param[1], digits=6)), b = $(round(fit_1.param[2], digits=6))")
plot!(X_1, Y, label=nothing, ls=:dash, c=:black, xerr=errs[:, 1])
plot!(X_2, model(X_2, fit_2.param), ribbon=stderror(fit_2)[1], c=:darkred, label="2nd: y = a + b.x, a = $(round(fit_2.param[1], digits=6)), b = $(round(fit_2.param[2], digits=6))")
plot!(X_2, Y, label=nothing, ls=:dash, c=:black, xerr=errs[:, 2])
plot!(X_3, model(X_3, fit_3.param), ribbon=stderror(fit_3)[1], c=:darkgreen, label="3rd: y = a + b.x, a = $(round(fit_3.param[1], digits=6)), b = $(round(fit_3.param[2], digits=6))")
plot!(X_3, Y, label=nothing, ls=:dash, c=:black, xerr=errs[:, 3])
plot!(X_4, model(X_4, fit_4.param), ribbon=stderror(fit_4)[1], c=:darkorange, label="4th: y = a + b.x, a = $(round(fit_4.param[1], digits=6)), b = $(round(fit_4.param[2], digits=6))")
plot!(X_4, Y, label=nothing, ls=:dash, c=:black, xerr=errs[:, 4])
plot!(ylabel="Height⁻¹ (1/cm)", xlabel="Frequency (Hz)", title="Frequency vs Height⁻¹ - Straw", legend=:top, frame=:box, ylims=(0.09, 0.24), size=(1000, 400), leftmargin=5mm, bottommargin=5mm)
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