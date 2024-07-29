using FileIO: load
using Plots, FFTW, LibSndFile
using CSV, DataFrames, LsqFit
using StatsBase: mean, std

DataDir = "../data/" # Directory of the data
FiguresDir = "../figure/" # Directory of the figures
cd(dirname(@__FILE__)) # Changing the current directory to the directory of the file

############################################## Exprience 6 Part 1 ##############################################
y, fs = load(DataDir * "Exp6-P1.wav")
fs = Int(fs)
y = mean(y, dims=2)
y = y .+ eps()

Heights = collect(1:10)
TLen = 1.5 # Length of each trial
NumTrials = 10 # Number of trials
H1FE = 1505.0
H2FE = 3850.0
PeakMags1 = []
PeakFreqs1 = []
PeakMags2 = []
PeakFreqs2 = []
for i ∈ 1:NumTrials
    Y = y[(i-1)*Int(TLen * fs)+1:i*Int(TLen * fs)]
    fft_result = fft(Y)
    freqaxis = fftfreq(length(Y), fs)
    N_half = div(length(Y), 2)
    magnitude_db = 20 * log10.(abs.(fft_result))
    v_, i_ = findmax(magnitude_db[1000:3000])
    push!(PeakMags1, v_)
    push!(PeakFreqs1, freqaxis[i_+1000])
    v_, i_ = findmax(magnitude_db[3000:6000])
    push!(PeakMags2, v_)
    push!(PeakFreqs2, freqaxis[i_+3000])
end

############################ Plots:
############################ Frequencies
plot(Heights, PeakFreqs1, m=:c, label="First Harmonic Frequencies", c=:steelblue)
plot!(Heights, PeakFreqs2, m=:s, label="Second Harmonic Frequencies", c=:purple)
plot!(title="Frequency of Tones, Percussion", frame=:box, xlabel="Height (CM)", ylabel="Frequency (Hz)", legend=:right)
savefig(FiguresDir * "Exp6-P1-Freqs.png")

############################ First Harmonic
X = (((H1FE ./ PeakFreqs1) .^ 2) .- 1) .^ 0.25
Y = Heights
p0 = [12.0, 10.0]
@. model(x, p) = p[1] - (p[1] / (p[2]^0.25)) * x
fit = curve_fit(model, X, Y, p0)
fit_err = stderror(fit)
p = coef(fit)
H☆ = p[1]
β = p[2]
Heights_err = 0.2
H☆_err = fit_err[1]
β_err = fit_err[2]
X_err = begin
    @. model(x, p) = p[1] + p[2] * x
    fit = curve_fit(model, X, Y, p0)
    p = coef(fit)
    fit_err = stderror(fit)
    X .* sqrt((fit_err[1] / p[1])^2 + (fit_err[2] / p[2])^2)
end
scatter(X, Y, label="Data", c=:red, yerr=fill(Heights_err, length(X)), xerr=X_err, ms=2)
plot!(0.2:0.01:1.25, model(0.2:0.01:1.25, p), c=:blue, label="H* = $(round(H☆, digits=2)),\nβ=$(round(β, digits=2))", xlabel="X", ylabel="Heights (cm)", legend=:bottomleft, frame=:box)
savefig(FiguresDir * "Exp6-P1-H1.png")

############################ Second Harmonic
X = (((H2FE ./ PeakFreqs2) .^ 2) .- 1) .^ 0.25
Y = Heights
p0 = [12.0, 10.0]
@. model(x, p) = p[1] - (p[1] / (p[2]^0.25)) * x
fit = curve_fit(model, X, Y, p0)
fit_err = stderror(fit)
p = coef(fit)
H☆ = p[1]
β = p[2]
Heights_err = 0.2
H☆_err = fit_err[1]
β_err = fit_err[2]
X_err = begin
    @. model(x, p) = p[1] + p[2] * x
    fit = curve_fit(model, X, Y, p0)
    p = coef(fit)
    fit_err = stderror(fit)
    X .* sqrt((fit_err[1] / p[1])^2 + (fit_err[2] / p[2])^2)
end
scatter(X, Y, label="Data", c=:red, yerr=fill(Heights_err, length(X)), xerr=X_err, ms=2)
plot!(0.1:0.01:1.25, model(0.1:0.01:1.25, p), c=:blue, label="H* = $(round(H☆, digits=2)) ± $(round(H☆_err, digits=3)),\nβ=$(round(β, digits=2))  ± $(round(β_err, digits=3))", xlabel="X", ylabel="Heights (cm)", legend=:bottomleft, frame=:box)
savefig(FiguresDir * "Exp6-P1-H2.png")

############################ saving the data
df = DataFrame(FH1=PeakFreqs1, FH2=PeakFreqs2, Heights=Heights)
CSV.write(DataDir * "Exp6-P1.csv", df)
############################################## Exprience 6 Part 2 ##############################################
y, fs = load(DataDir * "Exp6-P2.wav")
fs = Int(fs)
y = mean(y, dims=2)
y = y .+ eps()

function GetHarmonicFreqs(i_T, fs, y; window=50000, TLen=10)
    PeakFreqs1 = []
    PeakFreqs2 = []
    Y_ = y[(i_T-1)*Int(TLen * fs)+1:i_T*Int(TLen * fs)]
    for i in 1:window:length(Y_)-window
        Y = Y_[i:i+window]
        fft_result = fft(Y)
        freqaxis = fftfreq(length(Y), fs)
        N_half = div(length(Y), 2)
        magnitude_db = 20 * log10.(abs.(fft_result))
        v_, i_ = findmax(magnitude_db[800:2000])
        push!(PeakFreqs1, freqaxis[i_+800])
        v_, i_ = findmax(magnitude_db[2500+150*i_T:3000+200*i_T])
        push!(PeakFreqs2, freqaxis[i_+2500+150*i_T])
    end
    return PeakFreqs1, PeakFreqs2
end

Heights = collect(1:10)
H1FE = 1505.0
H2FE = 3850.0
ErrFreqs1 = []
PeakFreqs1 = []
ErrFreqs2 = []
PeakFreqs2 = []
for i ∈ 1:10
    TPF1, TPF2 = GetHarmonicFreqs(i, fs, y)
    push!(PeakFreqs1, mean(TPF1))
    push!(PeakFreqs2, mean(TPF2))
    push!(ErrFreqs1, std(TPF1))
    push!(ErrFreqs2, std(TPF2))
end

############################ Plots:
############################ Frequencies
plot(Heights, PeakFreqs1, m=:c, label="First Harmonic Frequencies", c=:steelblue)
plot!(Heights, PeakFreqs2, m=:s, label="Second Harmonic Frequencies", c=:purple)
plot!(title="Frequency of Tones, None-Percussion", frame=:box, xlabel="Height (CM)", ylabel="Frequency (Hz)", legend=:right)
savefig(FiguresDir * "Exp6-P2-Freqs.png")

############################ First Harmonic
X = (((H1FE ./ PeakFreqs1) .^ 2) .- 1) .^ 0.25
Y = Heights
p0 = [12.0, 10.0]
@. model(x, p) = p[1] - (p[1] / (p[2]^0.25)) * x
fit = curve_fit(model, X, Y, p0)
fit_err = stderror(fit)
p = coef(fit)
H☆ = p[1]
β = p[2]
Heights_err = 0.2
H☆_err = fit_err[1]
β_err = fit_err[2]
X_err = X .* sqrt.((ErrFreqs1 ./ PeakFreqs1) .^ 2 .+ (fit_err[1] / p[1])^2)
scatter(X, Y, label="Data", c=:red, yerr=fill(Heights_err, length(X)), xerr=X_err, ms=2)
plot!(0.2:0.01:1.25, model(0.2:0.01:1.25, p), c=:blue, label="H* = $(round(H☆, digits=2)) ± $(round(H☆_err, digits=3)),\nβ=$(round(β, digits=2))  ± $(round(β_err, digits=3))", xlabel="X", ylabel="Heights (cm)", legend=:bottomleft, frame=:box)
savefig(FiguresDir * "Exp6-P2-H1.png")

############################ Second Harmonic
X = (((H2FE ./ PeakFreqs2) .^ 2) .- 1) .^ 0.25
Y = Heights
p0 = [12.0, 10.0]
@. model(x, p) = p[1] - (p[1] / (p[2]^0.25)) * x
fit = curve_fit(model, X, Y, p0)
fit_err = stderror(fit)
p = coef(fit)
H☆ = p[1]
β = p[2]
Heights_err = 0.2
H☆_err = fit_err[1]
β_err = fit_err[2]
X_err = X .* sqrt.((ErrFreqs2 ./ PeakFreqs2) .^ 2 .+ (fit_err[1] / p[1])^2)
scatter(X, Y, label="Data", c=:red, yerr=fill(Heights_err, length(X)), xerr=X_err, ms=2)
plot!(0.2:0.01:1.25, model(0.2:0.01:1.25, p), c=:blue, label="H* = $(round(H☆, digits=2)) ± $(round(H☆_err, digits=3)),\nβ=$(round(β, digits=2))  ± $(round(β_err, digits=3))", xlabel="X", ylabel="Heights (cm)", legend=:bottomleft, frame=:box)
savefig(FiguresDir * "Exp6-P2-H2.png")

############################ saving the data
df = DataFrame(FH1=PeakFreqs1, FH2=PeakFreqs2, FH1_err=ErrFreqs1, FH2_err=ErrFreqs2, Heights=Heights)
CSV.write(DataDir * "Exp6-P2.csv", df)