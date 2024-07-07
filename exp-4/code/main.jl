using FileIO: load
using Plots, FFTW, StatsBase, LibSndFile
using CSV, DataFrames, LsqFit

DataDir = "../data/" # Directory of the data
FiguresDir = "../figure/" # Directory of the figures
cd(dirname(@__FILE__)) # Changing the current directory to the directory of the file

function GetMagnitudes(i_trial, fs, y; window=300)
    PeakMags = []
    times = []
    for i in 1:window:10*fs
        Y = y[i_trial+i:i_trial+window+i]
        fft_result = fft(Y)
        freqaxis = fftfreq(length(Y), fs)
        N_half = div(length(Y), 2)
        magnitude_db = 20 * log10.(abs.(fft_result[2:N_half]))
        v_, i_ = findmax(magnitude_db)
        push!(PeakMags, v_)
        push!(times, (i_trial + i) / fs)
    end
    return PeakMags, times
end

effective_area(T60, V) = 0.161 * V / T60

############################################## Exprience 4 Part 1 ##############################################
frequencies = [300, 450, 650, 900, 1300, 1800, 2500, 3500, 5000]
trialtimes = [10, 35, 60, 85, 110, 135, 160, 185, 210]
P1DataSets = [Dict("EffectiveArea" => [], "EAErr" => [], "T60" => [], "T60Err" => [], "Magnitude" => [], "Time" => []) for _ in 1:3]

for trial in 1:3
    y, fs = load(DataDir * "processed/Exp4-P1-T$trial.wav")
    fs = Int(fs)
    y = mean(y, dims=2)
    y = y .+ eps()
    for i in 1:9
        PeakMags, times = GetMagnitudes(trialtimes[i] * fs, fs, y, window=300)
        win = 10
        i_ = findfirst(x -> diff(PeakMags)[x] > 0 && mean(diff(PeakMags)[x:x+win]) > -1, 1:length(PeakMags)-win)
        Y = PeakMags[1:i_]
        X = times[1:i_]
        model(x, p) = p[1] * x .+ p[2]
        p0 = [1.0, 0.0]
        fit_ = curve_fit(model, X, Y, p0)
        p = coef(fit_)
        T60 = 60 / abs(p[1])
        T60Error = 60 * stderror(fit_)[1] / abs(p[1])^2
        V_room = 35
        EA = effective_area(T60, V_room)
        EAErr = EA * T60Error / T60
        push!(P1DataSets[trial]["Magnitude"], PeakMags)
        push!(P1DataSets[trial]["Time"], times)
        push!(P1DataSets[trial]["T60"], T60)
        push!(P1DataSets[trial]["T60Err"], T60Error)
        push!(P1DataSets[trial]["EffectiveArea"], EA)
        push!(P1DataSets[trial]["EAErr"], EAErr)
    end
end

############################################## Exprience 4 Part 2 ##############################################
frequencies = [300, 450, 650, 900, 1300, 1800, 2500, 3500, 5000]
trialtimes = [10, 35, 60, 85, 110, 135, 160, 185, 210]
P2DataSets = [Dict("EffectiveArea" => [], "EAErr" => [], "T60" => [], "T60Err" => [], "Magnitude" => [], "Time" => []) for _ in 1:3]

for trial in 1:3
    y, fs = load(DataDir * "processed/Exp4-P2-T$trial.wav")
    fs = Int(fs)
    y = mean(y, dims=2)
    y = y .+ eps()
    for i in 1:9
        PeakMags, times = GetMagnitudes(trialtimes[i] * fs, fs, y, window=300)
        MA_window = div(length(PeakMags), 20)
        MA_Mags = [mean(PeakMags[i:i+MA_window]) for i in 1:MA_window:length(PeakMags)-MA_window]
        MA_times = [mean(times[i:i+MA_window]) for i in 1:MA_window:length(PeakMags)-MA_window]
        i_ = findfirst(x -> diff(MA_Mags)[x] > 0 && mean(diff(MA_Mags)[x:x+3]) > -1, 1:length(MA_Mags))
        Y = MA_Mags[1:i_]
        X = MA_times[1:i_]
        model(x, p) = p[1] * x .+ p[2]
        p0 = [1.0, 0.0]
        fit_ = curve_fit(model, X, Y, p0)
        p = coef(fit_)
        T60 = 60 / abs(p[1])
        T60Error = 60 * stderror(fit_)[1] / abs(p[1])^2
        V_room = 60
        EA = effective_area(T60, V_room)
        EAErr = EA * T60Error / T60
        push!(P2DataSets[trial]["Magnitude"], MA_Mags)
        push!(P2DataSets[trial]["Time"], MA_times)
        push!(P2DataSets[trial]["T60"], T60)
        push!(P2DataSets[trial]["T60Err"], T60Error)
        push!(P2DataSets[trial]["EffectiveArea"], EA)
        push!(P2DataSets[trial]["EAErr"], EAErr)
    end
end

# Ploting the results

# Reverberation Time
Avgs1 = mean([P1DataSets[1]["T60"] P1DataSets[2]["T60"] P1DataSets[3]["T60"]], dims=2)
Err1 = mean([P1DataSets[1]["T60Err"] P1DataSets[2]["T60Err"] P1DataSets[3]["T60Err"]], dims=2)
Avgs2 = mean([P2DataSets[1]["T60"] P2DataSets[2]["T60"] P2DataSets[3]["T60"]], dims=2)
Err2 = mean([P2DataSets[1]["T60Err"] P2DataSets[2]["T60Err"] P2DataSets[3]["T60Err"]], dims=2)
plot(frequencies, Avgs1, c=:darkblue, label="Anechoic Chamber")
scatter!(frequencies, Avgs1, err=Err1, c=:blue, msc=:darkblue, label=nothing, ms=3)
plot!(frequencies, Avgs2, c=:darkred, label="Reverberation Chamber")
scatter!(frequencies, Avgs2, err=Err2, c=:red, msc=:darkred, label=nothing, ms=3)
plot!(title="Reverberation Time", xlabel="Frequency (Hz)", ylabel="T60 (s)",
    xticks=(frequencies, string.(frequencies)), xscale=:log2, yscale=:log10, legend=:left, size=(800, 600), frame=:box)
savefig(FiguresDir * "ReverberationTime.png")

# Effective Area
Avgs1 = mean([P1DataSets[1]["EffectiveArea"] P1DataSets[2]["EffectiveArea"] P1DataSets[3]["EffectiveArea"]], dims=2)
Err1 = mean([P1DataSets[1]["EAErr"] P1DataSets[2]["EAErr"] P1DataSets[3]["EAErr"]], dims=2)
Avgs2 = mean([P2DataSets[1]["EffectiveArea"] P2DataSets[2]["EffectiveArea"] P2DataSets[3]["EffectiveArea"]], dims=2)
Err2 = mean([P2DataSets[1]["EAErr"] P2DataSets[2]["EAErr"] P2DataSets[3]["EAErr"]], dims=2)
plot(frequencies, Avgs1, c=:darkblue, label="Anechoic Chamber")
scatter!(frequencies, Avgs1, err=Err1, c=:blue, ls=:dash, msc=:darkblue, label=nothing, ms=3)
plot!(frequencies, Avgs2, c=:darkred, label="Reverberation Chamber")
scatter!(frequencies, Avgs2, err=Err2, c=:red, ls=:dash, msc=:darkred, label=nothing, ms=3)
plot!(title="Effective Area", xlabel="Frequency (Hz)", ylabel="Area (m²)",
    xticks=(frequencies, string.(frequencies)), xscale=:log2, yscale=:log10, legend=:left, size=(800, 600), frame=:box)
savefig(FiguresDir * "EffectiveArea.png")

# Magnitude Time Series
T1Times = mean([P1DataSets[1]["Time"][1] P1DataSets[2]["Time"][1] P1DataSets[3]["Time"][1]], dims=2)
T1Mags = mean([P1DataSets[1]["Magnitude"][1] P1DataSets[2]["Magnitude"][1] P1DataSets[3]["Magnitude"][1]], dims=2)
T2Times = mean([P1DataSets[1]["Time"][2] P1DataSets[2]["Time"][2] P1DataSets[3]["Time"][2]], dims=2)
T2Mags = mean([P1DataSets[1]["Magnitude"][2] P1DataSets[2]["Magnitude"][2] P1DataSets[3]["Magnitude"][2]], dims=2)
T3Times = mean([P1DataSets[1]["Time"][3] P1DataSets[2]["Time"][3] P1DataSets[3]["Time"][3]], dims=2)
T3Mags = mean([P1DataSets[1]["Magnitude"][3] P1DataSets[2]["Magnitude"][3] P1DataSets[3]["Magnitude"][3]], dims=2)
T4Times = mean([P1DataSets[1]["Time"][4] P1DataSets[2]["Time"][4] P1DataSets[3]["Time"][4]], dims=2)
T4Mags = mean([P1DataSets[1]["Magnitude"][4] P1DataSets[2]["Magnitude"][4] P1DataSets[3]["Magnitude"][4]], dims=2)
T5Times = mean([P1DataSets[1]["Time"][5] P1DataSets[2]["Time"][5] P1DataSets[3]["Time"][5]], dims=2)
T5Mags = mean([P1DataSets[1]["Magnitude"][5] P1DataSets[2]["Magnitude"][5] P1DataSets[3]["Magnitude"][5]], dims=2)
T6Times = mean([P1DataSets[1]["Time"][6] P1DataSets[2]["Time"][6] P1DataSets[3]["Time"][6]], dims=2)
T6Mags = mean([P1DataSets[1]["Magnitude"][6] P1DataSets[2]["Magnitude"][6] P1DataSets[3]["Magnitude"][6]], dims=2)
T7Times = mean([P1DataSets[1]["Time"][7] P1DataSets[2]["Time"][7] P1DataSets[3]["Time"][7]], dims=2)
T7Mags = mean([P1DataSets[1]["Magnitude"][7] P1DataSets[2]["Magnitude"][7] P1DataSets[3]["Magnitude"][7]], dims=2)
T8Times = mean([P1DataSets[1]["Time"][8] P1DataSets[2]["Time"][8] P1DataSets[3]["Time"][8]], dims=2)
T8Mags = mean([P1DataSets[1]["Magnitude"][8] P1DataSets[2]["Magnitude"][8] P1DataSets[3]["Magnitude"][8]], dims=2)
T9Times = mean([P1DataSets[1]["Time"][9] P1DataSets[2]["Time"][9] P1DataSets[3]["Time"][9]], dims=2)
T9Mags = mean([P1DataSets[1]["Magnitude"][9] P1DataSets[2]["Magnitude"][9] P1DataSets[3]["Magnitude"][9]], dims=2)

plot(T1Times[1:35] .- trialtimes[1], T1Mags[1:35], label="f = $(frequencies[1])Hz")
plot!(T2Times[1:35] .- trialtimes[2], T2Mags[1:35], label="f = $(frequencies[2])Hz")
plot!(T3Times[1:35] .- trialtimes[3], T3Mags[1:35], label="f = $(frequencies[3])Hz")
plot!(T4Times[1:35] .- trialtimes[4], T4Mags[1:35], label="f = $(frequencies[4])Hz")
plot!(T5Times[1:35] .- trialtimes[5], T5Mags[1:35], label="f = $(frequencies[5])Hz")
plot!(T6Times[1:35] .- trialtimes[6], T6Mags[1:35], label="f = $(frequencies[6])Hz")
plot!(T7Times[1:35] .- trialtimes[7], T7Mags[1:35], label="f = $(frequencies[7])Hz")
plot!(T8Times[1:35] .- trialtimes[8], T8Mags[1:35], label="f = $(frequencies[8])Hz")
plot!(T9Times[1:35] .- trialtimes[9], T9Mags[1:35], label="f = $(frequencies[9])Hz")
plot!(title="Anechoic Chamber, Magnitude Time Series", xlabel="Time (s)", ylabel="Magnitude (dB)", size=(800, 600), frame=:box)
savefig(FiguresDir * "Exp4-P1-Magnitude-Time-Series.png")

T1Times = mean([P2DataSets[1]["Time"][1] P2DataSets[2]["Time"][1] P2DataSets[3]["Time"][1]], dims=2)
T1Mags = mean([P2DataSets[1]["Magnitude"][1] P2DataSets[2]["Magnitude"][1] P2DataSets[3]["Magnitude"][1]], dims=2)
T2Times = mean([P2DataSets[1]["Time"][2] P2DataSets[2]["Time"][2] P2DataSets[3]["Time"][2]], dims=2)
T2Mags = mean([P2DataSets[1]["Magnitude"][2] P2DataSets[2]["Magnitude"][2] P2DataSets[3]["Magnitude"][2]], dims=2)
T3Times = mean([P2DataSets[1]["Time"][3] P2DataSets[2]["Time"][3] P2DataSets[3]["Time"][3]], dims=2)
T3Mags = mean([P2DataSets[1]["Magnitude"][3] P2DataSets[2]["Magnitude"][3] P2DataSets[3]["Magnitude"][3]], dims=2)
T4Times = mean([P2DataSets[1]["Time"][4] P2DataSets[2]["Time"][4] P2DataSets[3]["Time"][4]], dims=2)
T4Mags = mean([P2DataSets[1]["Magnitude"][4] P2DataSets[2]["Magnitude"][4] P2DataSets[3]["Magnitude"][4]], dims=2)
T5Times = mean([P2DataSets[1]["Time"][5] P2DataSets[2]["Time"][5] P2DataSets[3]["Time"][5]], dims=2)
T5Mags = mean([P2DataSets[1]["Magnitude"][5] P2DataSets[2]["Magnitude"][5] P2DataSets[3]["Magnitude"][5]], dims=2)
T6Times = mean([P2DataSets[1]["Time"][6] P2DataSets[2]["Time"][6] P2DataSets[3]["Time"][6]], dims=2)
T6Mags = mean([P2DataSets[1]["Magnitude"][6] P2DataSets[2]["Magnitude"][6] P2DataSets[3]["Magnitude"][6]], dims=2)
T7Times = mean([P2DataSets[1]["Time"][7] P2DataSets[2]["Time"][7] P2DataSets[3]["Time"][7]], dims=2)
T7Mags = mean([P2DataSets[1]["Magnitude"][7] P2DataSets[2]["Magnitude"][7] P2DataSets[3]["Magnitude"][7]], dims=2)
T8Times = mean([P2DataSets[1]["Time"][8] P2DataSets[2]["Time"][8] P2DataSets[3]["Time"][8]], dims=2)
T8Mags = mean([P2DataSets[1]["Magnitude"][8] P2DataSets[2]["Magnitude"][8] P2DataSets[3]["Magnitude"][8]], dims=2)
T9Times = mean([P2DataSets[1]["Time"][9] P2DataSets[2]["Time"][9] P2DataSets[3]["Time"][9]], dims=2)
T9Mags = mean([P2DataSets[1]["Magnitude"][9] P2DataSets[2]["Magnitude"][9] P2DataSets[3]["Magnitude"][9]], dims=2)

plot(T1Times .- trialtimes[1], T1Mags, label="f = $(frequencies[1])Hz")
plot!(T2Times .- trialtimes[2], T2Mags, label="f = $(frequencies[2])Hz")
plot!(T3Times .- trialtimes[3], T3Mags, label="f = $(frequencies[3])Hz")
plot!(T4Times .- trialtimes[4], T4Mags, label="f = $(frequencies[4])Hz")
plot!(T5Times .- trialtimes[5], T5Mags, label="f = $(frequencies[5])Hz")
plot!(T6Times .- trialtimes[6], T6Mags, label="f = $(frequencies[6])Hz")
plot!(T7Times .- trialtimes[7], T7Mags, label="f = $(frequencies[7])Hz")
plot!(T8Times .- trialtimes[8], T8Mags, label="f = $(frequencies[8])Hz")
plot!(T9Times .- trialtimes[9], T9Mags, label="f = $(frequencies[9])Hz")
plot!(title="Reverberation Chamber, Magnitude Time Series", xlabel="Time (s)", ylabel="Magnitude (dB)", size=(800, 600), frame=:box)
savefig(FiguresDir * "Exp4-P2-Magnitude-Time-Series.png")

# Saving the results
# Part 1 Trial 1:
CSV.write(DataDir * "processed/Exp4-P1-T1.csv", DataFrame(Frequency=frequencies, T60=P1DataSets[1]["T60"], T60Err=P1DataSets[1]["T60Err"], EffectiveArea=P1DataSets[1]["EffectiveArea"], EAErr=P1DataSets[1]["EAErr"]))
# Part 1 Trial 2:
CSV.write(DataDir * "processed/Exp4-P1-T2.csv", DataFrame(Frequency=frequencies, T60=P1DataSets[2]["T60"], T60Err=P1DataSets[2]["T60Err"], EffectiveArea=P1DataSets[2]["EffectiveArea"], EAErr=P1DataSets[2]["EAErr"]))
# Part 1 Trial 3:
CSV.write(DataDir * "processed/Exp4-P1-T3.csv", DataFrame(Frequency=frequencies, T60=P1DataSets[3]["T60"], T60Err=P1DataSets[3]["T60Err"], EffectiveArea=P1DataSets[3]["EffectiveArea"], EAErr=P1DataSets[3]["EAErr"]))
# Part 2 Trial 1:
CSV.write(DataDir * "processed/Exp4-P2-T1.csv", DataFrame(Frequency=frequencies, T60=P2DataSets[1]["T60"], T60Err=P2DataSets[1]["T60Err"], EffectiveArea=P2DataSets[1]["EffectiveArea"], EAErr=P2DataSets[1]["EAErr"]))
# Part 2 Trial 2:
CSV.write(DataDir * "processed/Exp4-P2-T2.csv", DataFrame(Frequency=frequencies, T60=P2DataSets[2]["T60"], T60Err=P2DataSets[2]["T60Err"], EffectiveArea=P2DataSets[2]["EffectiveArea"], EAErr=P2DataSets[2]["EAErr"]))
# Part 2 Trial 3:
CSV.write(DataDir * "processed/Exp4-P2-T3.csv", DataFrame(Frequency=frequencies, T60=P2DataSets[3]["T60"], T60Err=P2DataSets[3]["T60Err"], EffectiveArea=P2DataSets[3]["EffectiveArea"], EAErr=P2DataSets[3]["EAErr"]))