using Plots, Statistics
FigureDirectory = "../figure/"
cd(dirname(@__FILE__)) # Changing the current directory to the directory of the file

Ferequency = [125, 250, 500, 1000, 2000, 4000, 8000, 16000]
Y_L_Increasing = [[-2, -3, -2], [-6, -8, -10], [-8, -8, -10], [-22, -24, -23], [-27, -29, -28], [-28, -29, -28], [-33, -33, -34], [-2, -3, -3]]
Y_R_Increasing = [[-1, -2, -3], [-5, -7, -9], [-8, -8, -9], [-21, -25, -23], [-26, -28, -28], [-27, -28, -28], [-34, -32, -33], [-1, -2, -2]]
Y_L_Decreasing = [[-4, -3, -6], [-9, -10, -11], [-15, -16, -18], [-29, -31, -31], [-29, -30, -31], [-27, -28, -28], [-33, -33, -34], [0, -3, -5]]
Y_R_Decreasing = [[-4, -4, -5], [-9, -11, -11], [-14, -16, -19], [-29, -30, -31], [-28, -29, -30], [-28, -27, -29], [-33, -34, -35], [-1, -4, -5]]
F_L_Increasing = [[18, 20, 18], [15, 14, 14], [7, 8, 7], [-5, -5, -5], [4, 3, 3], [-22, -19, -20], [11, 10, 10], [36, 36, 36]]
F_R_Increasing = [[18, 21, 18], [15, 13, 15], [6, 9, 8], [-4, -4, -6], [5, 4, 4], [-21, -20, -19], [12, 11, 9], [36, 36, 36]]
F_L_Decreasing = [[19, 20, 19], [15, 13, 15], [8, 7, 6], [-4, -6, -7], [5, 4, 4], [-22, -21, -18], [11, 10, 11], [36, 36, 36]]
F_R_Decreasing = [[19, 19, 19], [16, 14, 14], [8, 8, 7], [-4, -5, -6], [6, 5, 5], [-21, -20, -17], [11, 9, 12], [36, 36, 36]]

Y_L_Increasing_mean = mean.(Y_L_Increasing)
Y_R_Increasing_mean = mean.(Y_R_Increasing)
Y_L_Decreasing_mean = mean.(Y_L_Decreasing)
Y_R_Decreasing_mean = mean.(Y_R_Decreasing)
F_L_Increasing_mean = mean.(F_L_Increasing)
F_R_Increasing_mean = mean.(F_R_Increasing)
F_L_Decreasing_mean = mean.(F_L_Decreasing)
F_R_Decreasing_mean = mean.(F_R_Decreasing)

Y_L_Increasing_std = std.(Y_L_Increasing)
Y_R_Increasing_std = std.(Y_R_Increasing)
Y_L_Decreasing_std = std.(Y_L_Decreasing)
Y_R_Decreasing_std = std.(Y_R_Decreasing)
F_L_Increasing_std = std.(F_L_Increasing)
F_R_Increasing_std = std.(F_R_Increasing)
F_L_Decreasing_std = std.(F_L_Decreasing)
F_R_Decreasing_std = std.(F_R_Decreasing)

Y_Plot = begin
    plot(Ferequency, Y_L_Increasing_mean, err=Y_L_Increasing_std, label="Increasing, Left Ear",
        title="Yaghoub's Data", xlabel="Frequency (Hz)", ylabel="Gain (dB)", xscale=:log10, xticks=(Ferequency, string.(Ferequency)))
    plot!(Ferequency, Y_L_Decreasing_mean, err=Y_L_Decreasing_std, label="Decreasing, Right Ear")
    plot!(Ferequency, Y_R_Increasing_mean, err=Y_R_Increasing_std, label="Increasing, Right Ear")
    plot!(Ferequency, Y_R_Decreasing_mean, err=Y_R_Decreasing_std, label="Decreasing, Right Ear", legend=:top, frame=:box)
end
savefig(Y_Plot, FigureDirectory * "Yaghoub-EB.png")

F_Plot = begin
    plot(Ferequency, F_L_Increasing_mean, err=F_L_Increasing_std, label="Increasing, Left Ear",
        title="Fereshteh's Data", xlabel="Frequency (Hz)", ylabel="Gain (dB)", xscale=:log10, xticks=(Ferequency, string.(Ferequency)))
    plot!(Ferequency, F_L_Decreasing_mean, err=F_L_Decreasing_std, label="Decreasing, Left Ear")
    plot!(Ferequency, F_R_Increasing_mean, err=F_R_Increasing_std, label="Increasing, Right Ear")
    plot!(Ferequency, F_R_Decreasing_mean, err=F_R_Decreasing_std, label="Decreasing, Right Ear", legend=:top, frame=:box)
end
savefig(F_Plot, FigureDirectory * "Fereshteh-EB.png")

Y_Plot = begin
    plot(Ferequency, Y_L_Increasing_mean, ribbon=Y_L_Increasing_std, label="Increasing, Left Ear",
        title="Yaghoub's Data", xlabel="Frequency (Hz)", ylabel="Gain (dB)", xscale=:log10, xticks=(Ferequency, string.(Ferequency)))
    plot!(Ferequency, Y_L_Decreasing_mean, ribbon=Y_L_Decreasing_std, label="Decreasing, Right Ear")
    plot!(Ferequency, Y_R_Increasing_mean, ribbon=Y_R_Increasing_std, label="Increasing, Right Ear")
    plot!(Ferequency, Y_R_Decreasing_mean, ribbon=Y_R_Decreasing_std, label="Decreasing, Right Ear", legend=:top, frame=:box)
end
savefig(Y_Plot, FigureDirectory * "Yaghoub-R.png")

F_Plot = begin
    plot(Ferequency, F_L_Increasing_mean, ribbon=F_L_Increasing_std, label="Increasing, Left Ear",
        title="Fereshteh's Data", xlabel="Frequency (Hz)", ylabel="Gain (dB)", xscale=:log10, xticks=(Ferequency, string.(Ferequency)))
    plot!(Ferequency, F_L_Decreasing_mean, ribbon=F_L_Decreasing_std, label="Decreasing, Left Ear")
    plot!(Ferequency, F_R_Increasing_mean, ribbon=F_R_Increasing_std, label="Increasing, Right Ear")
    plot!(Ferequency, F_R_Decreasing_mean, ribbon=F_R_Decreasing_std, label="Decreasing, Right Ear", legend=:top, frame=:box)
end
savefig(F_Plot, FigureDirectory * "Fereshteh-R.png")
