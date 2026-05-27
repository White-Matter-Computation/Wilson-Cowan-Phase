
using DrWatson, Test
@quickactivate "WilsonCowanRing"
using TOML
using WilsonCowanRing
using JLD2
using DifferentialEquations
using Plots
function convert_keys(d)
    Dict(
        :w_ee => d["w_ee"],
        :w_ei => d["w_ei"],
        :w_ie => d["w_ie"],
        :w_ii => d["w_ii"],
        :p_e => d["p_e"],
        :p_i => d["p_i"],
        :τ_e => d["tau_e"],
        :τ_i => d["tau_i"],
        :τ_0 => d["tau_0"],
        :β => d["beta"]
    )
end

println("Starting test")
ti = time()
all_params = TOML.parsefile("params.toml")
wc1 = WCParams(; convert_keys(all_params["set1"])...)
wc2 = WCParams(; convert_keys(all_params["set2"])...)





n_neurons = 11
wave_num = 0
T1_guess = 3.33
T2_guess = 10.0
#Pick different epsilon
ϵ = 0.001
c = 3
#Pertubate initial condition
#Use new one 
min_perturb = 0.0
max_perturb = 0.008
tend_mult = c / ϵ
pertub = min_perturb .+ (max_perturb - min_perturb) .* rand(2 * n_neurons)
#Solve the DDE network model

#For a single branch get the first and last stable points 
#! For ϵ = 0.001 we have branches in (2.9005,3.3247) and in (6.1960,6.6365)
#! For ϵ = 0.01 we have branches in (2.8342,3.3216) and in (6.0611,6.6141)
# stab_min = 2.8342
# stab_max = 3.3216
stab_min = 2.9005
stab_max = 3.3247
#How far from those values to look
stab_pertubs = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
#All values for which we expect it to be stable 
stable_range = vcat(stab_min .+ stab_pertubs, stab_max .- stab_pertubs)
unstable_range = vcat(stab_min .- stab_pertubs, stab_max .+ stab_pertubs)

all_τ1s = vcat(stable_range, unstable_range)
expected_stab = vcat(trues(length(stable_range)), falses(length(unstable_range)))


n_tests = length(all_τ1s)
observed_s = zeros(Float64, n_tests)

for (idx, τ_1) in enumerate(all_τ1s)
    u, v, t = run_full_network(ϵ, wc1, τ_1, n_neurons, T1_guess, wave_num, init_pertub=pertub, tend_mult=tend_mult)
    _, s = comp(u, v, t, wc1)     # discard R with _
    observed_s[idx] = s
    println("τ_1=$(round(τ_1,digits=4))  expected=$(expected_stab[idx] ? "stable" : "unstable")  s=$(round(s,digits=4))")
end

syn_threshold = max_perturb
observed_stab = observed_s .< syn_threshold    # s < 0.01 → stable

correct = observed_stab .== expected_stab
false_pos = .!expected_stab .&& observed_stab
false_neg = expected_stab .&& .!observed_stab

println("\nAccuracy: $(sum(correct)) / $n_tests ($(round(100*sum(correct)/n_tests, digits=1))%)")
println("False positives: $(sum(false_pos))")
println("False negatives: $(sum(false_neg))")

# Plot s vs τ_1
# p = scatter(title="Stability validation",
#     xlabel="τ₁", ylabel="s",
#     legend=:outertopright,
#     ylims=(-0.005, max(0.05, maximum(observed_s)) * 1.1),
#     size=(1200, 700), dpi=150,
#     margin=5Plots.mm,
#     tickfontsize=11, labelfontsize=13, legendfontsize=11)

# scatter!(p, all_τ1s[correct.&&expected_stab], observed_s[correct.&&expected_stab],
#     label="Correct stable", marker=:circle, color=:blue, ms=8)
# scatter!(p, all_τ1s[correct.&&.!expected_stab], observed_s[correct.&&.!expected_stab],
#     label="Correct unstable", marker=:square, color=:red, ms=8)
# scatter!(p, all_τ1s[false_neg], observed_s[false_neg],
#     label="False neg", marker=:xcross, color=:orange, ms=10, msw=3)
# scatter!(p, all_τ1s[false_pos], observed_s[false_pos],
#     label="False pos", marker=:xcross, color=:purple, ms=10, msw=3)

# hline!(p, [syn_threshold], color=:black, linestyle=:dash, label="Threshold (s=0.01)")
# vline!(p, [stab_min, stab_max], color=:grey, linestyle=:dot, label="Predicted boundary")
zoom_pad = 0.05   # how much to show either side of the boundary

p_main = scatter(title="For ϵ=$ϵ and threshold = $syn_threshold and max pertubation = $max_perturb after $tend_mult periods pass",
    xlabel="τ₁", ylabel="s",
    legend=:outertopright,
    size=(1200, 400), dpi=150, margin=5Plots.mm,
    tickfontsize=11, labelfontsize=13, legendfontsize=11)

p_left = scatter(title="Zoom: stab\\_min ($(round(stab_min, digits=3)))",
    xlabel="τ₁", ylabel="s",
    xlims=(stab_min - zoom_pad, stab_min + zoom_pad),
    legend=false, dpi=150, margin=5Plots.mm,
    tickfontsize=11, labelfontsize=13)

p_right = scatter(title="Zoom: stab\\_max ($(round(stab_max, digits=3)))",
    xlabel="τ₁", ylabel="s",
    xlims=(stab_max - zoom_pad, stab_max + zoom_pad),
    legend=false, dpi=150, margin=5Plots.mm,
    tickfontsize=11, labelfontsize=13)

# Plot onto all three panels
for p in (p_main, p_left, p_right)
    scatter!(p, all_τ1s[correct.&&expected_stab], observed_s[correct.&&expected_stab],
        label="Correct stable", marker=:circle, color=:blue, ms=3)
    scatter!(p, all_τ1s[correct.&&.!expected_stab], observed_s[correct.&&.!expected_stab],
        label="Correct unstable", marker=:square, color=:red, ms=3)
    scatter!(p, all_τ1s[false_neg], observed_s[false_neg],
        label="False neg", marker=:xcross, color=:orange, ms=4, msw=3)
    scatter!(p, all_τ1s[false_pos], observed_s[false_pos],
        label="False pos", marker=:xcross, color=:purple, ms=4, msw=3)
    hline!(p, [syn_threshold], color=:black, linestyle=:dash, label="Threshold")
    vline!(p, [stab_min, stab_max], color=:grey, linestyle=:dot, label="Boundary")
end

# Combine into one figure: main on top, two zooms below
p_combined = plot(p_main, p_left, p_right,
    layout=@layout([a{0.6h}; b c]),   # main takes top 60%, zooms split bottom
    size=(1200, 800),
    dpi=150)

display(p_combined)
savefig(p_combined, "compare_$(ϵ)_$max_perturb.svg")