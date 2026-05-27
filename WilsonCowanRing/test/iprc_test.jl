#using Revise
using DrWatson, Test
@quickactivate "WilsonCowanRing"
ENV["GKSwstype"] = "nul"
using Plots
gr()
using TOML
using WilsonCowanRing
using JLD2
using DifferentialEquations



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

e0 = 0.8
i0 = 0.8
t_end = 100.0
#modes 
M = 30
#points
N = 2 * M + 1
#total points
N_total = 2 * N


u0 = [e0, i0]
t_span = (0.0, t_end)

his(p, t) = u0

sol_dde1 = compute_orbit(u0, his, t_span, wc1)
sol_dde2 = compute_orbit(u0, his, t_span, wc2)

plot_orbit(sol_dde1, wc1)
plot_orbit(sol_dde2, wc2)

X0_guess1, T1_guess = find_init_guess(sol_dde1, N)
X0_guess2, T2_guess = find_init_guess(sol_dde2, N)

# plot_init_guess(X0_guess1, wc1)
# plot_init_guess(X0_guess2, wc2)

#fourier params 
fp = FourierParams(M, 2)

X1_sol, T1, sol1 = compute_or_load_hb(X0_guess1, T1_guess, fp, wc1)
X2_sol, T2, sol2 = compute_or_load_hb(X0_guess2, T2_guess, fp, wc2)

A_iprc1 = compute_or_load_iprc(X1_sol, T1, fp, wc1)
A_iprc2 = compute_or_load_iprc(X2_sol, T2, fp, wc2)

A_phase1 = fp.S_p_inv * X1_sol
A_phase2 = fp.S_p_inv * X2_sol
# plot_iprc(A_iprc1, T1, fp, wc1)
# plot_iprc(A_iprc2, T2, fp, wc2)

H_coeff1 = compute_H_coeff(wc1, fp, A_iprc1, A_phase1)
H_coeff2 = compute_H_coeff(wc2, fp, A_iprc2, A_phase2)


###################################################################
###### Below compares fourier and integral computation of H #######
###################################################################
#Compare H computation
# Xs = range(0, 2π; length=100)
# Ys = range(0, 2π; length=100)

# Hx_fourier = zeros(length(Xs), length(Ys))
# Hx_integral = zeros(length(Xs), length(Ys))
# H_zeros = zeros(length(Ys))


# for (i, X) in enumerate(Xs)
#     for (j, Y) in enumerate(Ys)
#         Hx_fourier[i, j] = real(compute_H_fourier(X, Y, fp, H_coeff2)[1])
#         Hx_integral[i, j] = real(compute_H_integral(X, Y, wc2, fp, A_iprc2, A_phase2))
#         H_zeros[j] = real(compute_H_fourier(0, Y, fp, H_coeff1))
#     end
# end

# #Uncomment to compare to integral solution

# surface(Xs, Ys, Hx_fourier'; alpha=0.6, label="Fourier")
# surface!(Xs, Ys, Hx_integral'; alpha=0.6, label="Integral")

# savefig("H_comparison_3d_wc2.svg")

# Hx_diff = abs.(Hx_integral .- Hx_fourier)
# heatmap(Xs, Ys, Hx_diff'; alpha=0.6, color=:plasma, label="Integral")

# savefig("H_compare_heatmap_wc2.svg")

# plot(
#     heatmap(Xs, Ys, Hx_fourier'; color=:plasma, title="Fourier"),
#     heatmap(Xs, Ys, Hx_integral'; color=:plasma, title="Integral"),
#     layout=(1, 2)
# )
# savefig("H_compare_heatmap_side_wc2.svg")

n_neurons = 11
ψ = 2π / n_neurons#0.0
ω1 = 2π / T1
ω2 = 2π / T2
ϵ1 = 0.1
ϵ2 = 0.01
params_Ω1 = ParamsΩ(ϵ1, ω1, wc1.τ_0, ψ, n_neurons, fp, H_coeff1)
params_Ω2 = ParamsΩ(ϵ2, ω2, wc2.τ_0, ψ, n_neurons, fp, H_coeff2)

#contour_phase_lock(params_Ω1, wc1, v_min=1.0, v_max=3.0, w_min=0.0, w_max=5.0)
#contour_phase_lock(params_Ω2, wc2, v_min=0.0, v_max=5.0, w_min=0.0, w_max=5.0)

Ωs_adap, τs_adap, X0, V = continuation_phase_lock(params_Ω1, wc1)
plot_continuation(Ωs_adap, τs_adap, X0, V, params_Ω1, wc1)
plot_cont_cont(params_Ω1, wc1, Ωs_adap, τs_adap, v_min=1.0, v_max=3.0, w_min=0.0, w_max=5.0)
Ωs_adap2, τs_adap2, X02, V2 = continuation_phase_lock(params_Ω2, wc2, t_0=0.1, Ω_guess=0.7)
plot_continuation(Ωs_adap2, τs_adap2, X02, V2, params_Ω2, wc2)
plot_cont_cont(params_Ω2, wc2, Ωs_adap2, τs_adap2, v_min=0.0, v_max=1.0, w_min=0.0, w_max=15.0)




# diag_A, diag_C, B = continuation_stability(0.8, 0.2, 0.01, n_neurons, ψ, wc1, fp, H_coeff1)
# save_for_ddebiftool(diag_A, diag_C, B, wc1.τ_0, "dde_coeffs.mat")

#batch_continuation_stability(Ωs_adap2[1:2:end], τs_adap2[1:2:end], ϵ2, n_neurons, ψ, wc2, fp, H_coeff2)
#batch_continuation_stability(Ωs_adap[1:2:end], τs_adap[1:2:end], ϵ1, n_neurons, ψ, wc1, fp, H_coeff1, filename="dde_batch2.mat")
#diagnose_single_point(Ωs_adap2[2], τs_adap2[2], ϵ2, n_neurons, ψ, wc2, fp, H_coeff2)

#Forward Evolve Phase reduced Equation 
p_forward1 = [ϵ1, ω1, n_neurons, wc1, 0.2, fp, H_coeff1]
p_forward2 = [ϵ2, ω2, n_neurons, wc2, 0.2, fp, H_coeff2]

wave_num = 1

# h_0 = let sol = sol1, T = T1, wave_num = wave_num
#     function (p, t)
#         n_neurons = p[3]
#         x = zeros(n_neurons)
#         for i in 1:n_neurons
#             φ = wave_num * (i - 1) / n_neurons * T
#             ti = sol1.t[end] - T + mod(t + φ, T)
#             x[i] = sol(ti)
#         end
#         return x
#     end
# end
θ0_guess = zeros(n_neurons)
function h_0(p, t)
    return θ0_guess
end
prob_red1 = DDEProblem(wc_reduced, θ0_guess, h_0, (0.0, 300), p_forward1)
sol_red1 = solve(prob_red1, MethodOfSteps(Tsit5()), abstol=1e-9, reltol=1e-8)
plot(sol_red1')
savefig("asdas.svg")
ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")
