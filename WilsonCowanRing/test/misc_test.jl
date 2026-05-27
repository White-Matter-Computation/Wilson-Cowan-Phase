using DrWatson, Test
@quickactivate "WilsonCowanRing"
using TOML
using WilsonCowanRing
using JLD2
using DifferentialEquations
using Plots
using BenchmarkTools


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

M = 30
n_neurons = 11
wave_num = 0
T1_guess = 3.33
T2_guess = 10.0
#Pick different epsilon
ϵ1 = 0.1
ϵ2 = 0.01
#Pertubate initial condition
#pertub = fill(0.01, 2 * n_neurons) #a faster version of [0.01 for _ in 2*n_neurons]

t_end = 200
N = 2 * M + 1
N_total = 2 * N
fp = FourierParams(M, 2)
τ_0 = wc1.τ_0
u0 = [0.8, 0.8]
t_span = (0.0, t_end)

his(p, t) = u0
sol_dde1 = compute_orbit(u0, his, t_span, wc1)
X0_guess, T1_guess = find_init_guess(sol_dde1, N)
X_sol, T, sol = compute_or_load_hb(X0_guess, T1_guess, fp, wc1)
A_iprc = compute_or_load_iprc(X_sol, T, fp, wc1)
A_phase = fp.S_p_inv * X_sol
H_coeff = compute_H_coeff(wc1, fp, A_iprc, A_phase)

vals = range(0.0, 100.0; length=100)


function show_diff(x, y, z)
    if x ≈ y ≈ z
        #println("Equal")
    else
        println("Differences x,y,z are $(abs(x-y)) and $(abs(y-z))")
    end

end

# for x in vals
#     for y in vals
#         show_diff(compute_H_fourier_old(x, y, fp, H_coeff), compute_H_fourier(x, y, fp, H_coeff), compute_H_fourier2(x, y, fp, H_coeff))
#     end
# end

# x = 6.5
# y = 2.5
# @btime compute_H_fourier_old($x, $y, $fp, $H_coeff)
# @btime compute_H_fourier($x, $y, $fp, $H_coeff)
# @btime compute_H_fourier2($x, $y, $fp, $H_coeff)


#r1 = compute_G_coeff(wc1, fp, A_phase)
#r2 = compute_G_coeff2(wc1, fp, A_phase)

#println(findmax(abs.(r1 - r2)))

#@btime compute_G_coeff(wc1, fp, A_phase)
#@btime compute_G_coeff2(wc1, fp, A_phase)

#Weight display
#
idxs = collect(1:n_neurons)
res = zeros(n_neurons,n_neurons)
res = [w_k(i-j, n_neurons) for i in idxs, j in idxs]

heatmap(idxs, idxs, res'; color=:plasma,colorbar_title="w_k")
xlabel!("i")
ylabel!("j")
title!("Weight for a=2")
savefig("weight.svg")
#w_k(|i-j|,n_neurons)
#Test H_fourier
ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")
