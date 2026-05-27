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

#Expect is 
#τ_0=0.2 ψ=0.0 Stable τ_1∈ {2.5,3,6} and Unstable τ_1 ∈ {1,1.5,4}
#τ_0=1.5 ψ=0.0 Stable τ_1∈ {0.2,8,9} and Unstable τ_1 ∈ {2,5}




n_neurons = 11
wave_num = 0
T1_guess = 3.33
T2_guess = 10.0
#Pick different epsilon
ϵ1 = 0.001
#ϵ2 = 0.01
#Pertubate initial condition
#pertub = fill(0.01, 2 * n_neurons) #a faster version of [0.01 for _ in 2*n_neurons]
#pertub = collect(range(0, 0.02; length=2 * n_neurons))
min_perturb = 0.0
max_perturb = 0.03

pertub = min_perturb .+ (max_perturb - min_perturb) .* rand(2 * n_neurons)
#Solve the DDE network model
println("Using a perturbation \n $pertub")
u1_stab1, v1_stab1, t1_stab1 = run_full_network(ϵ1, wc1, 3.0, n_neurons, T1_guess, wave_num, init_pertub=pertub, T_sample=10000, tend_mult=300)
#u1_stab1, v1_stab1, t1_stab1 = run_full_network(ϵ1, wc1, 0.5, n_neurons, T1_guess, wave_num, init_pertub=pertub)
#u1_stab2, v1_stab2, t1_stab2 = run_full_network(ϵ1, wc1, 3.3, n_neurons, T1_guess, wave_num, init_pertub=pertub, T_sample=10000, tend_mult=300)
# u1_stab3, v1_stab3, t1_stab3 = run_full_network(ϵ1, wc1, 6.0, n_neurons, T1_guess, wave_num, init_pertub=pertub)

u1_unstab1, v1_unstab1, t1_unstab1 = run_full_network(ϵ1, wc1, 2.0, n_neurons, T1_guess, wave_num, init_pertub=pertub, T_sample=10000, tend_mult=300)
#u1_unstab1, v1_unstab1, t1_unstab1 = run_full_network(ϵ1, wc1, 1.5, n_neurons, T1_guess, wave_num, init_pertub=pertub)
#u1_unstab2, v1_unstab2, t1_unstab2 = run_full_network(ϵ1, wc1, 3.38, n_neurons, T1_guess, wave_num, init_pertub=pertub, T_sample=10000, tend_mult=300)
#u1_unstab3, v1_unstab3, t1_unstab3 = run_full_network(ϵ1, wc1, 3.4, n_neurons, T1_guess, wave_num, init_pertub=pertub, tend_mult=1000)

# u2_stab1, v2_stab1, t2_stab1 = run_full_network(ϵ2, wc2, 0.2, n_neurons, T2_guess, wave_num, init_pertub=pertub)
# u2_stab2, v2_stab2, t2_stab2 = run_full_network(ϵ2, wc2, 8.0, n_neurons, T2_guess, wave_num, init_pertub=pertub)
# u2_stab3, v2_stab3, t2_stab3 = run_full_network(ϵ2, wc2, 9.0, n_neurons, T2_guess, wave_num, init_pertub=pertub)

# u2_unstab1, v2_unstab1, t2_unstab1 = run_full_network(ϵ2, wc2, 2.0, n_neurons, T2_guess, wave_num, init_pertub=pertub)
# u2_unstab2, v2_unstab2, t2_unstab2 = run_full_network(ϵ2, wc2, 5.0, n_neurons, T2_guess, wave_num, init_pertub=pertub)



#plot_stefan2(u1_stab1, v1_stab1, t1_stab1, 2.91, wc1, n_neurons, wave_num, ϵ1, max_plot=50.0)
#plot_stefan2(u1_stab1, v1_stab1, t1_stab1, 3.3, wc1, n_neurons, wave_num, ϵ1, max_plot=50.0)
# plot_stefan2(u1_stab2, v1_stab2, t1_stab2, 3.0, wc1, n_neurons, wave_num, ϵ1)
# plot_stefan2(u1_stab3, v1_stab3, t1_stab3, 6.0, wc1, n_neurons, wave_num, ϵ1)

#plot_stefan2(u1_unstab1, v1_unstab1, t1_unstab1, 2.88, wc1, n_neurons, wave_num, ϵ1, max_plot=50.0)
#plot_stefan2(u1_unstab2, v1_unstab2, t1_unstab2, 3.38, wc1, n_neurons, wave_num, ϵ1, max_plot=50.0)
#plot_stefan2(u1_unstab2, v1_unstab2, t1_unstab2, 3.4, wc1, n_neurons, wave_num, ϵ1, max_plot=50.0)
# plot_stefan2(u1_unstab3, v1_unstab3, t1_unstab3, 4.0, wc1, n_neurons, wave_num, ϵ1)

comp(u1_stab1, v1_stab1, t1_stab1, wc1)
#comp(u1_stab2, v1_stab2, t1_stab2, wc1)
comp(u1_unstab1, v1_unstab1, t1_unstab1, wc1)
#comp(u1_unstab2, v1_unstab2, t1_unstab2, wc1)

plot_stefan2(u1_stab1, v1_stab1, t1_stab1, 3.0, wc1, n_neurons, wave_num, ϵ1,max_plot=50.0)
# plot_stefan2(u2_stab2, v2_stab2, t2_stab2, 8.0, wc2, n_neurons, wave_num, ϵ2)
# plot_stefan2(u2_stab3, v2_stab3, t2_stab3, 9.0, wc2, n_neurons, wave_num, ϵ2)

plot_stefan2(u1_unstab1, v1_unstab1, t1_unstab1, 2.0, wc1, n_neurons, wave_num, ϵ1,max_plot=50.0)
# plot_stefan2(u2_unstab2, v2_unstab2, t2_unstab2, 5.0, wc2, n_neurons, wave_num, ϵ2)


ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")
