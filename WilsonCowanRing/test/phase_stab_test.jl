using DrWatson, Test
@quickactivate "WilsonCowanRing"
ENV["GKSwstype"] = "nul"
using Plots
gr()
using TOML
using WilsonCowanRing
using JLD2

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

params1 = all_params["set1"]
params2 = all_params["set2"]
wc1 = WCParams(; convert_keys(all_params["set1"])...)
wc2 = WCParams(; convert_keys(all_params["set2"])...)
#mby round those so your save data is actually nice 
#τ0_range = range(0.1, 2.0; length=50)
τ0_range = 0.1:0.1:2.0 #Instead of range so we can have nice τ_0 number (more usefull when we save and load data too)

#stability_mult_tau0(params1, τ0_range)

n_neurons = 11

ψ_0 = 0.0
ψ_1 = 2π / n_neurons
ψ_2 = 4π / n_neurons
ψ_3 = 6π / n_neurons
ψ_4 = 8π / n_neurons
ψ_5 = 10π / n_neurons

ϵ1 = 0.1

stability_single_tau0(wc1; n_neurons=n_neurons, ψ=ψ_0, ϵ=ϵ1, τ1_max=20.0, filename="dde_batch_$(wc1.τ_0)_$(round(ψ_0;digits=2))_eps_$(ϵ1)_20.mat", N_max=20000)
# stability_single_tau0(wc1; n_neurons=n_neurons, ψ=ψ_1, ϵ=ϵ1, τ1_max=10.0)
# stability_single_tau0(wc1; n_neurons=n_neurons, ψ=ψ_2, ϵ=ϵ1, τ1_max=10.0)
# stability_single_tau0(wc1; n_neurons=n_neurons, ψ=ψ_3, ϵ=ϵ1, τ1_max=10.0)
# stability_single_tau0(wc1; n_neurons=n_neurons, ψ=ψ_4, ϵ=ϵ1, τ1_max=10.0)
# stability_single_tau0(wc1; n_neurons=n_neurons, ψ=ψ_5, ϵ=ϵ1, τ1_max=10.0)

# stability_single_tau0(wc2; n_neurons=n_neurons, ψ=ψ_0)
# stability_single_tau0(wc2; n_neurons=n_neurons, ψ=ψ_1)
# stability_single_tau0(wc2; n_neurons=n_neurons, ψ=ψ_2)
# stability_single_tau0(wc2; n_neurons=n_neurons, ψ=ψ_3)
# stability_single_tau0(wc2; n_neurons=n_neurons, ψ=ψ_4)
# stability_single_tau0(wc2; n_neurons=n_neurons, ψ=ψ_5)
