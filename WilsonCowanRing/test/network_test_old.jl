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


τ_1 = 0.2#1.3  and 0.2

n_neurons = 11
wave_num = 1
T = 3.33
#Pick different epsilon
ϵ1 = 0.1
ϵ2 = 0.01
#Solve the DDE network model


alg = MethodOfSteps(Tsit5())
##########################
###  Find Limit Cycle  ###
#########################
X0_guess = [0.4, 0.6]
h_0(p, t) = X0_guess
prob_1 = DDEProblem(wc_network, X0_guess, h_0, (0.0, 100 * T), [ϵ1, 1, wc1, τ_1])
sol_1 = solve(prob_1, alg, abstol=1e-9, reltol=1e-8)
x0_1, T_1, last_peak_1 = find_network_guess(sol_1, n_neurons)


prob_2 = DDEProblem(wc_network, X0_guess, h_0, (0.0, 100 * T), [ϵ2, 1, wc2, τ_1])
sol_2 = solve(prob_2, alg, abstol=1e-9, reltol=1e-8)
x0_2, T_2, last_peak_2 = find_network_guess(sol_2, n_neurons)
####################
###### Plots #######
####################

h_1 = let sol = sol_1, T = T_1, last_peak = last_peak_1, wave_num = wave_num
    function (p, t)
        n_neurons = p[2]
        x = zeros(2 * n_neurons)
        for i in 1:n_neurons
            φ = wave_num * (i - 1) / n_neurons * T
            ti = last_peak - T + mod(t + φ, T)
            uv = sol(ti)
            x[2i-1] = uv[1]
            x[2i] = uv[2]
        end
        return x
    end
end

h_2 = let sol = sol_2, T = T_2, last_peak = last_peak_2, wave_num = wave_num
    function (p, t)
        n_neurons = p[2]
        x = zeros(2 * n_neurons)
        for i in 1:n_neurons
            φ = wave_num * (i - 1) / n_neurons * T
            ti = last_peak - T + mod(t + φ, T)
            uv = sol(ti)
            x[2i-1] = uv[1]
            x[2i] = uv[2]
        end
        return x
    end
end


P1 = [ϵ1, n_neurons, wc1, τ_1]
P2 = [ϵ2, n_neurons, wc2, τ_1]
#Using a guess from the h function is a lot better 
x0_1 = h_1(P1, 0.0)
x0_2 = h_2(P2, 0.0)

#Check if the init condition is on the cycle 
t_range = range(last_peak_1 - T_1, last_peak_1, length=200)
lc_full = hcat(sol_1.(t_range)...)
u_lc = lc_full[1, :]
v_lc = lc_full[2, :]

u_init = x0_1[1:2:end]
v_init = x0_1[2:2:end]
plot(u_lc, v_lc, label="limit cycle", legend=false)
#scatter!(u_init, v_init, color=:red, markersize=4, markerstrokewidth=0)
scatter!(u_init, v_init,
    marker_z=1:n_neurons,         # color by neuron index
    color=:viridis,
    colorbar=true,
    colorbar_title="neuron i",
    markersize=6,
    markerstrokewidth=0,
    label=false)
savefig("network_plot/limit_cycle_$(wave_num)_$(n_neurons)_$τ_1.svg")


tend = 100 * T_1
tspan = (0., tend)
# initalise history function



prob = DDEProblem(wc_network, x0_1, h_1, tspan, P1)
prob2 = DDEProblem(wc_network, x0_2, h_2, tspan, P2)
Tsample = 2000
network_sol = solve(prob, alg, saveat=tend / Tsample, abstol=1e-9, reltol=1e-8)
network_sol2 = solve(prob2, alg, saveat=tend / Tsample, abstol=1e-9, reltol=1e-8)

X = Array(network_sol)   # size = (2*n_neurons, Tsample+1)
t = network_sol.t
u = X[1:2:end, :]   # all u_i
v = X[2:2:end, :]   # all v_i

X2 = Array(network_sol2)   # size = (2*n_neurons, Tsample+1)
t2 = network_sol2.t
u2 = X2[1:2:end, :]   # all u_i
v2 = X2[2:2:end, :]   # all v_i

# plot(t, u', xlabel="t", ylabel="u_i", title="τ_1 = $τ_1, τ_0 = $(wc1.τ_0)", legend=false)
# savefig("network_plot/us_neurons_$(n_neurons)_$τ_1.svg")

# plot(t, v', xlabel="t", ylabel="v_i", title="τ_1 = $τ_1, τ_0 = $(wc1.τ_0)", legend=false)
# savefig("network_plot/vs_neurons_$(n_neurons)_$τ_1.svg")

# plot(t, u2', xlabel="t", ylabel="u_i", title="τ_1 = $τ_1, τ_0 = $(wc2.τ_0)", legend=false)
# savefig("network_plot/us_neurons_$(n_neurons)_wc2_$τ_1.svg")

# plot(t, v2', xlabel="t", ylabel="v_i", title="τ_1 = $τ_1, τ_0 = $(wc2.τ_0)", legend=false)
# savefig("network_plot/vs_neurons_$(n_neurons)_wc2_$τ_1.svg")


#Show history 
τ_k_max = τ_1 * floor(n_neurons / 2)
t_hist = range(-τ_k_max, 0.0, length=200)
H1 = hcat([h_1(P1, t) for t in t_hist]...)
u_hist1 = H1[1:2:end, :]
v_hist1 = H1[2:2:end, :]

u0_1 = X[1:2:end, 1]
v0_1 = X[2:2:end, 1]

un_1 = X[1:2:end, end]
vn_1 = X[2:2:end, end]

plot(u', v', xlabel="u", ylabel="v", title="τ_1 = $τ_1, τ_0 = $(wc1.τ_0), ϵ=$(ϵ1)", label=false)
plot!(u_hist1', v_hist1', color=:black, alpha=0.5, linestyle=:dash, label=false)
#scatter!(u0_1, v0_1, color=:red, markersize=4, markerstrokewidth=0)
#scatter!(un_1, vn_1, color=:green, markersize=4, markerstrokewidth=0)
scatter!(u0_1, v0_1,
    marker_z=1:n_neurons,         # color by neuron index
    color=:viridis,
    colorbar=true,
    colorbar_title="neuron i",
    markershape=:circle,
    markersize=6,
    markerstrokewidth=0,
    label=false)
scatter!(un_1, vn_1,
    marker_z=1:n_neurons,         # color by neuron index
    color=:viridis,
    markershape=:utriangle,
    markersize=6,
    markerstrokewidth=0,
    label=false)
savefig("network_plot/n1_neurons_$(wave_num)_$(n_neurons)_wc1_$(τ_1)_epsilon_$ϵ1.svg")


#Show history 
H2 = hcat([h_2(P2, t) for t in t_hist]...)
u_hist2 = H2[1:2:end, :]
v_hist2 = H2[2:2:end, :]

u0_2 = X2[1:2:end, 1]
v0_2 = X2[2:2:end, 1]

un_2 = X2[1:2:end, end]
vn_2 = X2[2:2:end, end]


plot(u2', v2', xlabel="u", ylabel="v", title="τ_1 = $τ_1, τ_0 = $(wc2.τ_0), ϵ=$(ϵ2)", label=false)
plot!(u_hist2', v_hist2', color=:black, alpha=0.5, linestyle=:dash, label=false)
#scatter!(u0_2, v0_2, color=:red, markersize=4, markerstrokewidth=0)
scatter!(u0_2, v0_2,
    marker_z=1:n_neurons,         # color by neuron index
    color=:viridis,
    colorbar=true,
    colorbar_title="neuron i",
    markershape=:circle,
    markersize=6,
    markerstrokewidth=0,
    label=false)
#scatter!(un_2, vn_2, color=:green, markersize=4, markerstrokewidth=0)
scatter!(un_2, vn_2,
    marker_z=1:n_neurons,         # color by neuron index
    color=:viridis,
    markershape=:utriangle,
    markersize=6,
    markerstrokewidth=0,
    label=false)

savefig("network_plot/n1_neurons$(wave_num)_$(n_neurons)_wc2_$(τ_1)_epsilon_$ϵ2.svg")
plot(u2[10, :], v2[10, :], xlabel="u", ylabel="v", title="u_10 vs v_10 for τ_1 = $τ_1, τ_0 = $(wc2.τ_0), ϵ=$(ϵ2)", label=false)
savefig("network_plot/$(wave_num)_u10_v10.svg")
plot_stefan(u, v, t, τ_1, wc1, n_neurons, wave_num, ϵ1)
plot_stefan(u2, v2, t2, τ_1, wc2, n_neurons, wave_num, ϵ2)

plot_stefan2(u, v, t, τ_1, wc1, n_neurons, wave_num, ϵ1)
plot_stefan2(u2, v2, t2, τ_1, wc2, n_neurons, wave_num, ϵ2)

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")