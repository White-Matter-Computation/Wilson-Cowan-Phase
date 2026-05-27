"""
Arguments
ϵ::Float64
wc::WCParams
τ_1::Float64
n_neurons::Int
T_guess::Float64
wave_num::Int - corresponds to k in ψ=2πk/N

(Optional)
T_sample=2000 - How many samples to do 
tend_mult=100 - how many periods to simulate
X0_guess=[0.4, 0.6] - Guess for the orbit of the system with ϵ = 0 
abstol=1e-9
reltol=1e-8
plot_network=true - If true will do phase plot and a space time plot
init_pertub=nothing - can be used to peterb the initial condition for the network off the limit cycle
"""
function run_full_network(ϵ::Float64, wc::WCParams, τ_1::Float64, n_neurons::Int, T_guess::Float64, wave_num::Int; T_sample=2000, tend_mult=100, X0_guess=[0.4, 0.6], abstol=1e-9, reltol=1e-8, plot_network=true, init_pertub=nothing)
    #tend_mult is how many periods to simulate 

    alg = MethodOfSteps(Tsit5())

    ##########################
    ###  Find Limit Cycle  ###
    ##########################

    h_0(p, t) = X0_guess
    #For N=1 neruon
    prob_single = DDEProblem(wc_network, X0_guess, h_0, (0.0, 100 * T_guess), [ϵ, 1, wc, τ_1])
    sol_single = solve(prob_single, alg, abstol=abstol, reltol=reltol)
    x0_raw, T, last_peak = find_network_guess(sol_single, n_neurons)
    if isnothing(init_pertub)
        init_pertub = zeros(2 * n_neurons)
    else
        @assert length(init_pertub) == 2 * n_neurons
    end
    println("Found Limit cycle with T = $T, where the guess was $T_guess")
    ##########################
    ###  History Function  ###
    ##########################
    h = let sol = sol_single, T = T, last_peak = last_peak, wave_num = wave_num, init_pertub = init_pertub
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
            return x + init_pertub
        end
    end

    P = [ϵ, n_neurons, wc, τ_1]
    x0 = h(P, 0.0)
    tend = tend_mult * T
    tspan = (0.0, tend)

    prob = DDEProblem(wc_network, x0, h, tspan, P)
    network_sol = solve(prob, alg, saveat=tend / T_sample, abstol=1e-9, reltol=1e-8)

    X = Array(network_sol)
    t = network_sol.t
    u = X[1:2:end, :]
    v = X[2:2:end, :]
    if plot_network
        τ_k_max = τ_1 * floor(n_neurons / 2)
        t_hist = range(-τ_k_max, 0.0, length=200)
        H = hcat([h(P, t) for t in t_hist]...)
        u_hist = H[1:2:end, :]
        v_hist = H[2:2:end, :]
        u0 = X[1:2:end, 1]
        v0 = X[2:2:end, 1]
        un = X[1:2:end, end]
        vn = X[2:2:end, end]

        plot(u', v', xlabel="u", ylabel="v", title="τ_1 = $τ_1, τ_0 = $(wc.τ_0), ϵ=$(ϵ)", label=false)
        plot!(u_hist', v_hist', color=:black, alpha=0.5, linestyle=:dash, label=false)
        scatter!(u0, v0,
            marker_z=1:n_neurons,         # color by neuron index
            color=:viridis,
            colorbar=true,
            colorbar_title="neuron i",
            markershape=:circle,
            markersize=6,
            markerstrokewidth=0,
            label=false)
        scatter!(un, vn,
            marker_z=1:n_neurons,         # color by neuron index
            color=:viridis,
            markershape=:utriangle,
            markersize=6,
            markerstrokewidth=0,
            label=false)
        savefig("network_plot/phase_$(wave_num)_$(n_neurons)_$(τ_1)_$(wc.τ_0)_epsilon_$ϵ.svg")
    end

    return u, v, t
end
