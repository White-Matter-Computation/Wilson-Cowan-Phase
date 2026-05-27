using DifferentialEquations, DelayDiffEq
using Statistics

function compute_orbit(u0, history, tspan, wc::WCParams)
    prob = DDEProblem(wc_dde!, u0, history, tspan, wc; constant_lags=[wc.τ_0])
    sol_dde = solve(prob, MethodOfSteps(Tsit5()))

    return sol_dde
end

function find_network_guess(sol_dde, N; t_lim=50.0, th=0.98)

    t_end = sol_dde.t[end]
    t_lim_range = range(t_end - t_lim, t_end, length=1001)

    # use u variable to detect peaks
    x_lim = collect(sol_dde(t_lim_range, idxs=1))

    x_max = maximum(x_lim)
    x_th = th * x_max

    hit = false
    hit_points = Int[]

    for i in 1:length(x_lim)-1
        x = x_lim[i]

        if x > x_th && !hit && x > x_lim[i+1]
            push!(hit_points, i)
            hit = true
        elseif hit && x < x_th
            hit = false
        end
    end

    # estimate period
    t_guesses = [
        t_lim_range[hit_points[i]] - t_lim_range[hit_points[i-1]]
        for i in 2:length(hit_points)
    ]

    T = median(t_guesses)

    # last peak
    last_peak = t_lim_range[hit_points[end]]

    # sample along last cycle
    t_samples = range(last_peak - T, last_peak, length=N)

    u_samples = sol_dde(t_samples, idxs=1)
    v_samples = sol_dde(t_samples, idxs=2)

    # build interleaved state
    x0 = zeros(2N)

    for i in 1:N
        x0[2i-1] = u_samples[i]
        x0[2i] = v_samples[i]
    end

    return x0, T, last_peak
end

function find_init_guess(sol_dde, N, t_lim=50.0, th=0.98)
    #Find T
    N_total = 2 * N
    t_end = sol_dde.t[end]
    t_lim_range = range(t_end - t_lim, t_end, length=1001)
    #x_lim = [sol_dde(t)[1] for t in t_lim_range]
    x_lim = collect(sol_dde(t_lim_range, idxs=1))
    x_max = findmax(x_lim)[1]
    #Threshold to detect peak, note it must me lower than peak so we can check if function is increasing i.e. checking x[i] > x[i+1]
    x_th = th * x_max

    hit = false
    hit_points = []
    for (i, x) in enumerate(x_lim)
        if x > x_th && hit == false && x > x_lim[i+1]
            push!(hit_points, i)
            hit = true
        elseif hit == true && x < x_th
            hit = false
        end
    end
    #Get periods 
    t_guesses = [t_lim_range[hit_points[i]] - t_lim_range[hit_points[i-1]] for i in range(2, length(hit_points))]
    #use median guess but all should be the same
    if isempty(t_guesses)
        println("t_guess is Empty t_lim is from $(t_end-t_lim) to $(t_end) \n with hit points $(hit_points)")
    end
    T = median(t_guesses)

    last_peak = t_lim_range[hit_points[end]]
    t_samples = range(last_peak - T, last_peak, length=N)

    e_samples = sol_dde(t_samples, idxs=1)
    i_samples = sol_dde(t_samples, idxs=2)
    #Above is a lot better computationally 
    #e_samples = [sol_dde(t)[1] for t in t_samples]
    #i_samples = [sol_dde(t)[2] for t in t_samples]

    X0_physical = zeros(Float64, N_total)
    for i in 1:N
        X0_physical[2*(i-1)+1] = e_samples[i]
        X0_physical[2*(i-1)+2] = i_samples[i]
    end

    return X0_physical, T
end
