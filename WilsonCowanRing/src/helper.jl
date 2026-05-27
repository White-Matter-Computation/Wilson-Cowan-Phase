using DrWatson
using JLD2

function paramsdict(p)
    Dict(name => getfield(p, name) for name in fieldnames(typeof(p)))
end

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


function compute_or_load_hb(X0, T0, fp, wc)

    fname = datadir("hb_sol_" * savename(wc) * ".jld2")

    if isfile(fname)
        println("Loading existing HB solution: ", fname)

        data = load(fname)
        X_sol = data["X_sol"]
        T = data["T"]
        sol = data["sol"]

    else
        println("Computing HB solution...")

        X_sol, T, sol = compute_hb_orbit(X0, T0, fp, wc)

        @save fname X_sol T sol wc

    end

    return X_sol, T, sol
end

function compute_or_load_iprc(X_sol, T, fp, wc)

    fname = datadir("iprc_" * savename(wc) * ".jld2")

    if isfile(fname)
        println("Loading existing iPRC solution: ", fname)

        data = load(fname)
        A_iprc = data["A_iprc"]

    else
        println("Computing iPRC...")

        A_iprc = compute_iprc(X_sol, T, fp, wc)

        @save fname A_iprc

    end

    return A_iprc
end


function fourier_recon(t, A, T, p::FourierParams)
    dim = p.dim
    M = p.M
    A_reshaped = reshape(A, dim, :) #dim x N Fourier coeff
    result = zeros(ComplexF64, dim)

    for (idx, m) in enumerate(-M:M)
        for d in 1:dim
            result[d] += A_reshaped[d, idx] * exp(im * 2π * m * t / T)
        end
    end
    return real.(result)
end

function fourier_prime_recon(t, A, T, p::FourierParams)
    #d/dt
    dim = p.dim
    M = p.M
    A_reshaped = reshape(A, dim, :) #dim x N Fourier coeff
    result = zeros(ComplexF64, dim)

    for (idx, m) in enumerate(-M:M)
        for d in 1:dim
            result[d] += 2π * im * m / T * A_reshaped[d, idx] * exp(im * 2π * m * t / T)
        end
    end
    return real.(result)
end

using DSP
using Statistics
using LinearAlgebra
function comp(us, vs, ts, wc::WCParams; window_size=50)

    # --- Hilbert transform ---
    u_h = hilbert(us)
    v_h = hilbert(vs)

    # --- analytic phases ---
    ϕ_u = angle.(us .+ im .* imag.(u_h))
    ϕ_v = angle.(vs .+ im .* imag.(v_h))

    # --- Kuramoto order parameter over time ---
    R_u_t = abs.(mean(exp.(im .* ϕ_u), dims=1)) |> vec
    R_v_t = abs.(mean(exp.(im .* ϕ_v), dims=1)) |> vec

    # --- define final time window ---
    n = length(ts)
    idx = max(1, n - window_size):n

    R_u_final = mean(R_u_t[idx])
    R_v_final = mean(R_v_t[idx])

    R_u_std = std(R_u_t[idx])
    R_v_std = std(R_v_t[idx])

    #println("Final u synchrony: $R_u_final ± $R_u_std")
    #println("Final v synchrony: $R_v_final ± $R_v_std")

    N = length(us[:, end])
    u_end = us[:, end]
    v_end = vs[:, end]
    D_final = zeros(N, N)
    for i in 1:N
        for j in 1:N
            D_final[i, j] = sqrt((u_end[i] - u_end[j])^2 + (v_end[i] - v_end[j])^2)
        end
    end
    #println("Distance Matrix is with Determinant $(det(D_final)) and")# eigenvals $(eigvals(D_final))")
    up_diag = D_final[diagind(D_final, 1)]
    println("Upper Diagonal mean is $(mean(up_diag)) and sum is $(sum(up_diag))")
    #Center of Mass
    #println("CM is $(mean(u_end)),$(mean(v_end))")

    return [R_u_final, R_v_final], sum(up_diag)
end