#Does not include the funcions as it will almost always be a sigmoid
Base.@kwdef struct WCParams
    w_ee::Float64
    w_ei::Float64
    w_ie::Float64
    w_ii::Float64
    p_e::Float64
    p_i::Float64
    τ_e::Float64
    τ_i::Float64
    τ_0::Float64
    β::Float64
end

function f(x, β)
    return 1.0 / (1.0 + exp(-β * x))
end

function f_prime(x, β)
    return β * exp(-β * x) / (1 + exp(-β * x))^2
end

function wc_rhs(wc::WCParams, u_t, u_delay)

    e_t, i_t = real.(u_t)
    e_d, i_d = real.(u_delay)

    du = similar(u_t)

    du[1] = (-e_t + f(wc.p_e + wc.w_ee * e_d - wc.w_ei * i_d, wc.β)) / wc.τ_e
    du[2] = (-i_t + f(wc.p_i + wc.w_ie * e_d - wc.w_ii * i_d, wc.β)) / wc.τ_i

    return du
end


function wc_dde!(du, u, h, wc::WCParams, t)

    e, i = u
    u_delay = h(wc, t - wc.τ_0)
    e_hist, i_hist = u_delay

    du[1] = (-e + f(wc.p_e + wc.w_ee * e_hist - wc.w_ei * i_hist, wc.β)) / wc.τ_e
    du[2] = (-i + f(wc.p_i + wc.w_ie * e_hist - wc.w_ii * i_hist, wc.β)) / wc.τ_i

    return nothing
end

function dde_rhs_fourier(X, X_delay, wc::WCParams, M::Int64)

    X_t = reshape(X, 2, :)
    X_d = reshape(X_delay, 2, :)

    dX = similar(X_t)
    N = 2 * M + 1
    for n in 1:N
        u_t = @view X_t[:, n]
        u_d = @view X_d[:, n]

        dX[:, n] .= wc_rhs(wc, u_t, u_d)
    end
    return vec(dX)
end

function w_k(k, n_neurons)
    if k == 0
        return 0.0
    end
    a = 2
    return exp(-min(abs(k), n_neurons - abs(k)) / a)
end

function τ_k(k, n_neurons, τ_1)
    return τ_1 * min(abs(k), n_neurons - abs(k))
end

function wc_network(du, u, h, p, t)

    ϵ, n_neurons, wc, τ_1 = p
    histτ0 = h(p, t - wc.τ_0)

    for i in 1:n_neurons

        ui = u[2i-1]
        vi = u[2i]

        uτ0 = histτ0[2i-1]
        vτ0 = histτ0[2i]

        coupling = 0.0
        for j in 1:n_neurons
            uτij = h(p, t - τ_k(i - j, n_neurons, τ_1))[2j-1]
            coupling += w_k(i - j, n_neurons) * uτij
        end

        du[2i-1] =
            (-ui +
             f(wc.p_e + wc.w_ee * uτ0 - wc.w_ei * vτ0 + ϵ * coupling, wc.β)) / wc.τ_e

        du[2i] =
            (-vi + f(wc.p_i + wc.w_ie * uτ0 - wc.w_ii * vτ0, wc.β)) / wc.τ_i
    end
end


function wc_reduced(du, u, h, p, t)
    ϵ, ω, n_neurons, wc, τ_1, fp, H_coeff = p
    histτ0 = h(p, t - wc.τ_0)
    for i in 1:n_neurons
        θi = u[i]
        θi_τ0 = histτ0[i]
        coupling = 0.0
        for j in 1:n_neurons
            θj_τij = h(p, t - τ_k(i - j, n_neurons, τ_1))[j]

            Δ1 = mod(θi_τ0 - θi, 2π)
            Δ2 = mod(θj_τij - θi, 2π)
            coupling += w_k(i - j, n_neurons) * compute_H_fourier(Δ1, Δ2, fp, H_coeff)
        end
        du[i] = ω + ϵ * coupling
    end
end