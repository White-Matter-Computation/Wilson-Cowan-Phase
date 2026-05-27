using NonlinearSolve, LinearAlgebra

struct FourierParams
    M::Int
    dim::Int
    N::Int
    k::Vector{Int}
    S_p::Matrix{ComplexF64}
    S_p_inv::Matrix{ComplexF64}
end
struct HBParams
    fourier::FourierParams
    wc::WCParams
end


function FourierParams(M, dim)

    k = collect(-M:M)
    N = length(k)

    S = [exp(2π * im * n * m / N) for n in k, m in k]
    S_inv = conj(S) / N

    S_p = kron(S, I(dim))
    S_p_inv = kron(S_inv, I(dim))

    return FourierParams(M, dim, N, k, S_p, S_p_inv)
end

function phase_condition(X, p::FourierParams)

    Y = p.S_p_inv * X
    A = reshape(Y, p.dim, :)

    A_m = A[1, :]

    return real(dot(p.k .* im, A_m))
end

function Γ_p(T, p::FourierParams, τ)

    Γ = Diagonal(exp.(-im * τ .* 2π .* p.k ./ T))
    return kron(Γ, I(p.dim))
end

function L_p(T, p::FourierParams)

    L = Diagonal(im * 2π .* p.k ./ T)
    return kron(L, I(p.dim))
end

function residual!(F, X, params::HBParams)
    p = params.fourier
    wc = params.wc

    N_xs = length(X) - 1
    X_vars = @view X[1:N_xs]
    T = X[end]

    Γ = Γ_p(T, p, wc.τ_0)
    L = L_p(T, p)

    delay_state = p.S_p * Γ * p.S_p_inv * X_vars

    main_eq =
        p.S_p * L * p.S_p_inv * X_vars -
        dde_rhs_fourier(X_vars, delay_state, wc, p.M)

    F[1:N_xs] .= real.(main_eq)

    phase_eq = phase_condition(X_vars, p)

    F[end] = 1e4 * phase_eq
end


function compute_hb_orbit(X0, T0, p::FourierParams, wc::WCParams)
    params = HBParams(p, wc)
    X0_extended = [X0; T0]
    prob = NonlinearProblem(residual!, X0_extended, params)

    sol = solve(prob, NewtonRaphson())

    X_sol_extended = sol.u

    N_total = length(X_sol_extended) - 1
    X_sol = X_sol_extended[1:N_total]
    T = X_sol_extended[end]

    return X_sol, T, sol
end