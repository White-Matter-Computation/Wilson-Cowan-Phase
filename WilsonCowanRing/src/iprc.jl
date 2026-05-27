using BlockDiagonals
using Integrals

function compute_iprc(X_sol, T, p::FourierParams, wc::WCParams)
    dim = p.dim
    X_reshaped = reshape(X_sol, dim, :)
    S_p = p.S_p
    S_p_inv = p.S_p_inv
    A_phase = S_p_inv * X_sol
    N = p.N
    M = p.M
    t_range = [n * T / N for n in (-M:M)]

    J0_block = Matrix{Float64}[]

    #J0_block = Vector{Matrix{Float64}}(undef, N) #this should be faster
    for n in 1:N
        x_val = real(X_reshaped[1, n])
        y_val = real(X_reshaped[2, n])

        df_dx = -1 / wc.τ_e
        df_dy = 0

        dg_dx = 0
        dg_dy = -1 / wc.τ_i

        DF0 = [df_dx df_dy; dg_dx dg_dy]
        push!(J0_block, DF0)
    end

    J0 = BlockDiagonal(J0_block)


    J1_block = Matrix{Float64}[]

    for n in 1:N
        x_val = real(X_reshaped[1, n])
        #x_val2 = real(fourier_recon(t_range[n],S_p_inv*X_sol,T)[1])
        #println(abs(x_val-x_val2))
        y_val = real(X_reshaped[2, n])
        e_val, i_val = fourier_recon(t_range[n], S_p_inv * X_sol, T, p)
        e_d, i_d = fourier_recon(t_range[n] - wc.τ_0, S_p_inv * X_sol, T, p)

        dE = f_prime(wc.p_e + wc.w_ee * e_d - wc.w_ei * i_d, wc.β)
        dI = f_prime(wc.p_i + wc.w_ie * e_d - wc.w_ii * i_d, wc.β)

        df_dx = wc.w_ee * dE / wc.τ_e
        df_dy = -wc.w_ei * dE / wc.τ_e

        dg_dx = wc.w_ie * dI / wc.τ_i
        dg_dy = -wc.w_ii * dI / wc.τ_i

        DF1 = [df_dx df_dy; dg_dx dg_dy]
        push!(J1_block, DF1)
    end
    J1 = BlockDiagonal(J1_block)

    L_μ_p(μ) = kron(Diagonal([μ + im * 2π * m / T for m in -M:M]), I(dim))
    M_matrix(μ) = S_p * L_μ_p(μ) * S_p_inv - J0 - J1 * (S_p * Γ_p(T, p, wc.τ_0) * S_p_inv) * exp(-μ * wc.τ_0)



    M_0 = M_matrix(0.0)

    # Right eigenvector for μ=0 is tangent to limit cycle
    # Left eigenvector is the iPRC
    F = eigen(M_0')
    iprc_idx = argmin(abs.(F.values))  # Find eigenvalue closest to 0
    iprc = real.(F.vectors[:, iprc_idx])
    A_iprc = S_p_inv * iprc

    #Normalize 
    γ_0 = fourier_recon(0, A_phase, T, p) #γ(0)
    γ_τ = fourier_recon(0 - wc.τ_0, A_phase, T, p) #γ(-τ)
    q_τ = fourier_recon(0, A_iprc, T, p) # Q(0)
    #Normalize for <Q(0),F(γ(0),γ(-τ))> = 1
    dot0 = dot(q_τ, real(wc_rhs(wc, γ_0, γ_τ)))



    domain = (-wc.τ_0, 0)
    prob_int = IntegralProblem(int_prob, domain, [T, A_phase, A_iprc, p, wc])
    sol_int = solve(prob_int, QuadGKJL(); abstol=1e-9, reltol=1e-9)
    #sol_int2 = solve(prob_int, HCubatureJL(); abstol=1e-9, reltol=1e-9)


    iprc_reshape = reshape(iprc, dim, :)
    iprc_reshape .*= 2π / T / (real(sol_int[1]) + dot0)


    #A_iprc
    return S_p_inv * vec(iprc_reshape)
end


function int_prob(t, p)
    T, A_phase, A_iprc, fp, wc = p
    τ = wc.τ_0
    qx, qy = fourier_recon(t + τ, A_iprc, T, fp)

    #X = (e,i)
    X = fourier_recon(t, A_phase, T, fp)
    X_d = fourier_recon(t - τ, A_phase, T, fp)

    e_val, i_val = X
    e_d, i_d = X_d

    dE = f_prime(wc.p_e + wc.w_ee * e_d - wc.w_ei * i_d, wc.β)
    dI = f_prime(wc.p_i + wc.w_ie * e_d - wc.w_ii * i_d, wc.β)

    df_dx = wc.w_ee * dE / wc.τ_e
    df_dy = -wc.w_ei * dE / wc.τ_e

    dg_dx = wc.w_ie * dI / wc.τ_i
    dg_dy = -wc.w_ii * dI / wc.τ_i

    fx, gx = wc_rhs(wc, X, X_d)

    #Q' x DF1 x (f(X) g(X))
    return qx * (df_dx * fx + df_dy * gx) + qy * (dg_dx * fx + dg_dy * gx)
end
