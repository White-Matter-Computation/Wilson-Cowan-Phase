using FastGaussQuadrature
using Roots

struct ParamsΩ
    ϵ::Float64
    ω::Float64
    τ_0::Float64
    ψ::Float64
    n_neurons::Int
    fp::FourierParams
    H_coeff::Matrix{ComplexF64}
end

function G(X, Y, wc::WCParams)
    gx = f_prime(wc.p_e + wc.w_ee * X[1] - wc.w_ei * X[2], wc.β) * Y[1]

    return [gx, 0]
end

function compute_G_coeff_old(wc::WCParams, fp::FourierParams, A_phase::Vector{ComplexF64})

    fname = datadir("G_coeff_" * savename(wc) * ".jld2")

    if isfile(fname)
        @load fname G_coeff
        return G_coeff
    end
    M = fp.M
    N = fp.N
    G_coeff = zeros(ComplexF64, 2, N, N)
    n_l = 100
    x_l, w_l = gausslegendre(n_l)

    #map [-1,1] to [0,2π]
    u_l = (x_l .+ 1) .* π
    w_l = w_l .* π
    println("\nStarting to compute G_coeff with $n_l Gauss-Legendre nodes\n")
    t1 = time()
    F_vals = [fourier_recon(u, A_phase, 2π, fp) for u in u_l]
    # Precompute exponentials
    exp_n = [exp(-im * n * u) for n in -M:M, u in u_l]
    exp_m = [exp(-im * m * u) for m in -M:M, u in u_l]

    for n in -M:M
        n_idx = n + M + 1

        for m in -M:M
            m_idx = m + M + 1

            acc = zeros(ComplexF64, 2)
            for i in 1:n_l
                for j in 1:n_l
                    gval = G(F_vals[i], F_vals[j], wc)
                    weight = w_l[i] * w_l[j]

                    acc .+= weight * gval * exp_n[n_idx, i] * exp_m[m_idx, j]
                end
            end

            G_coeff[:, n_idx, m_idx] .= acc / (2π)^2
        end
    end
    t2 = time()
    println("\nI have finished computing G_coeff\n")
    println("It took $(t2-t1) seconds to finish or $((t2-t1)/60) minutes")


    @save fname G_coeff

    return G_coeff
end


function compute_G_coeff(wc::WCParams, fp::FourierParams, A_phase::Vector{ComplexF64})

    fname = datadir("G_coeff_" * savename(wc) * ".jld2")

    if isfile(fname)
        @load fname G_coeff
        return G_coeff
    end
    M = fp.M
    N = fp.N
    G_coeff = zeros(ComplexF64, 2, N, N)
    n_l = 100
    x_l, w_l = gausslegendre(n_l)

    #map [-1,1] to [0,2π]
    u_l = (x_l .+ 1) .* π
    w_l = w_l .* π
    println("\nStarting to compute G_coeff with $n_l Gauss-Legendre nodes\n")

    t1 = time()

    F_vals = [fourier_recon(u, A_phase, 2π, fp) for u in u_l]
    # Precompute weighted G matrix: A_c[i,j] = w_i * w_j * G_c[i,j]
    # G returns a length-2 vector so G_mat is (2, n_l, n_l)
    G_mat = zeros(ComplexF64, 2, n_l, n_l)

    #! You have to launch Julia with more cores by doing e.g.
    #julia --threads 4
    Threads.@threads for i in 1:n_l
        for j in 1:n_l
            G_mat[:, i, j] = w_l[i] * w_l[j] * G(F_vals[i], F_vals[j], wc)
        end
    end

    exp_n = [exp(-im * n * u) for n in -M:M, u in u_l]
    exp_m = [exp(-im * m * u) for m in -M:M, u in u_l]

    # For each component: G_coeff[c,:,:] = exp_n * G_mat[c,:,:] * exp_m^T / (2π)^2
    G_coeff = zeros(ComplexF64, 2, N, N)
    for c in 1:2
        G_coeff[c, :, :] = (exp_n * G_mat[c, :, :] * transpose(exp_m)) / (2π)^2
    end
    t2 = time()
    println("\nI have finished computing G_coeff\n")
    println("It took $(t2-t1) seconds to finish or $((t2-t1)/60) minutes")


    @save fname G_coeff

    return G_coeff
end


#So G_nm = G_coeff[:,n,m] where n,m ∈ (1,2,...,N)
#H_nm = Z_{-n-m} * G_nm if -n-m ∈ (-M,...,M) otherwise 0 

function compute_H_coeff(wc::WCParams, fp::FourierParams, A_iprc::AbstractVector{ComplexF64}, A_phase::AbstractVector{ComplexF64})

    fname = datadir("H_coeff_" * savename(wc) * ".jld2")

    if isfile(fname)
        @load fname H_coeff
        return H_coeff
    end

    N = fp.N
    M = fp.M
    G_coeff = compute_G_coeff(wc, fp, A_phase)
    H_coeff = zeros(ComplexF64, N, N)
    Q_coeff = reshape(A_iprc, fp.dim, :) #Q_coeff[:,i] = Q_i = Z_i

    @inbounds for n in -M:M
        n_idx = n + M + 1
        for m in -M:M
            m_idx = m + M + 1
            #Compute -n-m in terms of (1,...,N) index
            j_idx = -n - m + M + 1

            #below is the same as j >= -M and j<= M
            if j_idx >= 1 && j_idx <= N
                #Careful here the multiplication is dot product and ' is Transpose + CONJUGATE in Julia
                #dot is weird and does congugate automatically
                H_coeff[n_idx, m_idx] = dot(Q_coeff[:, j_idx]', G_coeff[:, n_idx, m_idx])
            end
        end
    end

    @save fname H_coeff

    return H_coeff

end

#Old 
function compute_H_fourier_old(X, Y, fp::FourierParams, H_coeff)
    #result = 0.0
    result = ComplexF64(0.0)
    M = fp.M
    for n in -M:M
        n_idx = n + M + 1
        for m in -M:M
            m_idx = m + M + 1
            result += H_coeff[n_idx, m_idx] * exp(im * (n * X + m * Y))
        end
    end

    return real(result)
end

# function compute_H_fourier(X, Y, fp::FourierParams, H_coeff)
#     M = fp.M

#     expX = [exp(im * n * X) for n in -M:M]
#     expY = [exp(im * m * Y) for m in -M:M]

#     result = 0.0 + 0.0im

#     for n in 1:(2M+1)
#         inner = 0.0 + 0.0im
#         for m in 1:(2M+1)
#             inner += H_coeff[n, m] * expY[m]
#         end
#         result += inner * expX[n]
#     end

#     return real(result)
# end
#The BLAS way
function compute_H_fourier(X, Y, fp::FourierParams, H_coeff)
    M = fp.M

    expX = [exp(im * n * X) for n in -M:M]
    expY = [exp(im * m * Y) for m in -M:M]

    result = transpose(expX) * H_coeff * expY

    return real(result)
end

"""Just for comparison, when doing computation use compute_H_fourier"""
function H_integral(t, X, Y, wc::WCParams, fp::FourierParams, A_iprc, A_phase)
    #I think this should be 2π periodic as we integrate over the phase 
    q_t = fourier_recon(t, A_iprc, 2π, fp) #Q(t) or Z(t)

    γ_tx = fourier_recon(t + X, A_phase, 2π, fp) # γ(t+X)
    γ_ty = fourier_recon(t + Y, A_phase, 2π, fp) # γ(t+Y)


    g_t = G(γ_tx, γ_ty, wc)

    #Notice that gy = 0 thus H_y = 0, so we can focus only on the X 

    return dot(q_t', g_t) / 2π
end

function compute_H_integral(X, Y, wc::WCParams, fp::FourierParams, A_iprc, A_phase)
    domain_H = (0, 2π)
    prob_int_H = IntegralProblem((t, p) -> H_integral(t, X, Y, wc, fp, A_iprc, A_phase), domain_H)
    sol_int_H = solve(prob_int_H, QuadGKJL(); abstol=1e-9, reltol=1e-9)

    return sol_int_H[1]
end

function F1(Ω::Float64, t::Float64, params_Ω::ParamsΩ)
    (; ϵ, ω, τ_0, ψ, n_neurons, fp, H_coeff) = params_Ω


    s = 0.0
    for k in 0:(n_neurons-1)
        s += w_k(k, n_neurons) * real(compute_H_fourier(-Ω * τ_0, k * ψ - τ_k(k, n_neurons, t) * Ω, fp, H_coeff))
    end
    return ω - Ω + ϵ * s
end


function compute_H_x_fourier(X, Y, fp::FourierParams, H_coeff::Matrix{ComplexF64})
    #result = 0.0
    result = ComplexF64(0.0)
    M = fp.M
    for n in -M:M
        n_idx = n + M + 1
        for m in -M:M
            m_idx = m + M + 1
            result += H_coeff[n_idx, m_idx] * im * n * exp(im * (n * X + m * Y))
        end
    end

    return real.(result)
end


function compute_H_y_fourier(X, Y, fp::FourierParams, H_coeff::Matrix{ComplexF64})
    #result = 0.0
    result = ComplexF64(0.0)
    M = fp.M
    for n in -M:M
        n_idx = n + M + 1
        for m in -M:M
            m_idx = m + M + 1
            result += H_coeff[n_idx, m_idx] * im * m * exp(im * (n * X + m * Y))
        end
    end

    return real.(result)
end

function dF1dΩ(Ω::Float64, t::Float64, params_Ω::ParamsΩ)
    (; ϵ, ω, τ_0, ψ, n_neurons, fp, H_coeff) = params_Ω

    s = 0.0
    for k in 0:(n_neurons-1)
        τk = τ_k(k, n_neurons, t)
        s += w_k(k, n_neurons) * (-τ_0 * compute_H_x_fourier(-Ω * τ_0, k * ψ - τk * Ω, fp, H_coeff) - τk * compute_H_y_fourier(-Ω * τ_0, k * ψ - τk * Ω, fp, H_coeff))
    end
    return real(ϵ * s) - 1
end

function dF1dt(Ω::Float64, t::Float64, params_Ω::ParamsΩ)
    (; ϵ, ω, τ_0, ψ, n_neurons, fp, H_coeff) = params_Ω

    s = 0.0
    for k in 0:(n_neurons-1)
        τk = τ_k(k, n_neurons, t)
        s += w_k(k, n_neurons) * -min(abs(k), n_neurons - abs(k)) * Ω * compute_H_y_fourier(-Ω * τ_0, k * ψ - τk * Ω, fp, H_coeff)
    end
    return real(ϵ * s)
end


function contour_phase_lock(params_Ω::ParamsΩ, wc::WCParams; v_min=2.0, v_max=4.5, w_min=-5.0, w_max=5.0)
    vs = v_min:0.01:v_max
    ws = w_min:0.01:w_max
    (; ϵ, ω, τ_0, ψ, n_neurons, fp, H_coeff) = params_Ω
    fname = datadir("contour_omega$(ω)_psi_$(ψ)_n_neurons_$(n_neurons)_eps_$(ϵ)_vmin_$(v_min)_v_max_$(v_max)_w_min_$(w_min)_w_max_$(w_max)" * savename(wc) * ".jld2")

    if isfile(fname)
        @load fname Z
    else
        Z = [F1(x, p, params_Ω) for x in vs, p in ws]
        @save fname Z
    end


    pname = plotsdir("contour_omega$(ω)_psi_$(ψ)_n_neurons_$n_neurons" * savename(wc) * ".svg")
    contour(ws, vs, Z; levels=[0], colorbar=false, label="Contour plot", c=:red)
    xlabel!("τ1")
    ylabel!("Ω")
    title!("Contour plot for ψ = $(round(ψ,digits=3)), N = $n_neurons with τ0 = $(wc.τ_0), ϵ=$(ϵ)")
    savefig(pname)
end

function find_initial_guess_continuation(params_Ω::ParamsΩ;
    v_min=0.01, v_max=4.5,
    w_min=0.01, w_max=0.5,
    v_res=1000, w_res=100
)
    (; ϵ, ω, τ_0, ψ, n_neurons, fp, H_coeff) = params_Ω
    # vs = range(v_min, v_max, length=v_res)   # Ω range
    # ws = range(w_min, w_max, length=w_res)   # τ_1 range
    vs = v_min:0.01:v_max
    ws = w_min:0.01:w_max

    # Evaluate F1 on coarse grid
    Z = [F1(v, w, params_Ω) for v in vs, w in ws]

    # Find sign changes along each column (fixed τ_1, vary Ω)
    candidates = Tuple{Float64,Float64}[]
    for (j, w) in enumerate(ws)
        for i in 1:length(vs)-1
            if Z[i, j] * Z[i+1, j] < 0      # sign change between row i and i+1
                # Linear interpolation to refine Ω crossing
                Ω_interp = vs[i] - Z[i, j] * (vs[i+1] - vs[i]) / (Z[i+1, j] - Z[i, j])
                push!(candidates, (Ω_interp, Float64(w)))
            end
        end
    end

    if isempty(candidates)
        println("!!Plotting !!")
        pname = plotsdir("error_contour_omega$(ω)_psi_$(ψ)_n_neurons_$n_neurons.svg")
        contour(ws, vs, Z; levels=[0], colorbar=false, label="Contour plot", c=:red)
        xlabel!("τ1")
        ylabel!("Ω")
        title!("Contour plot for ψ = $(round(ψ,digits=3)), N = $n_neurons ϵ=$(ϵ)")
        savefig(pname)
    end
    isempty(candidates) && error("No zero crossing found in grid — expand v/w range")

    # Pick the candidate closest to the middle of the τ_1 range

    τ_mid = (w_min + w_max) / 2
    best = argmin(abs(c[2] - τ_mid) for c in candidates)
    Ω_0, τ_0_guess = candidates[best]

    println("  Found $(length(candidates)) candidates, using Ω_0=$(round(Ω_0,digits=4)), τ_1=$(round(τ_0_guess,digits=4))")
    return Ω_0, τ_0_guess
end

#function continuation_phase_lock(params_Ω::ParamsΩ, wc::WCParams; h_min=0.00001, h_max=0.01, h_con=0.001, N_max=5000, t_0=nothing, Ω_guess=nothing)
function continuation_phase_lock(params_Ω::ParamsΩ, wc::WCParams;
    h_min=0.00001, h_max=0.01, h_con=0.001, N_max::Int=5000,
    t_0=nothing, Ω_guess=nothing,      # now optional
    v_min=0.01, v_max=4.5,
    w_min=0.01, w_max=0.5, p_min=-Inf, p_max=Inf
)
    # Auto-detect guess if not provided
    if isnothing(t_0) || isnothing(Ω_guess)
        Ω_guess_auto, t_0_auto = find_initial_guess_continuation(params_Ω;
            v_min=v_min, v_max=v_max,
            w_min=w_min, w_max=w_max
        )
        t_0 = isnothing(t_0) ? t_0_auto : t_0
        Ω_guess = isnothing(Ω_guess) ? Ω_guess_auto : Ω_guess
    end

    f_root(x) = F1(x, t_0, params_Ω)
    Ω_0 = find_zero(f_root, Ω_guess)

    v_x = -dF1dt(Ω_0, t_0, params_Ω)
    v_y = dF1dΩ(Ω_0, t_0, params_Ω)
    V = [v_x; v_y]
    V /= -norm(V)
    Ωs_adap, τs_adap = cont(F1, Ω_0, t_0, h_con, N_max, V, params=params_Ω, dFdx=dF1dΩ, dFdp=dF1dt, adaptive_step=true, h_min=h_min, h_max=h_max, newton_maxitr=10, output_diverge=true, p_min=p_min, p_max=p_max)
    #Ωs_adap, τs_adap = cont(F1, Ω_0, t_0, h_con, N_max, V, params=params_Ω, dFdx=dF1dΩ, dFdp=dF1dt, adaptive_step=false, h_min=h_min, h_max=h_max, newton_maxitr=10, output_diverge=true)

    return Ωs_adap, τs_adap, [Ω_0, t_0], V
end


