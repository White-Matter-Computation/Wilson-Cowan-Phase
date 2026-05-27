using MAT

function continuation_stability(Ω::Float64, τ_1::Float64, ε::Float64,
    N::Int, ψ::Float64, wc::WCParams,
    fp::FourierParams, H_coeff::Matrix{ComplexF64})
    τ_0 = wc.τ_0

    Hx_map = Dict{Int,Float64}()
    Hy_map = Dict{Int,Float64}()
    w_map = Dict{Int,Float64}()
    τ_map = Dict{Int,Float64}()

    @inbounds for k in 0:(N-1)
        τd = round(τ_k(k, N, τ_1), digits=10)
        Hx_map[k] = compute_H_x_fourier(-Ω * τ_0, ψ * k - Ω * τd, fp, H_coeff)
        Hy_map[k] = compute_H_y_fourier(-Ω * τ_0, ψ * k - Ω * τd, fp, H_coeff)
        w_map[k] = w_k(k, N)
        τ_map[k] = τd
    end

    # Diagonals of A and C
    diag_A = zeros(N)
    diag_C = zeros(N)
    @inbounds for i in 1:N
        sx, sy = 0.0, 0.0
        #!can you skip loop
        for j in 1:N
            d = abs(i - j)
            sx += w_map[d] * Hx_map[d]
            sy += w_map[d] * Hy_map[d]
        end
        diag_C[i] = ε * sx
        diag_A[i] = -ε * (sx + sy)
    end

    # B matrices keyed by distinct delay value
    τ_vals = unique(values(τ_map))
    B = Dict(θ => zeros(N, N) for θ in τ_vals)
    @inbounds for i in 1:N, j in 1:N
        d = abs(i - j)#min(abs(i - j), N - abs(i - j))
        #! Possible issue with += instead of =
        B[τ_map[d]][i, j] += ε * w_map[d] * Hy_map[d]
    end

    return diag_A, diag_C, B
end



function save_for_ddebiftool(
    diag_A::Vector{Float64},
    diag_C::Vector{Float64},
    B::Dict{Float64,Matrix{Float64}},
    τ_0::Float64,
    filename::String
)
    N = length(diag_A)

    # Start with diagonal A and C as full matrices
    mat_A = Diagonal(diag_A) |> Matrix{Float64}
    mat_C = Diagonal(diag_C) |> Matrix{Float64}

    # Absorb any B[θ] whose delay coincides with 0 or τ_0 into A and C
    remaining = Dict{Float64,Matrix{Float64}}()
    for (θ, Bmat) in B
        if θ ≈ 0.0
            mat_A .+= Bmat          # x(t) contribution → merge into A
        elseif θ ≈ τ_0
            mat_C .+= Bmat          # x(t-τ_0) contribution → merge into C
        else
            remaining[θ] = Bmat
        end
    end

    # Sort remaining distinct delays
    θ_sorted = sort(collect(keys(remaining)))
    n_delays = 2 + length(θ_sorted)

    Tau = vcat(0.0, τ_0, θ_sorted)
    A_3d = zeros(Float64, N, N, n_delays)
    A_3d[:, :, 1] = mat_A
    A_3d[:, :, 2] = mat_C
    for (idx, θ) in enumerate(θ_sorted)
        A_3d[:, :, 2+idx] = remaining[θ]
    end

    matwrite(filename, Dict(
        "A_coeff" => A_3d,
        "Tau" => Tau,
        "N" => Float64(N),
        "tau_0" => τ_0,
    ))
    println("Saved: $(filename)")
    println("  N         = $N")
    println("  Tau       = $Tau")
    println("  A_coeff size = $(size(A_3d))")
end


function batch_continuation_stability(
    omegas::Vector{Float64},
    tau1s::Vector{Float64},
    ε::Float64,
    N::Int,
    ψ::Float64,
    wc::WCParams,
    fp::FourierParams,
    H_coeff::Matrix{ComplexF64};
    filename::String="dde_batch.mat"
)
    @assert length(omegas) == length(tau1s) "omegas and tau1s must have the same length"
    n_runs = length(omegas)

    A_cells = Vector{Array{Float64,3}}(undef, n_runs)
    Tau_cells = Vector{Vector{Float64}}(undef, n_runs)

    Threads.@threads for i in 1:n_runs
        diag_A, diag_C, B = continuation_stability(
            omegas[i], tau1s[i], ε, N, ψ, wc, fp, H_coeff
        )
        A_3d, Tau = _build_dde_arrays(diag_A, diag_C, B, wc.τ_0)
        A_cells[i] = A_3d
        Tau_cells[i] = Tau
    end

    A_cell_mat = Array{Any}(undef, 1, n_runs)   # 1×n_runs cell in MATLAB
    Tau_cell_mat = Array{Any}(undef, 1, n_runs)
    for i in 1:n_runs
        A_cell_mat[i] = A_cells[i]
        Tau_cell_mat[i] = Tau_cells[i]'            # row vector in MATLAB
    end

    filename = joinpath("matlab_data", filename)

    matwrite(filename, Dict(
        "A_coeff" => A_cell_mat,
        "Tau" => Tau_cell_mat,
        "omegas" => omegas',                       # save as row vectors
        "tau1s" => tau1s',
        "N" => Float64(N),
        "tau_0" => wc.τ_0,
    ))

    println("Saved $n_runs runs → $filename")
    #println("  Delay counts: ", [size(A_cells[i], 3) for i in 1:n_runs])
end

function _build_dde_arrays(
    diag_A::Vector{Float64},
    diag_C::Vector{Float64},
    B::Dict{Float64,Matrix{Float64}},
    τ_0::Float64;
    atol=1e-6
)
    N = length(diag_A)
    mat_A = Matrix(Diagonal(diag_A))
    mat_C = Matrix(Diagonal(diag_C))

    remaining = Dict{Float64,Matrix{Float64}}()
    #basically removes x(t) and x(t-τ0) from B
    for (θ, Bmat) in B
        if abs(θ) < atol
            mat_A .+= Bmat
        elseif abs(θ - τ_0) < atol
            mat_C .+= Bmat
        else
            remaining[θ] = Bmat
        end
    end

    θ_sorted = sort(collect(keys(remaining)))
    Tau = vcat(0.0, τ_0, θ_sorted)
    A_3d = zeros(Float64, N, N, 2 + length(θ_sorted))
    A_3d[:, :, 1] = mat_A
    A_3d[:, :, 2] = mat_C
    for (idx, θ) in enumerate(θ_sorted)
        A_3d[:, :, 2+idx] = remaining[θ]
    end

    return A_3d, Tau
end

function diagnose_single_point(Ω::Float64, τ_1::Float64, ε::Float64, N::Int, ψ::Float64, wc::WCParams, fp::FourierParams, H_coeff)
    diag_A, diag_C, B = continuation_stability(Ω, τ_1, ε, N, ψ, wc, fp, H_coeff)
    A_3d, Tau = _build_dde_arrays(diag_A, diag_C, B, wc.τ_0)

    println("=== Diagnostic for Ω=$Ω, τ_1=$τ_1 ===")
    println("diag(A) = ", round.(diag_A, digits=4))
    println("diag(C) = ", round.(diag_C, digits=4))
    println("Tau     = ", Tau)
    println("A[:,:,1] diagonal = ", round.(diag(A_3d[:, :, 1]), digits=4))
    println("A[:,:,2] diagonal = ", round.(diag(A_3d[:, :, 2]), digits=4))
    for k in 3:size(A_3d, 3)
        println("A[:,:,$k] (delay=$(Tau[k])) = ")
        display(round.(A_3d[:, :, k], digits=4))
    end

    # Quick eigenvalue check of just the A matrix (no delay, upper bound)
    ev = eigvals(A_3d[:, :, 1])
    println("Eigenvalues of A alone (no delay) = ", round.(ev, digits=4))
    println("Max Re(eig(A)) = ", maximum(real.(ev)))
    for (key, val) in B
        println("θ=$key = $val\n")
    end
end

"""
For a single wc (i.e. single τ_0 and τ_1) and ψ,
does the whole procedure to find the the continuation solution for the phase-locked solution of F(Ω,τ_1) = 0
and then saves it for matlab.

    If return_con == true, returns Ωs,τ_1s
"""
function stability_single_tau0(wc::WCParams; M::Int=30, n_neurons::Int=11, ψ::Float64=0.0, ϵ::Float64=0.01, return_con=false, filename::String="default", e0::Float64=0.8, i0::Float64=0.8, τ1_min::Float64=0.0, τ1_max::Float64=10.0, N_max::Int=5000)
    t_end = 200
    N = 2 * M + 1
    N_total = 2 * N
    fp = FourierParams(M, 2)
    τ_0 = wc.τ_0
    u0 = [e0, i0]
    t_span = (0.0, t_end)

    his(p, t) = u0

    #Orbit
    sol_dde1 = compute_orbit(u0, his, t_span, wc)
    X0_guess, T1_guess = find_init_guess(sol_dde1, N)
    X_sol, T, sol = compute_or_load_hb(X0_guess, T1_guess, fp, wc)
    A_iprc = compute_or_load_iprc(X_sol, T, fp, wc)
    A_phase = fp.S_p_inv * X_sol
    H_coeff = compute_H_coeff(wc, fp, A_iprc, A_phase)

    #Continuation
    ω = 2π / T
    params_Ω = ParamsΩ(ϵ, ω, wc.τ_0, ψ, n_neurons, fp, H_coeff)
    Ωs, τs, _, _ = continuation_phase_lock(params_Ω, wc;
        w_min=0.01, w_max=0.5,    # scale τ_1 range with τ_0
        v_min=0.01, v_max=5.0, p_min=τ1_min, p_max=τ1_max, N_max=N_max
    )

    if filename == "default"
        filename = "dde_batch_$(τ_0)_$(round(ψ;digits=2))_eps_$ϵ.mat"
    end

    batch_continuation_stability(Ωs[1:2:end], τs[1:2:end], ϵ, n_neurons, ψ, wc, fp, H_coeff, filename=filename)

    if return_con == true
        return Ωs, τs
    end
end

"""
Similar to stability_single_tau0 but computes the continuation of F(Ω,τ_1) = 0
for every τ_0 ∈ τ0_range. Recommended to use a range with nice and rounded numbers,
as many of the steps support loading precomputed data, so having τ_0 = 0.1,0.2 etc makes
it more likely to have been computed rather than having 0.1111. (To do this you can use e.g. τ0_range = 0.1:0.1:2.0 )
"""
function stability_mult_tau0(wc_params_dic::Dict, τ0_range; M::Int=30, n_neurons::Int=11, ψ::Float64=0.0, ϵ::Float64=0.01, filename::String="dde_batch_3d.mat")
    n_tau0 = length(τ0_range)

    # Outer cell arrays for MATLAB — filled per τ_0
    A_outer = Array{Any}(undef, 1, n_tau0)
    Tau_outer = Array{Any}(undef, 1, n_tau0)
    Ωs_per_tau0 = Vector{Vector{Float64}}(undef, n_tau0)
    τ1s_per_tau0 = Vector{Vector{Float64}}(undef, n_tau0)
    for (k, τ_0) in enumerate(τ0_range)
        println("\n=== τ_0 = $τ_0  ($k / $n_tau0) ===")
        wc_params_dic["tau_0"] = τ_0
        wc_k = WCParams(; convert_keys(wc_params_dic)...)

        e0 = 0.8
        i0 = 0.8
        t_end = 200
        N = 2 * M + 1
        N_total = 2 * N
        fp = FourierParams(M, 2)

        u0 = [e0, i0]
        t_span = (0.0, t_end)

        his(p, t) = u0

        #Orbit
        sol_dde1 = compute_orbit(u0, his, t_span, wc_k)
        X0_guess, T1_guess = find_init_guess(sol_dde1, N)
        X_sol, T, sol = compute_or_load_hb(X0_guess, T1_guess, fp, wc_k)
        A_iprc = compute_or_load_iprc(X_sol, T, fp, wc_k)
        A_phase = fp.S_p_inv * X_sol
        H_coeff = compute_H_coeff(wc_k, fp, A_iprc, A_phase)

        #Continuation
        ω = 2π / T
        params_Ω = ParamsΩ(ϵ, ω, wc_k.τ_0, ψ, n_neurons, fp, H_coeff)
        Ωs, τs, _, _ = continuation_phase_lock(params_Ω, wc_k;
            w_min=0.01, w_max=0.5,    # scale τ_1 range with τ_0
            v_min=0.01, v_max=5.0        # scale Ω range with natural frequency
        )
        Ωs_per_tau0[k] = Ωs
        τ1s_per_tau0[k] = τs

        # --- Build DDE matrices for each (Ω, τ_1) on this branch ---
        n_runs = length(Ωs)
        A_cells = Vector{Array{Float64,3}}(undef, n_runs)
        Tau_cells = Vector{Vector{Float64}}(undef, n_runs)

        Threads.@threads for i in 1:n_runs
            diag_A, diag_C, B = continuation_stability(
                Ωs[i], τs[i], ϵ, n_neurons, ψ, wc_k, fp, H_coeff
            )
            A_cells[i], Tau_cells[i] = _build_dde_arrays(diag_A, diag_C, B, τ_0)
        end

        A_row = Array{Any}(undef, 1, n_runs)
        Tau_row = Array{Any}(undef, 1, n_runs)
        for i in 1:n_runs
            A_row[i] = A_cells[i]
            Tau_row[i] = Tau_cells[i]'
        end

        A_outer[k] = A_row
        Tau_outer[k] = Tau_row
    end
    filename = joinpath("matlab_data", filename)
    matwrite(filename, Dict(
        "A_coeff" => A_outer,
        "Tau" => Tau_outer,
        "omegas" => [Ωs_per_tau0[k] for k in 1:n_tau0],   # cell array, no padding
        "tau1s" => [τ1s_per_tau0[k] for k in 1:n_tau0],
        "tau0s" => collect(τ0_range)',
        "n_neurons" => Float64(n_neurons),
        "epsilon" => ϵ,
    ))
    println("\nSaved → $filename")
end