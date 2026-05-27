using Plots
using DrWatson

function plot_orbit(sol_dde, p::WCParams)
    #es = first.(sol_dde.u)
    #is = last.(sol_dde.u)
    es = sol_dde[1, :]
    is = sol_dde[2, :]
    plot(sol_dde, label=["E(t)" "I(t)"], linewidth=2)
    params = paramsdict(p)
    filename = savename(params, "svg")
    savefig(plotsdir(join(["orbit_ei", filename])))
    plot(is, es, xlabel="I", ylabel="E", linewidth=2)
    savefig(plotsdir(join(["orbit_phase", filename])))
end

function plot_init_guess(X0, p::WCParams)
    X0_phs_reshape = reshape(X0, 2, :)
    xs0 = X0_phs_reshape[1, :]
    ys0 = X0_phs_reshape[2, :]

    params = paramsdict(p)
    filename = savename(params, "svg")

    plot(xs0, ys0)
    xlabel!("E(t)")
    ylabel!("I(t)")
    savefig(plotsdir(join(["X0_guess", filename])))
end




function plot_iprc(A_iprc, T, p::FourierParams, wc::WCParams)
    N = p.N
    M = p.M
    t_range = [n * T / N for n in (-M:M)]
    n_dense = 1000
    t_dense = range(t_range[1], t_range[end], length=n_dense)

    q_x_dense = Vector{Float64}(undef, n_dense)
    q_y_dense = Vector{Float64}(undef, n_dense)

    for (i, t) in enumerate(t_dense)
        q_x_dense[i], q_y_dense[i] = fourier_recon(t, A_iprc, T, p)
    end

    params = paramsdict(wc)
    filename = savename(params, "svg")

    θ_dense = 2π .* (t_dense ./ T) #(-π,π)
    plot(θ_dense, q_x_dense, linewidth=2, label="IPRC x")
    plot!(θ_dense, q_y_dense, linewidth=2, label="IPRC y")
    xlabel!("θ")
    ylabel!("iPRC")
    title!("iPRC Dense Plot")
    xlabel!("θ")
    ylabel!("iPRC")
    savefig(plotsdir(join(["iPRC", filename])))
end

function plot_continuation(xs, ps, X0, V, params_Ω::ParamsΩ, wc::WCParams)
    pname = plotsdir("continuation_omega$(params_Ω.ω)_psi_$(params_Ω.ψ)_n_neurons_$(params_Ω.n_neurons)" * savename(wc) * ".svg")

    plot(ps, xs, lw=2, linestyle=:dash, c=:blue)
    quiver!([X0[2]], [X0[1]],
        quiver=([V[2]], [V[1]]),
        label="initial tangent")
    xlabel!("τ_1")
    ylabel!("Ω")
    title!("Continuation plot for ψ = $(round(params_Ω.ψ, digits=2)), N = $(params_Ω.n_neurons) with τ_0 = $(params_Ω.τ_0)")
    savefig(pname)
end

function plot_cont_cont(params_Ω::ParamsΩ, wc::WCParams, Ωs_cont, τs_cont; v_min=2.0, v_max=4.5, w_min=-5.0, w_max=5.0)
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


    pname = plotsdir("cont_cont_omega$(round(ω,digits=2))_psi_$(round(ψ,digits=2))_n_neurons_$(n_neurons)_eps_$(ϵ)_vmin_$(v_min)_v_max_$(v_max)_w_min_$(w_min)_w_max_$(w_max)" * savename(wc) * ".svg")
    contour(ws, vs, Z; levels=[0], colorbar=false, label="Contour plot", c=:red)

    plot!(τs_cont, Ωs_cont, lw=2, linestyle=:dash, c=:green, alpha=0.5, label="Continuation")
    xlims!(0, w_max)
    ylims!(0.8 * v_min, 1.2 * v_max)
    xlabel!("τ1")
    ylabel!("Ω")
    title!("Plot for ψ = $(round(ψ,digits=3)), N = $n_neurons \n with τ0 = $(wc.τ_0), ϵ=$(ϵ)")
    savefig(pname)

end

function plot_stefan_phase(u, t, τ_1::Float64, wc::WCParams, n_neurons::Int, wave_num::Int, ϵ::Float64)
    pname_u = plotsdir("u_stefan_phase$(wave_num)_$(n_neurons)_$(τ_1)_$ϵ" * savename(wc) * ".svg")


    heatmap(t, 1:size(u, 1), u,
        xlabel="t",
        ylabel="i",
        color=:balance,
        #clims=(-maximum(abs, u), maximum(abs, u)),
        colorbar_title="u_i(t)")

    savefig(pname_u)
end


function plot_stefan(u, v, t, τ_1::Float64, wc::WCParams, n_neurons::Int, wave_num::Int, ϵ::Float64)

    pname_u = plotsdir("u_stefan$(wave_num)_$(n_neurons)_$(τ_1)_$ϵ" * savename(wc) * ".svg")
    pname_v = plotsdir("v_stefan$(wave_num)_$(n_neurons)_$(τ_1)_$ϵ" * savename(wc) * ".svg")


    heatmap(t, 1:size(u, 1), u,
        xlabel="t",
        ylabel="i",
        color=:balance,
        #clims=(-maximum(abs, u), maximum(abs, u)),
        colorbar_title="u_i(t)")

    savefig(pname_u)

    heatmap(t, 1:size(v, 1), v,
        xlabel="t",
        ylabel="i",
        color=:balance,
        colorbar_title="v_i(t)")

    savefig(pname_v)

end

function plot_stefan2(u, v, t, τ_1, wc::WCParams, n_neurons::Int, wave_num::Int, ϵ::Float64; max_plot::Union{Nothing,Float64}=nothing)
    pname_u = plotsdir("u_stefan2_$(wave_num)_$(n_neurons)_$(τ_1)_$ϵ" * savename(wc) * ".svg")
    pname_v = plotsdir("v_stefan2_$(wave_num)_$(n_neurons)_$(τ_1)_$ϵ" * savename(wc) * ".svg")

    if isnothing(max_plot)
        n = length(t)
        i_20 = round(Int, 0.2 * n)

        idx_early = 1:i_20
        idx_late = (n-i_20):n

        early_title = "Early (0–20%)"
        late_title = "Late (80–100%)"

    else
        idx_early = findall(t .<= max_plot)
        idx_late = findall(t .>= (t[end] - max_plot))
        early_title = "Early (0–$(max_plot))"
        late_title = "Late ($(round(t[end]-max_plot, digits=2))–$(round(t[end], digits=2)))"
    end

    for (pname, data, label) in [(pname_u, u, "u"), (pname_v, v, "v")]
        p1 = heatmap(
            t[idx_early], 1:n_neurons, data[:, idx_early],
            xlabel="t",
            ylabel="i",
            title=early_title,
            color=:balance,
            colorbar=false
        )

        p2 = heatmap(
            t[idx_late], 1:n_neurons, data[:, idx_late],
            xlabel="t",
            ylabel="i",
            title=late_title,
            color=:balance,
            colorbar_title="$(label)_i(t)"
        )

        plot(
            p1, p2,
            layout=(1, 2),
            size=(1200, 400),
            plot_title="τ_1=$τ_1, τ_0=$(wc.τ_0), wave num=$wave_num, ϵ=$ϵ"
        )

        savefig(pname)
    end

end
