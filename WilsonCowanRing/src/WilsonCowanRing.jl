module WilsonCowanRing


include("model.jl")
include("limitcycle.jl")
include("harmonic_balance.jl")
include("iprc.jl")
include("helper.jl")
include("h_fun.jl")
include("continuation.jl")
include("stability.jl")
include("plots.jl")
include("network.jl")

export WCParams, FourierParams, HBParams, f, f_prime, F1, wc_rhs, wc_dde!, compute_orbit, find_init_guess, plot_orbit, plot_init_guess, paramsdict, FourierParams, compute_hb_orbit, compute_or_load_hb, compute_iprc, compute_or_load_iprc, plot_iprc, wc_network, find_network_guess, compute_G_coeff, compute_H_coeff, compute_H_fourier, compute_H_integral, ParamsΩ, cont, contour_phase_lock, continuation_phase_lock, plot_continuation, plot_stefan, continuation_stability, save_for_ddebiftool, batch_continuation_stability, diagnose_single_point, plot_cont_cont, plot_stefan2, run_full_network, plot_stefan_phase, wc_reduced, stability_mult_tau0, stability_single_tau0, compute_G_coeff2, w_k, comp
end
