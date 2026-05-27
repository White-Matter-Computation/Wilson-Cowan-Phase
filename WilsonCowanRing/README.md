# WilsonCowanRing

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> WilsonCowanRing

It is authored by Emanuil.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "WilsonCowanRing"
```
which auto-activate the project and enable local path handling from DrWatson.

## Important Functions and Getting Started
In the `src` folder one can find all the code in different modules i.e. the code for the numerical continuation is in `src/continuation.jl`, the one of computing the iPRC is in `src/iprc.jl`, etc. In the `test` folder one can find examples and those are the files used to obtain all the plots. 
## Phase Reduced Model

The setup for the Wilson-Cowan model can be found in `src/model.jl` and parameters are saved in the struct `WCParams`. Example for parameters can be found in `params.toml` then one can load them up as 
```julia
all_params = TOML.parsefile("params.toml")
wc = WCParams(; convert_keys(all_params["set1"])...)
```
If one wishes to obtain the continuation results only, the function `stability_single_tau0` can be used for a single $\tau_0$ or `stability_mult_tau0(params1, τ0_range)` for multiple $\tau_0$, e.g. 
```julia
all_params = TOML.parsefile("params.toml")
params1 = all_params["set1"]
wc1 = WCParams(; convert_keys(all_params["set1"])...)
n_neurons = 11
ψ_0 = 0.0

stability_single_tau0(wc1; n_neurons=n_neurons, ψ=ψ_0)


#Multiple values
#Using round numbers for τ0 is recommended, as it is easier to make use of previously saved data
τ0_range = 0.1:0.1:2.0 

stability_mult_tau0(params1, τ0_range)
```
Those functions will save the results from contuations in a `.mat` file that can be loaded up in MATLAB using `DDEBifTool` (see the `matbad` folder outside of this project).
The functions have the following optional arguments 
```julia
#default will use filename = "dde_batch_$(τ_0)_$(round(ψ;digits=2)).mat"
stability_single_tau0(wc::WCParams; M::Int=30, n_neurons::Int=11, ψ::Float64=0.0, ϵ::Float64=0.01, return_con=false, filename::String="default", e0::Float64=0.8, i0::Float64=0.8, τ1_min::Float64=0.0, τ1_max::Float64=10.0, N_max::Int=5000)

stability_mult_tau0(wc_params_dic::Dict, τ0_range; M::Int=30, n_neurons::Int=11, ψ::Float64=0.0, ϵ::Float64=0.01, filename::String="dde_batch_3d.mat")
```
### Manual Stability
Alternatively one can do this manually by computing the iPRC and thus find the coefficients of the interaction function $H$ in 
```math
\frac{d\theta_i}{dt}=\omega+\varepsilon\sum_{j=1}^Nw_{ij} H(\theta_i(t-\tau_0)-\theta_i(t),\theta_j(t-\tau_{ij})-\theta_i(t))
```
Example
```julia
#Guess starting point
e0 = 0.8
i0 = 0.8
t_end = 200
#Fourier Mods
M  = 30
N = 2 * M + 1
N_total = 2 * N
fp = FourierParams(M, 2)
τ_0 = wc.τ_0
u0 = [e0, i0]
t_span = (0.0, t_end)

#Can use constant history for this simple case
his(p, t) = u0

#Orbit
sol_dde1 = compute_orbit(u0, his, t_span, wc)
#Obtain the Orbit
X0_guess, T1_guess = find_init_guess(sol_dde1, N)
#Use Harmonic Balance to obtain an expression for the limit cycle - X_sol
X_sol, T, sol = compute_or_load_hb(X0_guess, T1_guess, fp, wc)
#Fourier Coefficients for the iPRC
A_iprc = compute_or_load_iprc(X_sol, T, fp, wc)
#Fourier Coefficients for the Limit Cycle
A_phase = fp.S_p_inv * X_sol
#Compute the fourier coefficients of H 
H_coeff = compute_H_coeff(wc, fp, A_iprc, A_phase)

 ω = 2π / T
params_Ω = ParamsΩ(ϵ, ω, wc.τ_0, ψ, n_neurons, fp, H_coeff)
Ωs, τs, _, _ = continuation_phase_lock(params_Ω, wc;
        w_min=0.01, w_max=0.5,    # small τ_1 values
        v_min=0.01, v_max=5.0        # scale Ω range
    )

#Or alternatively one can use a contour plot e.g.
contour_phase_lock(params_Ω1, wc1, v_min=1.0, v_max=3.0, w_min=0.0, w_max=5.0)
```
### Continuation
The `continuation_phase_lock` function has the following default values, where `p_min` and `p_max` can be used to control how far your parameter can vary in continuation.
```julia
continuation_phase_lock(params_Ω::ParamsΩ, wc::WCParams;
    h_min=0.00001, h_max=0.01, h_con=0.001, N_max::Int=5000,
    t_0=nothing, Ω_guess=nothing,      # now optional
    v_min=0.01, v_max=4.5,
    w_min=0.01, w_max=0.5, p_min=-Inf, p_max=Inf
)
```
In case no guess is provided the function `find_initial_guess_continuation` is called, which uses a method similar to the countour plot where we look for sign changes on the grid provided i.e. 
```julia
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
```
For more on continuation examine  `src/continuation.jl`.

## Full Model
```math
\tau_e\frac{du_i(t)}{dt} = -u_i(t)  +f \left( I_u + w_{uu}u_i(t-\tau_0) - w_{vu}v_i(t-\tau_0) + \varepsilon \sum_{j=1}^Nw_{ij}u_j(t-\tau_{ij})\right)
```
```math
\tau_i\frac{v_i(t)}{dt} = -v_i(t) + f(I_v+w_{uv}u_i(t-\tau_0) - w_{vv}v_i(t-\tau_0)),
```
Now one can run the full model by just calling `run_full_network` which will give you the $u,v,t$ e.g.
```julia
n_neurons = 11
wave_num = 1
T1_guess = 3.33

#Pick  epsilon
ϵ1 = 0.1

#Solve the DDE network model
u1_stab1, v1_stab1, t1_stab1 = run_full_network(ϵ1, wc1, 2.5, n_neurons, T1_guess, wave_num)
```
The default values can be seen below (`tend_mult` indicates how many periods will be simulated) and it now support a pertubation `init_pertub` to be added to the history of each oscillator 
```julia
function run_full_network(ϵ::Float64, wc::WCParams, τ_1::Float64, n_neurons::Int, T_guess::Float64, wave_num::Int; T_sample=2000, tend_mult=100, X0_guess=[0.4, 0.6], abstol=1e-9, reltol=1e-8, plot_network=true, init_pertub=nothing)
```
If `plot_network=true` a phase plot, including starting and final points, as well as the history function, will be produced. After this one can make an ST(Space Time) plot by 
```julia
#For a single ST plot over the whole time interval
plot_stefan(u1_stab1, v1_stab1, t1_stab1, 2.5, wc1, n_neurons, wave_num, ϵ1)
#To make two sepated ST plots over the first and last 20% (better for larger t)
plot_stefan2(u1_stab1, v1_stab1, t1_stab1, 2.5, wc1, n_neurons, wave_num, ϵ1)
```
