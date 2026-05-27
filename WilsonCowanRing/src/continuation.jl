using LinearAlgebra
"""
    cont(F, x0, p0, h, N, V; kwargs...)

Pseudo-arclength continuation method.

# Arguments
- `F(x,p)`: nonlinear function defining F(x,p) = 0
- `x0`: initial solution (scalar)
- `p0`: initial parameter (scalar)
- `h`: step size
- `N`: number of continuation steps
- `V`: initial tangent vector (normalized)

# Keyword Arguments
- `params`: optional extra parameters passed to `F`
- `tol`: Newton tolerance (default 1e-8)
- `adaptive_step`: enable adaptive step size
- `dFdx`: derivative w.r.t x (optional)
- `dFdp`: derivative w.r.t p (optional)
- `h_der`: finite difference step (default 1e-8)
- `h_min`: minimum h if adaptive_step
- `h_max`: maximum h if adaptive_step
- `newton_maxitr`: maximum itteration for Newton (default 10)
- `output_diverge` : return result even if does not converge(default false)

# Returns
Continuation branch data xs,ps where xs = Vector{Any}; ps = Vector{Any}
"""
function cont(F, x0, p0, h::Float64, N::Int, V; params=nothing, tol=1e-8, adaptive_step=false, dFdx=nothing, dFdp=nothing, h_der=1e-8, h_min=nothing, h_max=nothing, newton_maxitr=10, output_diverge=false, p_min=-Inf, p_max=Inf)
    #F(x,p,params) is params != nothing else F(x,p)
    #If derivative is not provided do central difference 
    if adaptive_step
        if h_min === nothing || h_max === nothing
            error("If using adaptive step provide a h_min and h_max")
        end
    end


    Fwrap = params === nothing ?
            (x, p) -> F(x, p) :
            (x, p) -> F(x, p, params)

    # Default central differences
    if dFdx === nothing
        dFdx_wrap = (x, p) -> (Fwrap(x + h_der, p) - Fwrap(x - h_der, p)) / (2h_der)
    else
        #Make sure that if params are provided we call dFdx(x,p,params), note we do not need to be careful in the Central Difference case as Fwrap takes care of it 
        dFdx_wrap = params === nothing ?
                    (x, p) -> dFdx(x, p) :
                    (x, p) -> dFdx(x, p, params)
    end

    # Wrap dFdp
    if dFdp === nothing
        dFdp_wrap = (x, p) -> (Fwrap(x, p + h_der) - Fwrap(x, p - h_der)) / (2h_der)
    else
        dFdp_wrap = params === nothing ?
                    (x, p) -> dFdp(x, p) :
                    (x, p) -> dFdp(x, p, params)
    end

    #Xs[i] = [x_i,p_i] 
    X = [x0; p0]
    Xs = Vector{typeof(X)}()
    for n in 1:N
        push!(Xs, copy(X))

        X_guess = Xs[n] + h * V

        #for n = 1 X≡X1 = newton(X_guess,X0)
        if adaptive_step == false
            X = newton_itr(X_guess, Xs[n], h, V, Fwrap, dFdx_wrap, dFdp_wrap, adaptive_step, tol, h_min, h_max, newton_maxitr)
        elseif adaptive_step == true
            X, h = newton_itr(X_guess, Xs[n], h, V, Fwrap, dFdx_wrap, dFdp_wrap, adaptive_step, tol, h_min, h_max, newton_maxitr)
        else
            error("Invalid value for adaptive_step = $adaptive_step")
        end

        #X=nothing means newton did not converge
        if isnothing(X) || X[2] < p_min || X[2] > p_max
            if output_diverge
                println("Did not converge after $n continuation steps, but results are still returned")
                xs = first.(Xs)
                ps = last.(Xs)

                return xs, ps
            else
                error("Did not converge after $n continuation steps")
            end
        end

        #for n=1 W = X1-x0
        W = X - Xs[n]
        W /= norm(W)

        # if X[2] > 62
        #     println("\nX is $X and V = $V and W = $W and h =$h\n")
        # end

        if norm(W) < eps()
            error("Zero tangent encountered")
        end

        if dot(V, W) < 0
            println("Opposite direction or what ")
            W = -W
        end

        # if sign(W[2]) != sign(V[2])
        #     println("Parameter direction reversed at p = $(X[2]), likely near a fold")
        #     push!(Xs, copy(X))
        #     break
        # end

        V = W

    end

    push!(Xs, copy(X))
    xs = first.(Xs)
    ps = last.(Xs)

    return xs, ps
end

function newton_itr(X_guess, X0, h, V, F, dFdx, dFdp, adaptive_step, tol, h_min, h_max, maxitr)
    #No adaptive step 
    if adaptive_step == false
        for i in 1:maxitr
            G1 = F(X_guess[1], X_guess[2])
            G2 = dot(V, X_guess - X0) - h
            G = [G1; G2]

            if norm(G) < tol
                return X_guess
            end

            J = [
                dFdx(X_guess[1], X_guess[2]) dFdp(X_guess[1], X_guess[2]);
                V[1] V[2]
            ]

            #If ODE can use J to determine stability 
            # A \ b solves Ax = b 
            X_guess += J \ (-G)

            if i == maxitr
                println("Did not converge, after $maxitr newton iterations the norm is $(norm(G))")
                return nothing
                #error("Did not converge, after $maxitr itterations the norm is $(norm(G))")
            end

        end
        error("Got out of for loop")
    elseif adaptive_step == true
        X_predictor = copy(X_guess)  # save before iterating

        for c in 1:maxitr
            G1 = F(X_guess[1], X_guess[2])
            G2 = dot(V, X_guess - X0) - h
            G = [G1; G2]

            if norm(G) < tol
                #EXTRA TEST 
                #Idea is if the Newton makes too big of a jump it is probably going to a different branch so we reduce h if that happens
                if norm(X_guess - X_predictor) > 3 * h && h > h_min
                    println("Newton jumped too far (dist=$(norm(X_guess - X_predictor)), h=$h), halving h")
                    if h <= h_min
                        println("Already at h_min, cannot reduce further")
                        return nothing, nothing
                    end
                    new_h = max(h / 2, h_min)
                    X_new_guess = X0 + new_h * V
                    return newton_itr(X_new_guess, X0, new_h, V, F, dFdx, dFdp, adaptive_step, tol, h_min, h_max, maxitr)
                end

                #If fast converging increase h
                if c < 3 && 1.3 * h < h_max
                    return X_guess, 1.3 * h
                else
                    return X_guess, h
                end
            end

            J = [
                dFdx(X_guess[1], X_guess[2]) dFdp(X_guess[1], X_guess[2]);
                V[1] V[2]
            ]

            #If ODE can use J to determine stability 
            # A \ b solves Ax = b 
            X_guess += J \ (-G)

            if c == maxitr
                # == should suffice but it is a good practice to eer to side of caution in case somehow h turns out to be h_min - ε
                if h <= h_min
                    println("Did not converge, after $maxitr newton iterations the norm is $(norm(G)) and h is $h and X_guess is $X_guess")
                    return nothing, nothing
                else
                    new_h = max(h / 2, h_min)
                    #Update X_guess to as it was off 
                    X_new_guess = X0 + new_h * V
                    return newton_itr(X_new_guess, X0, new_h, V, F, dFdx, dFdp, adaptive_step, tol, h_min, h_max, maxitr)
                end
            end
        end
        error("Got out of for loop")

    else
        error("Invalid value for adaptive_step = $adaptive_step")
    end
end
