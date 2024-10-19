### --- Custom operators --- ###

# ------------------------------------------------------------------------------------------ #
#   Required Structure:                                                                      #
#       function operator_name(args...; state::State, system::System)::SystemSuperposition   #
#                              ^      ^--- required to split args... from kwargs...          #
#                              ^--- put as many arguments as you want (e.g. momentum k)      #
#                                                                                            #
#   Assumptions:                                                                             #
#       state::State ∈ Main.tJmodel1D.makeBasis(system::System)::Basis                       #
#       ^--- state::State is a single basis state from subspace defined by system::System    #
# ------------------------------------------------------------------------------------------ #


function ck_up_cq_down(k::Int64, q::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.spinsUp <= 0)
        return result
    end

    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N
    ik::Complex{Float64} = 2.0 * pi * im * k / system.size
    
    _, _, periodicity, _ = getStateInfo(state, system)

    # evaluate operator action on state and assign result to proper subspace
    siteValue = 1
    for R in 0:(system.size-1)
        alpha = 0
        if state.charges & siteValue > 0
            # calculate coefficient α_R^- + β_R^-
            if state.magnons & siteValue > 0
                alpha = 1
            end
            if iseven(R)
                alpha = 1 - alpha
            end
          
            if alpha > 0
                # creating a hole (include no hole & magnon constraint)
                newState = State(state.charges & ~siteValue, state.magnons & ~siteValue)
 
                newSystem = System(system, electrons = system.electrons - 1, spinsUp = system.spinsUp - 1, momentum = mod(p - k, N)) 
                # comment: p - k -> 2πi (p / N) - 2πi (2k / L)

                hasMomentum, repState, newPeriodicity, distance = getStateInfo(newState, newSystem)

                # assert periodicity-momentum match
                if hasMomentum
                    normalization = sqrt(periodicity / (newPeriodicity * system.size))
                    phase = exp(-ik * R - (ip - 2 * ik) * distance)
                    
                    coefficient = normalization * phase

                    # update result
                    if ~haskey(result, newSystem)
                        result[newSystem] = Superposition()
                    end
                    if haskey(result[newSystem], repState)
                        result[newSystem][repState] += coefficient
                    else
                        result[newSystem][repState] = coefficient
                    end
                    
                end # if hasMomentum

            end # if alpha

        end # if state.charges

        siteValue = siteValue << 1
    end # for R

    return result
end # ck_up

