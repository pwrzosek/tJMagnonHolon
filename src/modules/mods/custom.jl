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

### remove 2 electrons with momenta k, q and spin up, down respectively
function ck_up_cq_down(k::Int64, q::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.electrons < 2 || system.spinsUp == 0 || system.spinsUp == system.electrons)
        return result
    end

    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N
    iq::Complex{Float64} = 2.0 * pi * im * q / system.size
    ik::Complex{Float64} = 2.0 * pi * im * k / system.size
    
    _, _, periodicity, _ = getStateInfo(state, system)

    # evaluate operator action on state and assign result to proper subspace
    siteValue = 1
    for R in 0:(system.size-1)
        alpha = 0
        if state.charges & siteValue > 0
            # calculate coefficient α_R^+ + β_R^+
            if state.magnons & siteValue > 0
                alpha = 1
            end
            if isodd(R)
                alpha = 1 - alpha
            end
          
            if alpha > 0
                # creating a hole (include no hole & magnon constraint)
                newState = State(state.charges & ~siteValue, state.magnons & ~siteValue)
 
                newSystem = System(system, electrons = system.electrons - 1, momentum = mod(p - q, N)) 
                # comment: p - q -> 2πi (p / N) - 2πi (2q / L)

                hasMomentum, repState, newPeriodicity, distance = getStateInfo(newState, newSystem)

                # assert periodicity-momentum match for intermediate state
                if hasMomentum
                    state2 = repState

                    siteValue2 = 1
                    for R2 in 0:(system.size-1)
                        alpha = 0
                        if state2.charges & siteValue2 > 0
                            # calculate coefficient α_R^- + β_R^-
                            if state2.magnons & siteValue2 > 0
                                alpha = 1
                            end
                            if iseven(R2)
                                alpha = 1 - alpha
                            end
                          
                            if alpha > 0
                                # creating a hole (include no hole & magnon constraint)
                                newState2 = State(state2.charges & ~siteValue2, state2.magnons & ~siteValue2)
                 
                                newSystem2 = System(newSystem, electrons = newSystem.electrons - 1, spinsUp = newSystem.spinsUp - 1, momentum = mod(p - k - q, N)) 
                                # comment: p - k -> 2πi (p / N) - 2πi (2k / L)

                                hasMomentum, repState2, newPeriodicity2, distance2 = getStateInfo(newState2, newSystem2)

                                # assert periodicity-momentum match for final state
                                if hasMomentum
                                    normalization = sqrt(periodicity / newPeriodicity2) / system.size
                                    phase = exp(-iq * R - (ip - 2 * iq) * distance) * exp(-ik * R2 - (ip - 2 * iq - 2 * ik) * distance2)

                                    coefficient = normalization * phase

                                    # update result
                                    if ~haskey(result, newSystem2)
                                        result[newSystem2] = Superposition()
                                    end
                                    if haskey(result[newSystem2], repState2)
                                        result[newSystem2][repState2] += coefficient
                                    else
                                        result[newSystem2][repState2] = coefficient
                                    end
                                    
                                end # if hasMomentum

                            end # if alpha

                        end # if state2.charges

                        siteValue2 = siteValue2 << 1
                    end # for R2

                end # if hasMomentum

            end # if alpha

        end # if state.charges

        siteValue = siteValue << 1
    end # for R

    return result
end # ck_up_cq_down

