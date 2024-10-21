module Operators

using OrderedCollections

using Main.tJmodel1D: System, Basis, State, makeBasis, getStateInfo

"""
    Superposition === OrderedDict{State, Complex{Float64}}

Compact representation of wave function of a single subspace. Alternative to sparse vector.
"""
Superposition = OrderedDict{State, Complex{Float64}}

"""
    SystemSuperposition === OrderedDict{System, Superposition}

Compact representation of wave function overlaping with many subspaces of the model. 
"""
SystemSuperposition = OrderedDict{System, Superposition}

""" 
    SystemWaveFunction === OrderedDict{System, Vector{Complex{Float64}}}

Represents a wave function that overlaps with many orthogonal subspaces of the model.
"""
SystemWaveFunction = OrderedDict{System, Vector{Complex{Float64}}} 

"""
    applyOperator(system::System, waveFunction::Vector{Complex{Float64}}, operator, args...) -> SystemWaveFunction

Apply `operator` with arguments `args...` to `waveFunction` of `system` subspace.
In general, operator action on arbitrary wave function of a single subspace will overlap with many subspaces of the whole system.
Return hash table `SystemWaveFunction` representing overlaping subspaces and corresponding resulting subspace wave functions.
"""
function applyOperator(system::System, waveFunction::Vector{Complex{Float64}}, operator, args...)::SystemWaveFunction
    basis = makeBasis(system)
    basisSet = OrderedDict{System, Basis}()
    
    result = SystemWaveFunction()
    for (state, index) in basis
        systemSuperposition = operator(args..., state = state, system = system)
        for (newSystem, superposition) in systemSuperposition
            if ~haskey(basisSet, newSystem)
                basisSet[newSystem] = makeBasis(newSystem)
                result[newSystem] = zeros(Complex{Float64}, length(basisSet[newSystem]))
            end
            for (newState, coeff) in superposition
                result[newSystem][basisSet[newSystem][newState]] += coeff * waveFunction[index]
            end
        end
    end

    return result
end


function Sk_z(k::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()
    
    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N
    ik::Complex{Float64} = 2.0 * pi * im * k / system.size

    # evaluate operator action on state and assign result to proper subspace
    siteValue = 1
    for R in 0:(system.size-1)
        if state.charges & siteValue > 0
            alpha = 0.5
            if state.magnons & siteValue > 0
                alpha = -0.5
            end
            if isodd(R)
                alpha = -alpha
            end
           
            newSystem = System(system, momentum = mod(p - k, N))
            # comment: p - k -> 2πi (p / N) - 2πi (2k / L)
        
            hasMomentum, _, _, _ = getStateInfo(state, newSystem)

            # assert periodicity-momentum match
            if hasMomentum
                coefficient = alpha * exp(-ik * R) / sqrt(system.size)

                # update result
                if ~haskey(result, newSystem)
                    result[newSystem] = Superposition()
                end
                if haskey(result[newSystem], state)
                    result[newSystem][state] += coefficient
                else
                    result[newSystem][state] = coefficient
                end
                
            end # if hasMomentum

        end # if state.charges

        siteValue = siteValue << 1
    end # for R

    return result
end # Sk_z


function Sr_z(r::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()
    
    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N

    # evaluate operator action on state and assign result to proper subspace
    for q in 0:(system.size-1)
        iq::Complex{Float64} = 2.0 * pi * im * q / system.size
        siteValue = 1
        for R in 0:(system.size-1)
            if state.charges & siteValue > 0
                alpha = 0.5
                if state.magnons & siteValue > 0
                    alpha = -0.5
                end
                if isodd(R)
                    alpha = -alpha
                end
               
                newSystem = System(system, momentum = mod(p - q, N))
                # comment: p - q -> 2πi (p / N) - 2πi (2q / L)
            
                hasMomentum, _, _, _ = getStateInfo(state, newSystem)

                # assert periodicity-momentum match
                if hasMomentum
                    coefficient = alpha * exp(-iq * (R - r)) / system.size

                    # update result
                    if ~haskey(result, newSystem)
                        result[newSystem] = Superposition()
                    end
                    if haskey(result[newSystem], state)
                        result[newSystem][state] += coefficient
                    else
                        result[newSystem][state] = coefficient
                    end
                    
                end # if hasMomentum

            end # if state.charges

            siteValue = siteValue << 1
        end # for R

    end # for q

    return result
end # Sk_r


function Sk_plus(k::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.spinsUp >= system.electrons)
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
            # calculate coefficient α_R^+ + β_R^+
            if state.magnons & siteValue > 0
                alpha = 1
            end
            if isodd(R)
                alpha = 1 - alpha
            end
          
            if alpha > 0
                # apply magnon annihilation/creation operators
                newState = State(state.charges, xor(state.magnons, siteValue))
 
                newSystem = System(system, spinsUp = system.spinsUp + 1, momentum = mod(p - k, N)) 
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
end # Sk_plus


function Sr_plus(r::Int64; state::State, system::System)
    result = SystemSuperposition()

    if (system.spinsUp >= system.electrons)
        return result
    end

    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N
    
    _, _, periodicity, _ = getStateInfo(state, system)

    # evaluate operator action on state and assign result to proper subspace
    for q in 0:(system.size-1)
        iq::Complex{Float64} = 2.0 * pi * im * q / system.size
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
                    # apply magnon annihilation/creation operators
                    newState = State(state.charges, xor(state.magnons, siteValue))
     
                    newSystem = System(system, spinsUp = system.spinsUp + 1, momentum = mod(p - q, N)) 
                    # comment: p - q -> 2πi (p / N) - 2πi (2q / L)

                    hasMomentum, repState, newPeriodicity, distance = getStateInfo(newState, newSystem)

                    # assert periodicity-momentum match
                    if hasMomentum
                        normalization = sqrt(periodicity / newPeriodicity) / system.size
                        phase = exp(-iq * (R - r) - (ip - 2 * iq) * distance)
                        
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

    end # for q

    return result
end # Sr_plus


function Sk_minus(k::Int64; state::State, system::System)::SystemSuperposition
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
                # apply magnon annihilation/creation operators
                newState = State(state.charges, xor(state.magnons, siteValue))
 
                newSystem = System(system, spinsUp = system.spinsUp - 1, momentum = mod(p - k, N)) 
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
end # Sk_minus


function Sr_minus(r::Int64; state::State, system::System)
    result = SystemSuperposition()

    if (system.spinsUp <= 0)
        return result
    end

    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N
    
    _, _, periodicity, _ = getStateInfo(state, system)

    # evaluate operator action on state and assign result to proper subspace
    for q in 0:(system.size-1)
        iq::Complex{Float64} = 2.0 * pi * im * q / system.size
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
                    # apply magnon annihilation/creation operators
                    newState = State(state.charges, xor(state.magnons, siteValue))
     
                    newSystem = System(system, spinsUp = system.spinsUp - 1, momentum = mod(p - q, N)) 
                    # comment: p - q -> 2πi (p / N) - 2πi (2q / L)

                    hasMomentum, repState, newPeriodicity, distance = getStateInfo(newState, newSystem)

                    # assert periodicity-momentum match
                    if hasMomentum
                        normalization = sqrt(periodicity / newPeriodicity) / system.size
                        phase = exp(-iq * (R - r) - (ip - 2 * iq) * distance)
                        
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

    end # for q

    return result

end # Sr_minus


function ck_up(k::Int64; state::State, system::System)::SystemSuperposition
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


function cr_up(r::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.spinsUp <= 0)
        return result
    end

    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N
    
    _, _, periodicity, _ = getStateInfo(state, system)

    # evaluate operator action on state and assign result to proper subspace
    for q in 0:(system.size-1)
        iq::Complex{Float64} = 2.0 * pi * im * q / system.size
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
     
                    newSystem = System(system, electrons = system.electrons - 1, spinsUp = system.spinsUp - 1, momentum = mod(p - q, N)) 
                    # comment: p - q -> 2πi (p / N) - 2πi (2q / L)

                    hasMomentum, repState, newPeriodicity, distance = getStateInfo(newState, newSystem)

                    # assert periodicity-momentum match
                    if hasMomentum
                        normalization = sqrt(periodicity / newPeriodicity) / system.size
                        phase = exp(-iq * (R - r) - (ip - 2 * iq) * distance)
                        
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

    end # for q

    return result
end # cr_up


function ck_down(k::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.spinsUp >= system.electrons)
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
 
                newSystem = System(system, electrons = system.electrons - 1, momentum = mod(p - k, N)) 
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
end # ck_down


function cr_down(r::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.spinsUp >= system.electrons)
        return result
    end

    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N
    
    _, _, periodicity, _ = getStateInfo(state, system)

    # evaluate operator action on state and assign result to proper subspace
    for q in 0:(system.size-1)
        iq::Complex{Float64} = 2.0 * pi * im * q / system.size
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

                    # assert periodicity-momentum match
                    if hasMomentum
                        normalization = sqrt(periodicity / newPeriodicity) / system.size
                        phase = exp(-iq * (R - r) - (ip - 2 * iq) * distance)
                        
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

    end # for q

    return result
end # cr_down


function ck_up_dag(k::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.electrons >= system.size)
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
        if state.charges & siteValue == 0
            # annihilating a hole
            newState = State(state.charges | siteValue, state.magnons)
            
            # creating a magnon if on sublattice B
            if isodd(R)
                newState.magnons |= siteValue
            end
            
            newSystem = System(system, electrons = system.electrons + 1, spinsUp = system.spinsUp + 1, momentum = mod(p - k, N)) 
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

        end # if state.charges

        siteValue = siteValue << 1
    end # for R

    return result
end # ck_up_dag


function cr_up_dag(r::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.electrons >= system.size)
        return result
    end

    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N
    
    _, _, periodicity, _ = getStateInfo(state, system)

    # evaluate operator action on state and assign result to proper subspace
    for q in 0:(system.size-1)
        iq::Complex{Float64} = 2.0 * pi * im * q / system.size
        siteValue = 1
        for R in 0:(system.size-1)
            alpha = 0
            if state.charges & siteValue == 0
                # annihilating a hole
                newState = State(state.charges | siteValue, state.magnons)
                
                # creating a magnon if on sublattice B
                if isodd(R)
                    newState.magnons |= siteValue
                end
                
                newSystem = System(system, electrons = system.electrons + 1, spinsUp = system.spinsUp + 1, momentum = mod(p - q, N)) 
                # comment: p - q -> 2πi (p / N) - 2πi (2q / L)

                hasMomentum, repState, newPeriodicity, distance = getStateInfo(newState, newSystem)

                # assert periodicity-momentum match
                if hasMomentum
                    normalization = sqrt(periodicity / newPeriodicity) / system.size
                    phase = exp(-iq * (R - r) - (ip - 2 * iq) * distance)
                    
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

            end # if state.charges

            siteValue = siteValue << 1
        end # for R

    end # for q

    return result
end # cr_up_dag


function ck_down_dag(k::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.electrons >= system.size)
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
        if state.charges & siteValue == 0
            # annihilating a hole
            newState = State(state.charges | siteValue, state.magnons)
            
            # creating a magnon if on sublattice A
            if iseven(R)
                newState.magnons |= siteValue
            end
            
            newSystem = System(system, electrons = system.electrons + 1, momentum = mod(p - k, N)) 
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

        end # if state.charges

        siteValue = siteValue << 1
    end # for R

    return result

end # ck_down_dag


function cr_down_dag(r::Int64; state::State, system::System)::SystemSuperposition
    result = SystemSuperposition()

    if (system.electrons >= system.size)
        return result
    end

    N::Int64 = system.size / 2
    p::Int64 = system.momentum
    
    ip::Complex{Float64} = 2.0 * pi * im * p / N
    
    _, _, periodicity, _ = getStateInfo(state, system)

    # evaluate operator action on state and assign result to proper subspace
    for q in 0:(system.size-1)
        iq::Complex{Float64} = 2.0 * pi * im * q / system.size
        siteValue = 1
        for R in 0:(system.size-1)
            alpha = 0
            if state.charges & siteValue == 0
                # annihilating a hole
                newState = State(state.charges | siteValue, state.magnons)
                
                # creating a magnon if on sublattice A
                if iseven(R)
                    newState.magnons |= siteValue
                end
                
                newSystem = System(system, electrons = system.electrons + 1, momentum = mod(p - q, N)) 
                # comment: p - q -> 2πi (p / N) - 2πi (2q / L)

                hasMomentum, repState, newPeriodicity, distance = getStateInfo(newState, newSystem)

                # assert periodicity-momentum match
                if hasMomentum
                    normalization = sqrt(periodicity / newPeriodicity) / system.size
                    phase = exp(-iq * (R - r) - (ip - 2 * iq) * distance)
                    
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

            end # if state.charges

            siteValue = siteValue << 1
        end # for R

    end # for q

    return result
end # cr_down_dag


### custom operators
include("./mods/custom.jl")


end # module Operators
