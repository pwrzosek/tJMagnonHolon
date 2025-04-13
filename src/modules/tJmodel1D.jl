module tJmodel1D

using OrderedCollections
using SparseArrays
using KrylovKit
using JSON


"""
`struct System` immutable structure for system input parameters:
# Fields
*   `t::Float64`: value of the hopping constant t.
*   `J::Float64`: value of the coupling constant J. For ferromagnet `J < 0`, for antiferromagnet `J > 0`.
*   `λ::FLoat64`: parameter for scaling magnon-magnon interactions. For the pure t-J model `λ = 1.0`.
*   `α::Float64`: parameter for scaling XXZ anisotropy. For t-J model `α = 1.0`, for t-Jz model `α = 0.0`.
*   `size::Int64`: number of lattice sites.
*   `electrons::Int64`: number of electrons (must be in range `0:System.size` - Mott insulator constraint enforced).
*   `spinsUp::Int64`: index in range `0:electrons` indicating the number of spins up.
*   `momentum::Int64`: index of internal momentum sector in range `0:(System.size / 2 - 1)`. Index `momentum` corresponds to momentum `k = 2π * momentum / (System.size / 2)`.
"""
struct System
    t::Float64
    J::Float64
    λ::Float64
    α::Float64
    size::Int64
    electrons::Int64
    spinsUp::Int64
    momentum::Int64

    System(t, J, λ, α, size, electrons, spinsUp, momentum) = new(t, J, λ, α, size, electrons, spinsUp, momentum)
end

"""
    System(system::System = System(DEFAULT_ARGS...) [; t = system.t, J = system.J, ...,  momentum = system.momentum]) -> System

Create a new instance of `System` using `system::System` as a template with fields updated according to the given `kwargs`. Uses a default template if no template `system::System` is provided.
"""
System(
    system::System = System(
        1.0,    # DEFAULT_HOLE_HOPPING
        1.0,    # DEFAULT_SPIN_COUPLING
        1.0,    # DEFAULT_MAGNON_INTERACTION
        1.0,    # DEFAULT_ANISOTROPY
        16,     # DEFAULT_SYSTEM_SIZE
        16,     # DEFAULT_NUMBER_OF_ELECTRONS
        8,      # DEFAULT_NUMBER_OF_SPINS_UP
        0,      # DEFAULT_MOMENTUM_SUBSPACE
    ); 
    t           = system.t, 
    J           = system.J, 
    λ           = system.λ, 
    α           = system.α,
    size        = system.size, 
    electrons   = system.electrons, 
    spinsUp     = system.spinsUp, 
    momentum    = system.momentum
) = System(t, J, λ, α, size, electrons, spinsUp, momentum)


"""
`mutable struct State` mutable structure for storing the configuration of electrons and magnons in magnon-holon basis using binary representation
# Fields
*   `charges::Int64`: bit value 0 → hole,      1 → electron
*   `magnons::Int64`: bit value 0 → no magnon, 1 → magnon
"""
mutable struct State
    charges::Int64
    magnons::Int64
end


### comparison and hashing function (for State, see below)
import Base.isequal, Base.hash

"Compare magnon-holon basis states. Return `true` if `s1::State` and `s2::State` represent the same state, otherwise return `false`."
function isequal(s1::State, s2::State)
    return (s1.charges == s2.charges) && (s1.magnons == s2.magnons)
end

"Hashing function for magnon-holon basis states."
function hash(s::State)
    return hash(s.charges + ((1 << 32) * s.magnons))
end

"""
`Basis = OrderedDict{State, Int64}`

Basis is stored as a hash table with Robin Hood hashing algorithm. Key corresponds to state of the system and value contains position of that state in basis.
"""
Basis = OrderedDict{State, Int64}

"""
`mutable struct LinearCombination`: structure for storing result of Hamiltonian action on states from `basis::Basis`
# Fields
*   `state::State`: system configurations in binary representation written as decimal number
*   `coefficient::Vector{Complex{Float64}}`: coeffcient multiplying `state` in the linear combination
"""
mutable struct LinearCombination
    state::Vector{State}
    coefficient::Vector{Complex{Float64}}
end

"`Model = SparseMatrixCSC{Complex{Float64},Int64}`"
Model = SparseMatrixCSC{Complex{Float64},Int64}

"""
    run(input::Union{Missing, System} = missing; eigsolve = true, howmany = 1, kryldim = 30)

Run model diagonalization procedure. If no input provided, construct input using default arguments.
Return `(system::System, basis::Basis, model::Model, factorization)` where `factorization` is `nothing` if `eigsolve == false`.
If `eigsolve == true` then `factorization = (eigenvalues, eigenvectors, convergenceInfo)`.
Algorithm tries to converge at least `howmany` eigenvalues. 
Always print `convergenceInfo` to check accuracy since there is no warning if eigenvalues are poorly converged.
Keyword `kryldim` allows to specify the dimension of Krylov subspace for diagonalization algorithm (implicitely restarted Lanczos iteration).
"""
function run(input::Union{Missing, System} = missing; eigsolve = true, howmany = 1, kryldim = 30)
    system::System = if input === missing
        @info "No input detected. Applying default input arguments."
        System()
    else
        input
    end

    # run checks for input parameters
    checkSystem(system)
    
    println("\n", "Basis construction:")
    @time basis::Basis = makeBasis(system)
    
    println("\n", "Calculation of matrix coefficients:")
    @time model::Model = makeModel(basis, system)
    
    if eigsolve
        println("\n", "Factorization:")
        @time factorization = factorize(model, howmany = howmany, kryldim = kryldim) 
        return system, basis, model, factorization
    else
        return system, basis, model, nothing
    end
end

"Check if `system::System` describes a valid system. Throw `error` if any input parameter is invalid, otherwise return `nothing`."
function checkSystem(system::System)
    if (system.size < 2) || (system.size > 62) || isodd(system.size)
        error("`system size` must be even integer in range `[2, 62]`!")
    end
    if (system.electrons < 0) || (system.electrons > system.size)
        error("`number of electrons` must be integer in range `[0, system size]`!")
    end
    if (system.spinsUp < 0) || (system.spinsUp > system.electrons)
        error("`number of spins up` must be integer in range `[0, number of electrons]`!")
    end
    if (system.momentum < 0) || (2 * system.momentum >= system.size)
        error("`momentum subspace` must be integer in range `[0, (system size) / 2 - 1]`!")
    end
    return nothing 
end

"""
    sublatticeRotation(charges::Int64, spins::Int64, mask::Int64) -> Int64

Perform rotation of `spins` according to `mask` and return new `spins` configuration. Sites occupied by holes, according to `charges` parameter,  are not affected.
"""
function sublatticeRotation(charges::Int64, spins::Int64, mask::Int64)::Int64
    return xor(spins, mask) & charges
end

"""
    makeBasis(system::System) -> Basis

Return `Basis = OrderedDict{State, Int64}` hash table representing basis for given `System::system` parameters. 
Each index in `Basis` corresponds to `State`, and each value in `Basis` corresponds to position in the basis stored as `Int64`.
"""
function makeBasis(system::System)::Basis
    ### set sublattice rotation masks
    mask = sum(1 << k for k in 0:2:(system.size - 1))

    ### initialize the basis
    basis::Basis = Basis()

    ### get first state (i.e. with lowest index in binary basis)
    ### note: `1 << n == 2^n`, but former is faster
    charges::Int64 = system.electrons == 0 ? 0 : sum(n -> 1 << n, 0 : (system.electrons - 1)) # sum === (1 << system.electrons) - 1

    index = 0
    ### iterate over charge configurations
    while charges != -1
        spins::Int64 = system.spinsUp == 0 ? 0 : sum(n -> 1 << n, 0 : (system.spinsUp - 1))
        ### iterate over spin configurations
        while spins != -1
            ### locate spins at electrons and transform to magnon-holon basis
            magnons = sublatticeRotation(charges, locateSpins(charges, spins, system), mask)
            
            ### create state
            state = State(charges, magnons)

            ### obtain representative state and momentum match
            hasMomentum, representative, _, _ = getStateInfo(state, system)

            if hasMomentum
                if ~haskey(basis, representative)
                    push!(basis, representative => (index += 1))
                end
            end
            ### get next spin configuration
            spins = getNextConfiguration(spins, system.electrons)
        end
        ### get next charge configuration
        charges = getNextConfiguration(charges, system.size)
    end

    return basis
end

# helper function for makeBasis(system::System)
function getNextConfiguration(state::Int64, size::Int64)::Int64
    count = 0
    ### loop over <i,j> site pairs (bonds)
    for i in 0:(size - 2)
        j = i + 1

        ### get numeric value at i and j bit positions
        iValue, jValue = 1 << i, 1 << j

        ### check if there is bit 1 at site i
        if state & iValue > 0
            ### check if there is bit 1 at site j
            if state & jValue > 0
                ### clear bit at site i
                state &= ~iValue

                ### raise counter of cleared bits
                count += 1
            else
                ### swap bits at i and j
                state = xor(state, iValue + jValue)

                ### add cleared bits to lowest positions
                state += count == 0 ? 0 : sum(n -> 1 << n, 0:(count-1))

                ### return new state
                return state
            end
        end
    end
    ### if loop ends without returning new state we return -1
    ### as there is no next state in the requestes subspace
    return -1
end

# helper function for makeBasis(system::System)
function locateSpins(charges::Int64, spins::Int64, system::System)::Int64
    i = 0
    result = 0
    for j in 0:(system.size - 1) 
        jValue = 1 << j
        if charges & jValue > 0
            iValue = 1 << i
            if spins & iValue > 0
                result += jValue
            end
            i += 1
        end
    end
    return result
end

"""
    isGreater(a::State, b::State) -> Bool

Ordering function for basis states in magnon-holon representation. Returns `true` if `a > b`, and `false` otherwise.
"""
function isGreater(a::State, b::State)::Bool
    if (a.charges == b.charges) 
        return a.magnons > b.magnons
    end
    return a.charges > b.charges
end

"""
    getStateInfo(state::State, system::System) -> Tuple{Bool, State, Int64, Int64}

Return `(hasMomentum, representative, periodicity, distance)` where
*   `hasMomentum::Bool`: `true` if `state` matches momentum of `system`, and `false` otherwise
*   `representative::State`: representative state corresponding to `state`
*   `periodicity::Int64`: minimal `R` such that `state` shifted `2R` times with translation operator returns back to `state`
*   `distance::Int64`: distance between `state` and `representative` state in number of even translations
"""
function getStateInfo(state::State, system::System)::Tuple{Bool, State, Int64, Int64}
    ### initialize some constants for faster evaluation
    l::Int = system.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1

    ### initialize representative as state
    representative = state

    ### loop over translations of state
    newState = state
    distance::Int64 = 0
    periodicity::Int64 = 1
    while (
        newState = State(
            bitmov(newState.charges, l, false, hb = highestBit, hv = highestValue),
            bitmov(newState.magnons, l, false, hb = highestBit, hv = highestValue)
        );
        newState = State(
            bitmov(newState.charges, l, false, hb = highestBit, hv = highestValue),
            bitmov(newState.magnons, l, false, hb = highestBit, hv = highestValue)
        );
        (state.charges != newState.charges) || (state.magnons != newState.magnons)
    )
        if isGreater(representative, newState)
            representative = newState
            distance = periodicity
        end
        periodicity += 1
    end
 
    hasMomentum = rem(system.momentum * periodicity, system.size / 2) == 0
    return (hasMomentum, representative, periodicity, distance)
end

"""
    bitmov(s::Int, l::Int, f::Bool = false; hb::Int = 1 << (l - 1), hv::Int = (1 << l) - 1) -> Int

Arithmetic bit shift for calculating bit translations with periodic boundary conditions.

# Arguments
*   `s::Int` - value which binary representation will be shifted
*   `l::Int` - size of the cycle (total number of bits in the cycle)
*   `f::Bool` - `true ->` move forward (~mult by `2`); `false ->` backward (~div by `2`); ~mult and ~div mean multiplication and division with preservance of periodic boundary conditions within cycle size `l`
*   `hb::Int` [optional] - highest bit (for speed up put the value of 2^(l-1))
*   `hv::Int` [optional] - highest value (for speed put (2^l)-1)
"""
@inline bitmov(s::Int, l::Int, f::Bool = false; hb::Int = 1 << (l - 1), hv::Int = (1 << l) - 1) = f ? 2s - div(s, hb) * hv : div(s, 2) + rem(s, 2) * hb

"""
    hamiltonian(state::State, basis::Basis, system::System) -> LinearCombination

Apply Hamiltonian to `state` written in Sz momentum `basis` obtained for input `system` parameters. Returns `LinearCombination` representing resulting states with their coefficients.
"""
function hamiltonian(state::State, basis::Basis, system::System)::LinearCombination
    ### initialize result as linear combination
    result = LinearCombination(fill(state, system.size + 1), zeros(Complex{Float64}, system.size + 1))

    ### set sublattice rotation masks
    mask = sum(1 << k for k in 0:2:(system.size - 1))

    N::Int64 = system.size / 2

    ### check if initial state belongs to basis
    if haskey(basis, state)
        ### initialize some constants for faster evaluation
        l::Int = system.size
        highestBit::Int = 1 << (l - 1)
        highestValue::Int = (1 << l) - 1

        ### initialize ik for faster exponent calculations
        ik::Complex{Float64} = 2.0 * pi * im * system.momentum / N

        _, _, periodicity, _ = getStateInfo(state, system)

        charges = state.charges
        magnons = state.magnons

        ### loop over lattice sites
        for i in 1:system.size
            j = mod1(i + 1, system.size)

            ### get bit value at i and j bit positions
            iValue, jValue = (1 << (i - 1)), (1 << (j - 1))

            ### work out matrix coefficients [see Eq. (1.16) in PhD]

            iElectron, jElectron = div(charges & iValue, iValue), div(charges & jValue, jValue)

            if (iElectron == 1) && (jElectron == 1) # if sites i and j are occupied by electrons
                ### workout Hz
                iMagnon, jMagnon = div(magnons & iValue, iValue), div(magnons & jValue, jValue)

                coefficient = -1.0

                ### magnon cost 
                coefficient += iMagnon + jMagnon

                ### magnon-magnon interaction
                coefficient += -2.0 * system.λ * (iMagnon * jMagnon)

                coefficient *= 0.5 * system.J   

                result.coefficient[1] += coefficient

                ### workout Hxy (a * a + a^dagger * a^dagger terms) a.k.a. spin flip
                if iMagnon == jMagnon
                    ### perform spin flip === create/annihilate pair of magnons
                    newState = State(charges, xor(state.magnons, iValue + jValue))

                    ### get info about state after spin flip
                    hasMomentum, repState, newPeriodicity, distance = getStateInfo(newState, system)

                    ### check if it belongs to correct momentum subspace
                    ### if it does not, then ignore -- it will cancel out
                    ### with similar terms after summing over all the sites
                    if hasMomentum
                        ### calculate matrix coefficient
                        coefficient = 0.5 * system.J * system.α * exp(ik * distance) * sqrt(periodicity / newPeriodicity)

                        ### create a new entry in linear combination
                        ### and set its corresponding coeffcient
                        result.state[i + 1] = repState
                        result.coefficient[i + 1] = coefficient
                    end
                end
                
            elseif iElectron != jElectron # if one of the sites sites is occupied by a hole
                ### workout Ht terms
                newCharges = xor(charges, iValue + jValue) # hole hops between i and j

                ### hole annihilates or creates a magnon
                newMagnons = state.magnons
                iMagnon, jMagnon = div(newMagnons & iValue, iValue), div(newMagnons & jValue, jValue)
                if iMagnon == jMagnon # no magnons (other options are projcted out)
                    newMagnons = xor(newMagnons, iValue + jValue) & newCharges # hole creates a magnon behids
                else # magnon (only one, other options are projected out or taken care of before)
                    newMagnons = newMagnons & newCharges # the hole annihilates magnon
                end

                newState = State(newCharges, newMagnons)

                ### get info about state after electron/hole moved
                hasMomentum, repState, newPeriodicity, distance = getStateInfo(newState, system)

                ### check if it belongs to correct momentum subspace
                ### if it does not, then ignore -- it will cancel out
                ### with similar terms after summing over all the sites
                if hasMomentum
                    ### calculate matrix coefficient
                    coefficient = system.t * exp(ik * distance) * sqrt(periodicity / newPeriodicity)
                    # if fermion hops over the boundary add phase shift from anticommutation relations
                    if j < i && iseven(system.electrons)
                        coefficient = -coefficient
                    end

                    ### create a new entry in linear combination
                    ### and set its corresponding coeffcient
                    result.state[i + 1] = repState
                    result.coefficient[i + 1] = coefficient
                end
            end

        end
    end

    ### return resulting linear combination
    return result
end

"""
    makeModel(basis::Basis, system::System) -> Model

Calculate sparse matrix of the `model::Model` Hamiltonian, where `Model = SparseMatrixCSC{Complex{Float64},Int64}`.
"""
function makeModel(basis::Basis, system::System)::Model
    subspaceSize = length(basis)
    linearCombinationLength = system.size + 1
    I = Vector{Int64}(undef, linearCombinationLength * subspaceSize)
    J = Vector{Int64}(undef, linearCombinationLength * subspaceSize)
    V = Vector{Complex{Float64}}(undef, linearCombinationLength * subspaceSize)
    for (state, index) in basis
        linearCombination::LinearCombination = hamiltonian(state, basis, system)
        for it in 1:linearCombinationLength
            I[(index - 1) * linearCombinationLength + it] = index
            J[(index - 1) * linearCombinationLength + it] = basis[linearCombination.state[it]]
            V[(index - 1) * linearCombinationLength + it] = linearCombination.coefficient[it]
        end
    end
    return dropzeros!(sparse(I, J, V, subspaceSize, subspaceSize, +))
end

"""
    factorize(model::Model [; howmany = 1, which = :SR, kryldim = 30])

Compute eigenvalues (by default with smallest real part) and their corresponding eigenvectors.
"""
function factorize(model::Model; howmany = 1, which = :SR, kryldim = 30) 
    if length(model) != 0
        return eigsolve(model, howmany, which, ishermitian = true, krylovdim = kryldim)
    else
        return (missing, missing, missing)
    end
end

end
