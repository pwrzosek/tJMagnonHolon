module SpectralFunction

using OrderedCollections
using LinearAlgebra
using SparseArrays
using KrylovKit

using Main.tJmodel1D: Model

"""
`struct Krylov` immutable structure for passing parameters to Krylov subspace methods.
# Fields
*   `dimension::Int64` - maximum dimension of Krylov subspace.
"""
struct Krylov
    dimension::Int64
end

"""
    run(ωRange, δ, system, ψ, operator, args...; krylovDimension = 400, returnGreensFunction = false) -> Union{Vector{Float64}, Vector{ComplexF64}}

Apply `operator` with arguments `args...` to `system` wave function `ψ`. Then calculate spectral function over `ωRange` with artificial broadening `δ` using Lanczos tridiagonalization with initial state `Ô | ψ >`.
Return Vector{Float64} of size `length(ωRange)` with calculated values of the spectral function. A keyword argument `krylovDimension` defines maximum dimension of Krylov subspace for Lanczoas tridiagonalization algorithm. Typically `krylovDimension` ∈ [100, 1000] gives best performance.
Small values of `krylovDimension` speed up the calculation process but reduce the quality of spectral function. It is not advised to set Krylov dimension too high - performance eventually drops due to orthogonality loss.
Return Greens function instead if keyword argument `returnGreensFunction` is set to `true`.
"""
function run(ωRange, δ, system, ψ, operator, args...; krylovDimension = 400, returnGreensFunction::Bool = false)
    println("Evaluating step for operator args: ", args...)    
    spectrum = if returnGreensFunction
        zeros(ComplexF64, length(ωRange))
    else
        zeros(Float64, length(ωRange))
    end
    
    println("> Applying operator to ψ...")    
    @time initialStateSystem = Main.Operators.applyOperator(system, ψ, operator, args...)
    
    println("> Calculating spectral function...")    
    @time for (system, initialState) in initialStateSystem
        basis = Main.tJmodel1D.makeBasis(system)
        model = Main.tJmodel1D.makeModel(basis, system)

        spectrum += Main.SpectralFunction.calculate(ωRange, 1.0im * δ, initialState, model, krylovDimension, returnGreensFunction)
    end
    println()

    return spectrum
end


"""
    calculate(ωRange::Vector{Float64}, iDelta::Complex{Float64}, initialState::Vector{Complex{Float64}}, model::Model, krylovDimension::Int64 = 400, returnGreensFunction = false)

Run Lanczos tridiagonalization of a subspace `model` starting from `initialState` including up to `krylovDimension` dimensions.
Return spectral function calculated over range `ωRange` with data broadening `iDelta`.
"""
function calculate(ωRange::Vector{Float64}, iDelta::Complex{Float64}, initialState::Vector{Complex{Float64}}, model::Model, krylovDimension::Int64 = 400, returnGreensFunction::Bool = false)
    krylov = Krylov(krylovDimension)
    diagonal, offDiagonal, size = calculateLanczos(initialState, model::Model, krylov::Krylov)
    if returnGreensFunction 
        return greensFunction(ωRange, iDelta, initialState, diagonal, offDiagonal, size)
    end
    return spectralFunction(ωRange, iDelta, initialState, diagonal, offDiagonal, size)
end

"""
    calculateLanczos(initialState::Vector{Complex{Float64}}, model::Model, krylov::Krylov)

Calculate tridiagonal form of hermitian `model` for given `krylov` parameters of Krylov subspace starting from `initialState`.
Return compact form `(diagonal, offDiagonal, size)` of symmetric tridiagonal matrix as two vectors of `diagonal` and `offDiagonal` terms with sizes `size` and `size-1` respectively. 
"""
function calculateLanczos(initialState::Vector{Complex{Float64}}, model::Model, krylov::Krylov)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    diagonal = Float64[]
    offDiagonal = Float64[]
    krylovSpaceDimensions = 1

    ### normalization of the inital state
    initialNorm = norm(initialState)
    if initialNorm == 0.0
        return ([0.0], [], 1)
    end

    state1 = initialState / initialNorm
    state2 = model * state1

    push!(diagonal, real(dot(state1, state2)))
    axpy!(-diagonal[end], state1, state2)
    push!(offDiagonal, Float64(norm(state2)))

    maxDimension = min(krylov.dimension, length(initialState))

    while krylovSpaceDimensions < maxDimension
        if offDiagonal[end] == 0.0
            break
        end

        state1 .*= -offDiagonal[end]
        state2 ./= offDiagonal[end]
        state1, state2 = state2, state1

        state2 += model * state1
        push!(diagonal, real(dot(state1, state2)))
        axpy!(-diagonal[end], state1, state2)
        push!(offDiagonal, Float64(norm(state2)))

        krylovSpaceDimensions += 1
    end

    return (diagonal, offDiagonal[1:(end-1)], krylovSpaceDimensions)
end

"""
    greensFunction(ω::Complex{Float64}, initialState, diagonal, offDiagonal, size) -> Complex{Float64}

Calculate the value of Greens function for given `ω` point.
Internally evaluates a continued fraction of depth `size` based on results `diagonal` and `offDiagonal` of Lanczos tridiagonalization procedure.
The result is re-scaled by the norm of `initialState`.
"""
function greensFunction(ω::Complex{Float64}, initialState, diagonal, offDiagonal, size)::Complex{Float64}
    result = ω - diagonal[size]
    for it in (size - 1):-1:1
        result = ω - diagonal[it] - offDiagonal[it]^2 / result
    end
    return dot(initialState, initialState) / result
end

"""
    greensFunction(ωRange::Vector{Float64}, iDelta::Complex{Float64}, initialState::Vector{Complex{Float64}}, diagonal::Vector{Float64}, offDiagonal::Vector{Float64}, size::Int64) -> Vector{Complex{Float64}}

Calculate Greens function over range of points `ωRange` applying data broadening `iDelta = iδ`.
"""
function greensFunction(ωRange::Vector{Float64}, iDelta::Complex{Float64}, initialState::Vector{Complex{Float64}}, diagonal::Vector{Float64}, offDiagonal::Vector{Float64}, size::Int64)::Vector{Complex{Float64}}
    result = Vector{ComplexF64}(undef, length(ωRange))
    for (it, ω) in enumerate(ωRange)
        result[it] = greensFunction(ω + iDelta, initialState, diagonal, offDiagonal, size)
    end
    return result
end

"""
    spectralFunction(ωRange::Vector{Float64}, iDelta::Complex{Float64}, initialState::Vector{Complex{Float64}}, diagonal::Vector{Float64}, offDiagonal::Vector{Float64}, size::Int64) -> Vector{Float64}

Calculate spectral function over range of points `ωRange` applying data broadening `iDelta = iδ`.
"""
function spectralFunction(ωRange::Vector{Float64}, iDelta::Complex{Float64}, initialState::Vector{Complex{Float64}}, diagonal::Vector{Float64}, offDiagonal::Vector{Float64}, size::Int64)::Vector{Float64}
    result = Vector{Float64}(undef, length(ωRange))
    for (it, ω) in enumerate(ωRange)
        result[it] = -imag(greensFunction(ω + iDelta, initialState, diagonal, offDiagonal, size)) / π
    end
    return result
end


end


### Q: Greens function templates?
