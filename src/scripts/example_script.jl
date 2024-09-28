### --- Example script for spectral function calculation --- ### 

using OrderedCollections

SpectrumData = OrderedDict{String, Union{Int64, Float64, Array{Float64}, Array{ComplexF64}}}

function calculate()::SpectrumData
# 1. Define system    
    system = Main.tJmodel1D.System(
        t           = 1.0,      # hole hopping
        J           = 1.0,      # spin coupling
        λ           = 1.0,      # magnon interaction
        α           = 1.0,      # XXZ anisotropy scaling 
        size        = 16,       # number of lattice sites
        electrons   = 16,       # number of electrons
        spinsUp     = 8,        # number of spins up
        momentum    = 0         # internal momentum subspace
    )
    
# 2. Define bare wave function | ψ > and energy shift E0 (e.g. ground state of above defiend system)
    _, _, _, factorization = Main.tJmodel1D.run(system)
    eigenvalues, eigenvectors, info = factorization
    
    ψ   = eigenvectors[1]       # GS wave function
    E0  = real(eigenvalues[1])  # GS energy
    
    println("\n", info) # convergence info (good to check)
 
# 3. Define operator
    operator = Main.Operators.Sk_plus
    
# 4. Spectrum generation
    iδ = 0.02im
    ωRange = collect(-3:0.002:7)
    kRange = collect(0:system.size)
    
    spectrum = zeros(Float64, length(ωRange), length(kRange))

    # iterate over operator indices/args...
    for k in kRange
        spectrum[:, k + 1] = Main.SpectralFunction.run(ωRange .+ E0, iδ, system, ψ, operator, k)
    end

# 5. Return whatever you need
    return SpectrumData(
        "coupling"      => system.J,
        "interaction"   => system.λ,
        "size"          => system.size,
        "momentum"      => (2 / system.size) .* kRange, # momenta in units of π
        "energy"        => ωRange,
        "spectrum"      => spectrum
    )
end


### Run calculation (script above)
spectrumData = calculate()


### Saving Spectrum Data
Main.Utils.saveData([spectrumData], name = "test")


### Plotting [optional; only for quick lookup]
using Plots

k = spectrumData["momentum"]
ω = spectrumData["energy"]
A = spectrumData["spectrum"]

mapfig = heatmap(k, ω, A, clim = (0, 1), xlabel = "k / π", ylabel = "ω");
display(mapfig) # image window will disappear if julia process stops

