# Guide

The purpose of this [Guide](@ref) is to provide a simple introduction to [tJMagnonHolon](@ref) features.
It will allow you to start producing research relevant data right away.

## Introduction

To start using [tJMagnonHolon](@ref) you need a distribution of Julia.
You can download recent stable version from its official website [julialang.org](https://julialang.org).

Following packages are required for [tJMagnonHolon](@ref) code to work.
```
OrderedCollections
LinearAlgebra
SparseArrays
KrylovKit
DelimitedFiles
Dates
JSON
```
To download them, open terminal and move to the following directory.
```
.../tJMagnonHolon/src/
``` 
In terminal type following command.
```
julia setup.jl
```
Required packages will be automatically downloaded and added to Julia (you need an internet connection). This step has to be done only once. 
You can also add packages by yourself (e.g. if you already have some of them installed and don't wish to upgrade). 

## Tutorial

Let us start with the basics.

### Running the code
 
- open terminal and move to `.../tJMagnonHolon/src/` directory.
    - in terminal type `julia` to start julia process, then in Julia REPL type `include("run.jl")`.
    - alternatively type `julia run.jl` (to run without graphical interface)


File `run.jl` will execute `example_script.jl` from `.../tJMagnonHolon/src/scripts/` directory.
You can add more scripts to `.../tJMagnonHolon/src/scripts/` (or to some other location). 
To execute one or more of those scripts include them in `run.jl` file.
For instance, to run `.../tJMagnonHolon/scripts/my_script.jl` add fillowing line to `run.jl`.
```Julia
include("./scripts/my_script.jl")
```

If you don't want `example_script.jl` to be executed, remove or comment out from `run.jl` following line. 
```Julia
include("./scripts/example_script.jl")
```

### Basic structure

Let us now discuss the features of the [tJMagnonHolon](@ref).
Open `example_script.jl` on a side and compare how the features we are going to discuss are put together in a single script.

!!! tip "[optional] Basic structure"
    It is convenient to define output type for results we want to collect. Let us call it `Data`.
    ```Julia
    Data = OrderedDict{String, Union{Int64, Float64, Array{Float64}, Array{ComplexF64}}}
    ```
    You can write your script in the main scope of Julia but wrapping your script in a function (or few functions) helps in organizing and reusing the code.
    ```Julia
    function calculate()::Data
        ...
    end
    ```
    Above defined function `calculate()` is supposed to return output of type `Data`.
    You can also add some arguments to function `calculate()` or use different function name.


### System definition

To perform any calculations we need to first define parameters of our system.
All the system parameters are handled by [`System`](@ref Main.tJmodel1D.System) structure from `tJmodel1D` module.
Let us use `system` for a corresponding variable name.
There are multiple ways to define `system` but probably most readble is to use 
[`System`](@ref Main.tJmodel1D.System(::Main.tJmodel1D.System)) 
function and provide a set of keyword arguments as in example below. 
The order in which keyword arguments are listed is not important.

```@example
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
```

The first four parameters ``\textcolor{orange}{t}, \textcolor{orange}{J}, \textcolor{orange}{\lambda}, \textcolor{orange}{\alpha}`` are parameters of the Hamiltonian ``\hat{H}`` (in magnon-holon basis).
```math
\hat{H} = \hat{H}_t + \hat{H}_{xy} + \hat{H}_z
```
```math
\hat{H}_t = \textcolor{orange}{t}\sum_{\langle i,j \rangle} \hat{P}_i \hat{h}_i^\dag \hat{h}_j \left( \hat{a}_i + \hat{a}_j^\dag \right) \hat{P}_j + \textrm{H.c.}
```
```math
\hat{H}_{xy} = \frac{\textcolor{orange}{\alpha} \textcolor{orange}{J}}{2} \sum_{\langle i,j \rangle} \hat{h}_i \hat{h}_i^\dag \left(\hat{P}_i \hat{P}_j \hat{a}_i \hat{a}_j + \mathrm{H.c} \right) \hat{h}_j \hat{h}_j^\dag
```
```math
\hat{H}_{z} = \frac{\textcolor{orange}{J}}{2} \sum_{\langle i,j \rangle} \hat{h}_i \hat{h}_i^\dag \left(\hat{a}_i^\dag \hat{a}_i + \hat{a}_j^\dag \hat{a}_j - 2 \textcolor{orange}{\lambda} \hat{a}_i^\dag \hat{a}_i \hat{a}_j^\dag \hat{a}_j - 1 \right) \hat{h}_j \hat{h}_j^\dag
```

The last four parameters ``\mathrm{\textcolor{orange}{size}, \textcolor{orange}{electrons}, \textcolor{orange}{spinsUp}, \textcolor{orange}{momentum}}`` are conserved by the Hamiltonian and point to a single subspace of the model (orthogonal to other subspaces).
See section [Advanced] for detailed discussion.

#### Other ways to define System

Used in the above example values for system parameters are equal to the default values defined inside the `tJmodel1D` module. 
To create a `system` with default values it is enaugh to call 
[`System`](@ref Main.tJmodel1D.System(::Main.tJmodel1D.System)) 
function without arguments.
```@example
system = Main.tJmodel1D.System()
```

You can also provide any subset of keyword arguments. Not provided arguments will take default values.
```@example
system = Main.tJmodel1D.System(spinsUp = 0)
```

Additionally, if you are alergic to explicitely writing arguments names, you can use a constructor of [`System`](@ref Main.tJmodel1D.System) structure. 
In this case the order of arguments matters and it follows: `(t, J, λ, α, size, electrons, spinsUp, momentum)`. 
```@example 1
system = Main.tJmodel1D.System(1.0, 1.0, 1.0, 1.0, 16, 16, 8, 0)
```

[`System`](@ref Main.tJmodel1D.System) structure is immutable meaning that once you define it you cannot change values of its fields.
But sometimes you may want to create a new instance of [`System`](@ref Main.tJmodel1D.System) 
(let's call it `newSystem`) that has only one or few fields changed with respect to some previously defined `system`.
In such case you can provide `system` as an argument to [`System`](@ref Main.tJmodel1D.System(::Main.tJmodel1D.System)) function, see below.
```@example 1
newSystem = Main.tJmodel1D.System(system, electrons = system.electrons - 1, spinsUp = system.spinsUp - 1)
```
Defined above `newSystem` has one electron less and one spin up less than `system` (e.g. one electron with spin up was removed). Values of other parameters are copied from `system`.

### Model generation and factorization

Once `system` is assigned you can call defined in `tJmodel1D` module function [`run`](@ref Main.tJmodel1D.run) to generate Hamiltonian matrix for subspace described by `system` and factorize it.
```Julia
system, basis, model, factorization = Main.tJmodel1D.run(system)
```
This function returns 3 objects and a tuple:
- [`System`](@ref Main.tJmodel1D.System) - parameters of the system subspace.
- [`Basis`](@ref Main.tJmodel1D.Basis) - basis of representative magnon-holon states for the system subspace.
- [`Model`](@ref Main.tJmodel1D.Model) - generated sparse matrix of the Hamiltonian for the system subspace.
- `factorization = (eigenvalues, eigenvectors, convergenceInfo)` - results of diagonalization procedure.

!!! note "Type and value note"
    - Note that `eltype(eigenvalues) <: ComplexF64`. Complex numbers may be returned if `Model` matrix has complex coefficients (even if imaginary parts are numerically zero). It is safe to assume that:
    ```Julia        
    eigenvalues::Vector{ComplexF64}
    eigenvectors::Vector{Vector{ComplexF64}}
    ```
    - Note that eigenvectors are determined up to a random (complex) phase that changes from run to run.

You can easily access calculated eigenvalues and eigenvectors. For example, the first value and its corresponding vector:
```Julia
ψ = eigenvectors[1]
E0 = real(eigenvalues[1])
```

!!! note
    Eigenvalues are sorted with respect to real part from smallest to largest.

By default only 1 eigenvalue with smallest real part will be calculated. To calculate more eigenvalues use `howmany` keyword argument.
For example, code below finds singlet and triplet energy of 2-site antiferromagnetic spin 1/2 Heisenberg chain for spin coupling J = 1.
```@example 2
system = Main.tJmodel1D.System(J = 1.0, size = 2, electrons = 2, spinsUp = 1)
system, basis, model, factorization = Main.tJmodel1D.run(system, howmany = 13) 
                            ### we pretend we don't know it should be 2 ---^
eigenvalues, eigenvectors, convergenceInfo = factorization
real.(eigenvalues)
[-2.0, 0.0] # hide
```

!!! warn "Convergence"
    To make sure that there are no poorly converged values, it is always good to check `convergenceInfo`.
    ```Julia
    println(convergenceInfo)
    ```
    Norms of residuals close to zero indicate well converged values. 

If you calculate only few eigenvalues, it is unlikely to happen, but in case you cannot converge desired number of eigenvalues, you may need to increase dimension of Krylov subspace.
You can achieve this by setting keyword argument `kryldim` to value higher than 30 (which is default value).
You can also set value smaller than 30 to save memory e.g. when looking for the lowest eigenvalue (and eigenvector) of a large system.
Remember that `kryldim` cannot be smaller than `howmany`.
```Julia
system, basis, model, factorization = Main.tJmodel1D.run(system, howmany = 40, kryldim = 100)
```

#### Skipping diagonalization

If you just need Hamiltonian matrix and its basis tell the [`run`](@ref Main.tJmodel1D.run) function to skip the factorization procedure by setting keyword argument `eigsolve = false`.
In such case `factorization` will recieve `missing` value from [`run`](@ref Main.tJmodel1D.run).
```Julia
system, basis, model, factorization = Main.tJmodel1D.run(system, eigsolve = false)
```

!!! tip
    Use dummy variables if you don't need some of the returned objects.
    ```Julia
    _, _, _, factorization = Main.tJmodel1D.run(system)
    system, basis, model, _ = Main.tJmodel1D.run(system, eigsolve = false)
    ```

### List of operators

Operators are defined in module `Operators`. You can use any of the predefiend operators listed below.
You can also add your own operators to `Operators` module.
See section [Operators](@ref) to learn about operator functions design.

- ``\hat{S}_{k}^{z}``, ``\hat{S}_{r}^{z}``, ``\hat{S}_{k}^{+}``, ``\hat{S}_{r}^{+}``, ``\hat{S}_{k}^{-}``, ``\hat{S}_{r}^{-}``
```Julia
Sk_z(k::Int64; state::State, system::System)
Sr_z(r::Int64; state::State, system::System)
Sk_plus(k::Int64; state::State, system::System)
Sr_plus(r::Int64; state::State, system::System)
Sk_minus(k::Int64; state::State, system::System)
Sr_minus(r::Int64; state::State, system::System)
```

- ``\hat{\tilde{c}}_{k\uparrow}``, ``\hat{\tilde{c}}_{r\uparrow}``, ``\hat{\tilde{c}}_{k\downarrow}``, ``\hat{\tilde{c}}_{r\downarrow}``
```Julia
ck_up(k::Int64; state::State, system::System)
cr_up(r::Int64; state::State, system::System)
ck_down(k::Int64; state::State, system::System)
cr_down(r::Int64; state::State, system::System)
```

- ``\hat{\tilde{c}}_{k\uparrow}^{\dag}``, ``\hat{\tilde{c}}_{r\uparrow}^{\dag}``, ``\hat{\tilde{c}}_{k\downarrow}^{\dag}``, ``\hat{\tilde{c}}_{r\downarrow}^{\dag}``
```Julia
ck_up_dag(k::Int64; state::State, system::System)
cr_up_dag(r::Int64; state::State, system::System)
ck_down_dag(k::Int64; state::State, system::System)
cr_down_dag(r::Int64; state::State, system::System)
```

### Spectral function and Greens function



To generate the spectral function, define your resolution parameters,
- artificial broadening ``\delta`` of the peaks,
- set of ``\omega`` points at which spectral function should be calculated.
For example:
```Julia
δ = 0.02
ωRange = collect(-3:0.002:7)
```
Smaller values of ``\delta`` make peaks sharper. But to actually see the effect you need to set small enough step in ``\omega`` to resove it.
Otherwise there will be too few points per peak to properly cover its shape. Step ``\delta / 5`` is usually small enough.

If the operator you use takes arguments (i.e. it has some indices), define a set of arguments to iterate over. 
For example, if you use ``\hat{S}_{k}^{+}``, define range of momenta ``k`` you want to evaluate.
```Julia
kRange = collect(0:system.size)
```
The actual values depend on how the operators arguments are defined. Each ``k`` in the above set corresponds to momentum `2πk / system.size`. See section [Operators](@ref) for more details.

We can now specify dimensions for `spectrum` variable where we are going to store the specral function results. We fill it with zeros to make sure there are no undefined values in it.
```Julia
spectrum = zeros(Float64, length(ωRange), length(kRange))
```

Now we can calculate the spectral function. For that we call function [`run`](@ref Main.SpectralFunction.run) from `SpectralFunction` module and supply it with necessary arguments.
```Julia
for k in kRange
    spectrum[:, k + 1] = Main.SpectralFunction.run(ωRange .+ E0, δ, system, ψ, operator, k)
end
```
Above we shifted the set of energy points `ωRange` by the energy `E0` of the bare wave function `ψ` (this is optional, but usually that's what you want to do).
The `operator` with argument `k` is applied to `ψ` internally. Proper interpratation of `ψ` (which is just a vector of complex numbers) is allowed by `system` argument.

!!! tip
    On HPC and for large systems, instead of a single loop, consider separate runs for different operator arguments (or subsets of those).
    For example, you can run parallel calculations on separate nodes to obtain your results faster.


## Operators

### Applying operators to arbitrary wave functions

- where to put new operators
- operator function design
- send to advanced for maths

