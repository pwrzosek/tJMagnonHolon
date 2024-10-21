# Guide

The purpose of this [Guide](@ref) is to provide a simple introduction to [tJMagnonHolon](@ref) features.
It will allow you to start producing research relevant data right away.

```@contents
Pages = [
    "guide.md"
]
Depth = 3
```

## Introduction

To start using [tJMagnonHolon](@ref) you need a distribution of Julia.
You can download recent stable version from its official website [julialang.org](https://julialang.org).

!!! tip
    Check out the [Manual](https://docs.julialang.org/en/v1/manual/getting-started/) section of the Julia documentation
    if you're not familiar with Julia programming language.

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
 
Open terminal and move to the following directory.
```
.../tJMagnonHolon/src/
``` 
- in terminal type `julia` to start julia process, then in Julia REPL type `include("run.jl")`.
- alternatively type `julia run.jl` (to run without graphical interface).

This will execute `example_script.jl` from `.../tJMagnonHolon/src/scripts/` directory with all the [tJMagnonHolon](@ref) modules included.
You can add more scripts to `.../tJMagnonHolon/src/scripts/` (or to any location of your choice). 
To execute one or more of those scripts include them in `run.jl` file.
```Julia
include("./scripts/my_script.jl")
```

---

Let us now learn about the features of the [tJMagnonHolon](@ref).

!!! tip
    Open `example_script.jl` on a side and compare how the discussed features are put together into a single script.

!!! note "[optional] Basic structure"
    It is convenient to define output type for results we want to collect. Let us call it `Data`.
    ```Julia
    Data = OrderedDict{String, Union{Int64, Float64, Array{Float64}, Array{ComplexF64}}}
    ```
    You can write your script in the main scope of Julia but wrapping it in a function (or few functions) helps in organizing and reusing the code.
    ```Julia
    function script()::Data
        ...
    end
    ```
    Above defined function `script()` is supposed to return output of type `Data`.
    You can also add some arguments to function `script()` or use different function name.


### System definition

To perform any calculations we need to first define parameters of our system.
All the system parameters are handled by [`System`](@ref Main.tJmodel1D.System) structure from `tJmodel1D` module.
Let us use `system` for a corresponding variable name.
There are multiple ways to define `system` but probably most readable is to use 
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

##### Hamiltonian

The first four parameters ``\textcolor{orange}{t}, \textcolor{orange}{J}, \textcolor{orange}{\lambda}, \textcolor{orange}{\alpha}`` are parameters of the Hamiltonian ``\hat{H}`` (in magnon-holon basis).
```math
\hat{H} = \hat{H}_t + \hat{H}_{xy} + \hat{H}_z
```
```math
\hat{H}_t = \textcolor{orange}{t}\sum_{\langle i,j \rangle} \hat{h}_i^\dag \hat{h}_j \left( \hat{a}_i + \hat{a}_j^\dag \right) + \textrm{H.c.}
```
```math
\hat{H}_{xy} = \frac{\textcolor{orange}{\alpha} \textcolor{orange}{J}}{2} \sum_{\langle i,j \rangle} \hat{h}_i \hat{h}_i^\dag \left( \hat{a}_i \hat{a}_j + \hat{a}_i^\dag \hat{a}_j^\dag \right) \hat{h}_j \hat{h}_j^\dag
```
```math
\hat{H}_{z} = \frac{\textcolor{orange}{J}}{2} \sum_{\langle i,j \rangle} \hat{h}_i \hat{h}_i^\dag \left(\hat{a}_i^\dag \hat{a}_i + \hat{a}_j^\dag \hat{a}_j - 2 \textcolor{orange}{\lambda} \hat{a}_i^\dag \hat{a}_i \hat{a}_j^\dag \hat{a}_j - 1 \right) \hat{h}_j \hat{h}_j^\dag
```

The last four parameters ``\mathrm{\textcolor{orange}{size}, \textcolor{orange}{electrons}, \textcolor{orange}{spinsUp}, \textcolor{orange}{momentum}}`` are conserved by the Hamiltonian and point to a single subspace of the model (orthogonal to other subspaces).
See section [Advanced](@ref) for detailed discussion.

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

Once `system` is known, you can call defined in `tJmodel1D` module function [`run`](@ref Main.tJmodel1D.run) to generate Hamiltonian matrix for subspace described by `system` and factorize it.
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
eigenvalues, eigenvectors, info = factorization
ψ = eigenvectors[1]
E0 = real(eigenvalues[1])
```

!!! note
    Eigenvalues are sorted with respect to real part from smallest to largest.

By default only 1 eigenvalue with smallest real part will be calculated. To calculate more eigenvalues use `howmany` keyword argument.
For example, code below finds singlet and triplet energy of 2-site antiferromagnetic spin 1/2 Heisenberg chain (periodic) for spin coupling J = 1.
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

There should be no problems with convergence if you calculate just few eigenvalues. But in case you cannot converge desired number of eigenvalues, you may need to increase dimension of Krylov subspace.
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

### Operators

Operators are defined in module `Operators`. You can use any of the predefiend operators listed below.

#### List of operators

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

If you're interested in calculations of correlation functions (Greens / Spectral function), 
you can skip further subsections of [Operators](@ref) section and move straight to section [Correlation functions](@ref).

#### Applying operators to wave-functions

For given wave-function ``\psi`` from `system` subspace and operator ``\hat{O}_I`` (where ``I`` - ordered collection of indices)
the result of ``\hat{O}_I \vert \psi \rangle`` can be calculated with [`applyOperator`](@ref Main.Operators.applyOperator) function.
For example, this is how to remove an electron with spin up and with momentum ``k = 0``,
```Julia
operator = Main.Operators.ck_up
result = Main.Operators.applyOperator(system, ψ, operator, 0)
```
The first 3 arguments are always `system`, wave function `ψ` and `operator`. Futher arguments are considered `operator` indices and will be passed to the `operator` function. 
It is important to remember that the above `result` may in general overlap with many different subspaces of the model.
The `result` is a hash table with keys of type [`System`](@ref Main.tJmodel1D.System) describing the resulting subspaces and values of type `Vector{ComplexF64}` representing corresponding subspace wave-functions.
```@example
result = Main.Operators.SystemWaveFunction(Main.tJmodel1D.System() => [1.0im, 0.0]) # hide
typeof(result) 
```
Once an operator is applied, one can for example calculate its expectation value on `ψ`.
```Julia
E = if haskey(result, system)
    dot(ψ, result[system])
else
    0
end
```

#### Custom Operators

You can extend the original code base with your own operators. All you have to do is to put your operator function 
in `custom.jl` file located in `.../tJMagnonHolon/src/modules/mods/` directory. 
Every operator function follows the same design.
```Julia
function operator_name(args...; state::State, system::System)::SystemSuperposition
    ...
end
```
It takes any number of arguments `args...` which may for example represent operator indices. The magnon-holon basis `state` and corresponding `system` are passed as keyword arguments.
For each basis state, the number of terms produced by the action of the operator is proportional to the `system.size` rather than `length(basis)`. For this reason,
the output type is a sparse representation of a wave-function that may overlap with more than one subspace of the model.
```Julia
SystemSuperposition = OrderedDict{System, Superposition}
Superposition = OrderedDict{State, ComplexF64} # state is paired with its coefficient
```

If you plan to add a custom operator, it may be useful to have a look at the operators derivations in the [Advanced](@ref) section. 

### Correlation functions

You can calculate Greens/correlation functions for the [Hamiltonian](@ref) ``\hat{H}`` with [tJMagnonHolon](@ref).
The following formula shows what kind of expression can be evaluated.
```math
\langle \psi \vert \hat{O}_{I}^{\dag} \frac{1}{\omega - \hat{H} + i\delta} \hat{O}_{I} \vert \psi \rangle
```
Above, ``\psi`` is any wave-function defined for a certain `system` subspace. Operator ``\hat{O}_{I}`` is any operator expressable in magnon-holon basis (i.e. with holon and magnon creation and annihilation operators).
This operator can depend on any arbitrary ordered collection of indices ``I``. For example ``I=(k,q)`` might represent momenta of two particles introduced to the system. 

To generate a correlation function, define your resolution parameters,
- artificial broadening ``\delta`` of the spectral features,
- set of ``\omega`` points at which spectral features should be calculated.
For example:
```Julia
δ = 0.02
ωRange = collect(-3:0.002:7)
```
Smaller values of ``\delta`` will make features sharper and thinner. Accordingly, you need to set small enough step in ``\omega`` to resolve them.
Otherwise there will be too few points per peak to properly cover its shape. Step ``\delta / 10`` is usually small enough.

If the operator you use takes arguments (i.e. it has some indices), define a set of arguments to iterate over. 
For example, if you use ``\hat{S}_{k}^{+}``, define range of momenta ``k`` you want to evaluate.
```Julia
kRange = collect(0:system.size)
```
The actual values of operators arguments depend on their definition in the code. Here each ``k`` in the above set corresponds to momentum `2πk / system.size`.

We can now specify dimensions for `correlations` variable where we are going to store the results. We fill it with zeros to make sure there are no undefined values in it.
```Julia
correlations = zeros(ComplexF64, length(ωRange), length(kRange))
```

Now we can calculate the correlation function. For that we call function [`calculate`](@ref Main.Correlations.calculate) from `Correlations` module and supply it with necessary arguments.
```Julia
for k in kRange
    correlations[:, k + 1] = Main.Correlations.calculate(ωRange .+ E0, δ, system, ψ, operator, k)
end
```
Above we shifted the set of energy points `ωRange` by the energy `E0` of the bare wave-function `ψ` for better alignment (the shift is optional).
The `operator` with argument `k` is applied to `ψ` internally. Proper interpretation of `ψ` (which is just a vector of complex numbers) is allowed by `system` argument.

!!! tip
    On HPC and for large systems, instead of a single loop, consider separate runs for different operator arguments (or subsets of those).
    For example, you can run parallel calculations on separate nodes to obtain your results faster.

Apart from spectrum resolution settings, one more parameter has an influence on the quality of generated results. It is the maximum depth of recursion in the used algorithm for
Greens function generation. You can change this parameter by setting `kryldim` keyword argument. The default value is `kryldim = 400`.
In general, values between 200 and 500 should be optimal for most calculations. You can set smaller values to speed up calculations for fast lookup, but the result will lose some of its details.
Use larger `kryldim > 500` only if you need to zoom in on a small ``\omega`` window with relatively small broadening ``\delta``. 
```Julia
Main.Correlations.calculate(ωRange, δ, system, ψ, operator, operatorArgs..., kryldim = 100)
```

