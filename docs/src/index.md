# tJMagnonHolon

*t-J model in magnon-holon basis.*

A set of modules for t-J model (related) Hamiltonians diagonalization and various spectral functions generation.

!!! note
    Please read through the [Manual](https://docs.julialang.org/en/v1/manual/getting-started/) section of the Julia documentation
    if you're not familiar with Julia programming language.

## Features

- Generate sparse matrix of t-J model Hamiltonian in magnon-holon basis.
- Generate full (small systems) or partial (large systems) set of eigen-values and related eigen-vectors of t-J model.
- Calculate operators action on wave functions in magnon-holon basis.
- Calculate Greens functions and spectral functions of t-J model.
- Save data to JSON files and use them as datasets in Mathematica, R, Matlab, etc.

The [Guide](@ref) provides a tutorial explaining how to use above features.
You will find there practical examples of scripts for running your calculations or how to extend the existing code base with your own operators.

See the [Documentation](@ref) for the complete list of documented modules and features.

In the [Advanced](@ref) section you can find technical discussion of the mathematics behind various features of the code.

## Guide Outline

```@contents
Pages = [
    "guide.md"
]
```

