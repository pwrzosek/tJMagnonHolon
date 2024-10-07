# Documentation

This page contains a complete documentation of implemented types and functions.
Each module has its own separate section for easier navigation.
Each module section ends with table of contents for that module listed in alphabetical order. 

## Outline
```@contents
Pages = ["documentation.md"]
Depth = 2
```

## Module tJmodel1D

### Features

```@docs
Main.tJmodel1D.System
```

```@docs
Main.tJmodel1D.System(::Main.tJmodel1D.System)
```

```@docs
Main.tJmodel1D.State
```

```@docs
hash
```

```@docs
isequal
```

```@docs
Main.tJmodel1D.isGreater
```

```@docs
Main.tJmodel1D.bitmov
```

```@docs
Main.tJmodel1D.getStateInfo
```

```@docs
Main.tJmodel1D.Basis
```

```@docs
Main.tJmodel1D.makeBasis
```

```@docs
Main.tJmodel1D.sublatticeRotation
```

```@docs
Main.tJmodel1D.Model
```

```@docs
Main.tJmodel1D.makeModel
```

```@docs
Main.tJmodel1D.hamiltonian
```

```@docs
Main.tJmodel1D.LinearCombination
```

```@docs
Main.tJmodel1D.run
```

```@docs
Main.tJmodel1D.checkSystem
```

```@docs
Main.tJmodel1D.factorize
```
!!! note "Algorithm details"
    For more details about diagonalization procedure see documentation of 
    [eigsolve](https://jutho.github.io/KrylovKit.jl/stable/man/eig/#KrylovKit.eigsolve) function from [KrylovKit.jl](https://jutho.github.io/KrylovKit.jl/stable) package.

### Index
```@index
Pages = ["documentation.md"]
Modules = [Main.tJmodel1D]
```


## Module Operators

### Features

```@docs
Main.Operators.applyOperator
```

```@docs
Main.Operators.SystemWaveFunction
```

```@docs
Main.Operators.SystemSuperposition
```

```@docs
Main.Operators.Superposition
```

### Index
```@index
Pages = ["documentation.md"]
Modules = [Main.Operators]
```


## Module Correlations

### Features

```@docs
Main.Correlations.calculate
```

```@docs
Main.Correlations.run
```

```@docs
Main.Correlations.greensFunction
```

```@docs
Main.Correlations.calculateLanczos
```

```@docs
Main.Correlations.Krylov
```

### Index
```@index
Pages = ["documentation.md"]
Modules = [Main.Correlations]
```


## Module Utils

### Features

```@docs
Main.Utils.saveData
```

### Index
```@index
Pages = ["documentation.md"]
Modules = [Main.Utils]
```

