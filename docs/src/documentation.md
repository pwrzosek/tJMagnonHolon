# Documentation

This page contains a list of implemented types, methods and functions.
All descriptions in the following sections are raw docstrings extracted from corresponding modules.
Sections are organized in logical order possible to read from top to bottom.
Each section ends with alphabetical index organized by feature type in order Type, Method, Function.


## Module tJmodel1D

### Features

```@docs
Main.tJmodel1D.System
```

```@docs
Main.tJmodel1D.System(::Main.tJmodel1D.System)
```

```@docs
Main.tJmodel1D.checkSystem
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
Main.tJmodel1D.getStateInfo
```

```@docs
Main.tJmodel1D.isGreater
```

```@docs
Main.tJmodel1D.bitmov
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
Main.tJmodel1D.LinearCombination
```

```@docs
Main.tJmodel1D.hamiltonian
```

```@docs
Main.tJmodel1D.Model
```

```@docs
Main.tJmodel1D.makeModel
```

```@docs
Main.tJmodel1D.factorize
```
!!! note "Algorithm details"
    For more details about factorization procedure see documentation of 
    [eigsolve](https://jutho.github.io/KrylovKit.jl/stable/man/eig/#KrylovKit.eigsolve) function from [KrylovKit.jl](https://jutho.github.io/KrylovKit.jl/stable) package.

```@docs
Main.tJmodel1D.run
```

### Index
```@index
Pages = ["documentation.md"]
Modules = [Main.tJmodel1D]
```


## Module Operators

### Features

```@docs
Main.Operators.Superposition
```

```@docs
Main.Operators.SystemSuperposition
```

```@docs
Main.Operators.SystemWaveFunction
```

```@docs
Main.Operators.applyOperator
```

### Index
```@index
Pages = ["documentation.md"]
Modules = [Main.Operators]
```


## Module Correlations

### Features

```@docs
Main.Correlations.Krylov
```

```@docs
Main.Correlations.calculateLanczos
```

```@docs
Main.Correlations.greensFunction
```

```@docs
Main.Correlations.run
```

```@docs
Main.Correlations.calculate
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

