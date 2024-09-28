# Documentation

This page contains a complete documentation of implemented types and functions.
Each module has its own separate section for easier navigation.
Each module section ends with table of contents for that module listed in alphabetical order. 

## Navigation
```@contents
Pages = ["documentation.md"]
Depth = 3
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


## Module SpectralFunction

### Features

```@docs
Main.SpectralFunction.run
```

!!! note "Returned type note"
    Returned type of `Main.SpectralFunction.run` is `Vector{Float64}` for spectral function and `Vector{Complex{Float64}}` for Greens function calculation.

```@docs
Main.SpectralFunction.calculate
```

```@docs
Main.SpectralFunction.spectralFunction
```

```@docs
Main.SpectralFunction.greensFunction
```

```@docs
Main.SpectralFunction.calculateLanczos
```

```@docs
Main.SpectralFunction.Krylov
```

### Index
```@index
Pages = ["documentation.md"]
Modules = [Main.SpectralFunction]
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

