include("./modules/tJmodel1D.jl")
include("./modules/operators.jl")
include("./modules/correlations.jl")
include("./modules/utils.jl")

### May be needed on Windows systems due to a bug in libopenblas64_.dll
if Sys.iswindows()
    LinearAlgebra.BLAS.set_num_threads(1)
end

### Run your script
include("./scripts/example_script.jl")

