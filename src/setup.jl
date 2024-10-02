using Pkg

pkglist = [
    "OrderedCollections",
    "LinearAlgebra",
    "SparseArrays",
    "KrylovKit",
    "DelimitedFiles",
    "Dates",
    "JSON"
]

for p in pkglist
    Pkg.add(p)
end

