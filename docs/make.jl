push!(LOAD_PATH, "../src/")

using Documenter

include("../src/modules/tJmodel1D.jl")
include("../src/modules/operators.jl")
include("../src/modules/correlations.jl")
include("../src/modules/utils.jl")

makedocs(
    modules = [Main.tJmodel1D, Main.Operators, Main.Correlations, Main.Utils],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    sitename = "tJMagnonHolon",
    authors = "Piotr Wrzosek",
    pages = [
        "Home" => "index.md",
        "Guide" => "guide.md",
        "Advanced" => "advanced.md",
        "Documentation" => "documentation.md"
    ]
)

deploydocs(
    repo = "github.com/pwrzosek/tJMagnonHolon.git",
)

