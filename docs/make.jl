push!(LOAD_PATH, realpath("../src"))

using Documenter
using Grid

makedocs(
    modules = [Grid],
    clean = true,
    sitename = "Julia Drift-Diffusion",
    pages = Any[
        "Home" => ["index.md"]
        "Grid" => ["SimplexGrid.md", "FiniteVolumeGrid.md"]
    ],
    checkdocs = :all,
    format = :html
)

