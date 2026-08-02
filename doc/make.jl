using Documenter


pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add LaMEM to environment stack
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
Pkg.activate(@__DIR__)
Pkg.instantiate()

using LaMEM_C

makedocs(;
    authors="Anton Popov, Boris Kaus",
    sitename="LaMEM",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "man/Home.md",
        "Quick start" => "man/Quickstart.md",
        "User Guide" => Any[
            "Installation" =>  "man/Installation.md",
            "Getting Started" =>  "man/GettingStarted.md",
            "Initial Model setup" =>  "man/InitialModelSetup.md",
        ],
        "Release Notes" => Any[
            "Upgrading from v2.2.1 to v3.0.0" => "man/Upgrade_v2.2.1_to_v3.0.0.md",
        ],
        "Development" => Any[
            "LaMEM Development" => "man/LaMEM_Development.md",
            "LaMEM Debugging" => "man/Debugging.md"
        ],
        "Code of conduct" => "man/CODE_OF_CONDUCT.md"
    ],  
)

deploydocs(;
    repo="github.com/UniMainzGeo/LaMEM.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "master",
    devurl = "dev",
    versions="v#",
    forcepush=true,
    push_preview = false
)
