using Documenter

makedocs(;
    authors="Anton Popov, Boris Kaus",
    sitename="LaMEM",
    root="docs",
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
        "Development" => Any[
            "LaMEM Development" => "man/LaMEM_Development.md",
            "LaMEM Debugging" => "man/Debugging.md"
        ],
        "Code of conduct" => "man/CODE_OF_CONDUCT.md"
    ],  
)

deploydocs(;
    root = "docs",
    repo="github.com/UniMainzGeo/LaMEM.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "master",
    devurl = "dev",
    versions="v#",
    forcepush=true,
    push_preview = false
)
