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
        "User Guide" => Any[
            "1. Installation" =>  "man/Installation.md",
            "2. Getting Started" =>  "man/GettingStarted.md",
            "3. Initial Model setup" =>  "man/InitialModelSetup.md",
            "4. Examples" =>  "man/Examples.md",
            "5. Features" =>  "man/Features.md",
            "6. LaMEM Development" => "man/LaMEM_Development.md",
            "7. LaMEM Debugging" => "man/Debugging.md",
        ],
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
