using Documenter

makedocs(;
    authors="Anton Popov, Boris Kaus",
    sitename="LaMEM",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "man/Home.md",
        "User Guide" => Any[
            "Installation" =>  "man/Installation.md",
            "Getting Started" =>  "man/InitialModelSetup.md",
            "Examples" =>  "man/Examples.md",
            "Features" =>  "man/Features.md",
            "LaMEM Development" => "man/LaMEM_Development.md",
            "LaMEM Debugging" => "man/Debugging.md",
        ],
    ],  
)

deploydocs(;
    repo="github.com/UniMainzGeo/LaMEM.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "master",
    devurl = "dev",
    forcepush=true,
    push_preview = true
)

