using Documenter

makedocs(;
    modules=[LaMEM_C],
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
    repo="github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "main",
    devurl = "dev",
    forcepush=true,
    push_preview = true
)

