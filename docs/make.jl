using ADCSSims
using Documenter

DocMeta.setdocmeta!(ADCSSims, :DocTestSetup, :(using ADCSSims); recursive=true)

makedocs(;
    modules=[ADCSSims],
    authors="Orestis Ousoultzoglou <orousoultzoglou@gmail.com>, Romanos Voulgarakis <voulgarakisromanos@gmail.com>, and contributors",
    repo="https://github.com/EMTech-AUTh/ADCSSims.jl/blob/{commit}{path}#{line}",
    sitename="ADCSSims.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EMTech-AUTh.github.io/ADCSSims.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EMTech-AUTh/ADCSSims.jl",
    devbranch="main",
)
