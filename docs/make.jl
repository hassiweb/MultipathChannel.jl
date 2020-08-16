using MultipathChannel
using Documenter

makedocs(;
    modules=[MultipathChannel],
    authors="hassiweb <nryk.hashimoto@gmail.com>",
    repo="https://github.com/hassiweb/MultipathChannel.jl/blob/{commit}{path}#L{line}",
    sitename="MultipathChannel.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hassiweb.github.io/MultipathChannel.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hassiweb/MultipathChannel.jl",
)
