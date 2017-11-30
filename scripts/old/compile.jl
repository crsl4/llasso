## julia script to compile Atom jmd files
## run in terminal as: julia compile.jl name for a file called name.jmd
## Claudia September 2017

isempty(ARGS) && error("need to include filename as argument")

using Weave
weave(string(ARGS[1],".jmd"),doctype="pandoc2pdf", fig_ext=".pdf")

