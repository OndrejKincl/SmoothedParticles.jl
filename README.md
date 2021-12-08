[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OndrejKincl.github.io/SPHLib.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://OndrejKincl.github.io/SPHLib.jl/dev)

# SPHLib.jl
A Julia framework for smoothed particle hydrodynamics (SPH) with a focus on (in this order):
1) robustness
2) performance
3) easy usage

Requires Julia 1.5 or newer.
This is still in development and the framework can be subject to changes.

For a quick start, you can try example "collapse_dry.jl". To this end, follow these steps
1) download this project into your desktop
2) open julia in the directory SPHLib.jl/examples (to use e.g. 8 threads open julia with command "julia -t 8")
3) in julia terminal, type "include("collapse_dry.jl")" and then "collapse_dry.main()".
4) you can be prompted to download some packages, so do that
5) once the simulation ends, this will create an output SPHLib.jl/examples/results/collapse_dry/result.pvd
6) open the file in [paraview](https://www.paraview.org/)
