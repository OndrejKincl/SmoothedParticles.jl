# Kernels

Smoothing kernels are the guts of SPH. They measure the strength of interaction
between neighbouring particles based on their distance. They are also used to
interpolate particle variables into continuous Eulerian fields.

```@autodocs
Modules = [SPHLib]
Pages = ["kernels.jl"]
```
