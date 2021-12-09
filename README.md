# SPHLib.jl
A lightweight framework for smoothed particle hydrodynamics (SPH). Requires Julia 1.5 or newer.
This is still in development and the API can be subject to changes. Currently, only 2d is supported. 

The tools SPHLib.jl provides are:
+ family of SPH kernels
+ neighbor list generation
+ multithreaded local particle-particle interactions
+ particle-particle sparse matrix assembly
+ particle initialization
+ I/O of simulation to .pvd file, which can be opened in [paraview](https://www.paraview.org/)

For a quick start, you can try example "collapse_dry.jl". To this end, follow these steps:

1: Clone or download whole project into your desktop

2: In the directory "[SPHLib.jl folder]/SPHLib.jl/examples" use command (replacing the number of threads with whatever you want and your computer can handle) 
```bash
    $ julia --threads 4
```
    
3: In julia terminal, type 
```@repl
    include("collapse_dry.jl")
```

4: If you get an error that a Package XY is missing, download it using Pkg and goto step 3

5: In julia terminal, type 
```@repl
    collapse_dry.main()
```

6: Wait for the simulation to end

7: Once the simulation ends, it will create an output in "[SPHLib.jl folder]/SPHLib.jl/examples/results/collapse_dry/result.pvd". Open the file using paraview client