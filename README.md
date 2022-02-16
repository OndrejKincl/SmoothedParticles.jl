# SmoothedParticles.jl
A lightweight framework for smoothed particle hydrodynamics (SPH) for 2d and 3d. Requires Julia 1.5 or newer.
This is still in development and the API can be subject to changes.

[Documentation](https://ondrejkincl.github.io/SmoothedParticles.jl/dev/index.html)



# Quick start
  
This package can be installed from Julia terminal using commands:
```
import Pkg
Pkg.add("SmoothedParticles")
```



For start, you can try to run one of the examples. Follow these steps:

1: Clone, download or copy file "collapse_dry.jl" from the [examples](https://github.com/OndrejKincl/SmoothedParticles.jl/tree/master/examples). 

2: Open Julia from the destination folder with command:

```
julia -t N
```

replacing `N` with number of cores that you wish to use in your simulation. 

3: Write:

```
include("collapse_dry.jl")
```

4: If this fails (of course it does) because package XY is missing, download it using 

```
import Pkg
Pkg.add("XY")
```

and repeat step 3.

5: Type

```
collapse_dry.main()
```

to run your simulation.

6: Wait for the simulation to end. Note that this can take several minutes. You should see info about time frames printed to the console.

7: Once it finishes, it creates a new file "results/collapse_dry/result.pvd" in the folder where the example was downloaded. Open it in [paraview](https://www.paraview.org/) to display the result. 
