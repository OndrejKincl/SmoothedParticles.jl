module SPHLib
import WriteVTK
import SparseArrays
import Match

include("structs.jl")
export AbstractParticle, 
       ParticleSystem, 
       ScalarField, 
       VectorField, 
       DataField, 
       @define_particle

include("kernels.jl")
export wendland2, 
       Dwendland2, 
       rDwendland2, 
       spline24, 
       Dspline24, 
       rDspline24,
       spline23, 
       Dspline23, 
       rDspline23

include("core.jl")
export apply!, 
       apply_unary!, 
       apply_binary!, 
       create_cell_list!, 
       dist, 
       assemble_matrix, 
       assemble_vector

include("geometry.jl")
export generate_particles!, 
       Shape, 
       Rectangle, 
       Circle, 
       Ellipse, 
       Polygon, 
       BooleanUnion, 
       BooleanIntersection, 
       BooleanDifference, 
       Specification, 
       BoundaryLayer, 
       is_inside, 
       snap_to_grid

include("IO.jl")
export save_frame!, 
       new_pvd_file,
       save_pvd_file, 
       capture_frame 
       #write_vtk, 
       #read_vtk!

end
