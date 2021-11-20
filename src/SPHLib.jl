module SPHLib
import WriteVTK
import SparseArrays
import StaticArrays
import Match

include("structs.jl")
export AbstractParticle, 
       ParticleSystem, 
       ParticleField, 
       DataField,
       Vec2,Mat2

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
       assemble_vector,
       dot,
       norm

include("grids.jl")
export Grid,
       Squaregrid,
       Hexagrid

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
       is_inside

include("IO.jl")
export save_frame!, 
       new_pvd_file,
       save_pvd_file

end