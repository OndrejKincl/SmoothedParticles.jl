module SmoothedParticles
import WriteVTK
import SparseArrays
import StaticArrays
import Match

include("algebra.jl")
export RealVector, VECX, VECY, VECZ, VEC0,
       RealMatrix, MAT0, MAT1,
       FlatMatrix, FMAT0, FMAT1,
       norm, dot

include("structs.jl")
export AbstractParticle, 
       ParticleSystem, 
       ParticleField, 
       DataField

include("kernels.jl")
export wendland2, Dwendland2, rDwendland2, 
       wendland3, Dwendland3, rDwendland3, DDwendland3,
       spline24, Dspline24, rDspline24,
       spline23, Dspline23, rDspline23

include("core.jl")
export apply!, 
       apply_unary!, 
       apply_binary!, 
       create_cell_list!, 
       dist, 
       assemble_matrix, 
       assemble_vector

include("grids.jl")
export Grid,
       Squaregrid,
       Hexagrid,
       CubicGrid

include("geometry.jl")
export generate_particles!, 
       Shape, 
       Rectangle, 
       Circle, 
       Ellipse, 
       Ball,
       Box,
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
