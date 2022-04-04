using StaticArrays

# Here I reimplemented some functions from StaticArrays. 
# By being less robust, I can usually save some time.

#VECTORS 
#======#

"""
	RealVector(x1::Float64, x2::Float64, x3::Float64)

Static Float64 vector with 3 elements.
"""
const RealVector = SArray{Tuple{3},Float64,1,3}

"""
	VECX

Static cartesian basis vector in the X direction. Equivalent to `RealVector(1.,0.,0.)`
"""
const VECX = RealVector(1., 0., 0.)

"""
	VECY

Static cartesian basis vector in the Y direction. Equivalent to `RealVector(0.,1.,0.)`.
"""
const VECY = RealVector(0., 1., 0.)

"""
	VECZ

Static cartesian basis vector in the Z direction. Equivalent to `RealVector(0.,0.,1.)`.
"""
const VECZ = RealVector(0., 0., 1.)

"""
	VEC0

Static zero vector. Equivalent to `zero(RealVector)`.
"""
const VEC0 = zero(RealVector)

"""
    dot(x::RealVector, y::RealVector)::Float64

Computes the dot product of two vectors.
"""
function dot(x::RealVector, y::RealVector)::Float64
	return x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
end

"""
    norm(x::RealVector)::Float64

Computes the norm of a vector.
"""
function norm(x::RealVector)::Float64
	return sqrt(dot(x,x))
end

#MATRICES 3x3
#===========#

"""
	RealMatrix

Static Float64 matrix with 3x3 elements.
"""
const RealMatrix = SArray{Tuple{3,3},Float64,2,9}


"""
	MAT0

Static 3x3 zero matrix.
"""
const MAT0 = zero(RealMatrix)

"""
	MAT1

Static 3x3 identity matrix.
"""
const MAT1 = RealMatrix(1., 0., 0., 
                        0., 1., 0.,
						0., 0., 1.)

"""
    trace(A::RealMatrix)::Float64

Trace of 3x3 matrix.
"""
function trace(A::RealMatrix)::Float64
    @inbounds return A[1] + A[5] + A[9]
end                        

"""
    dev(A::RealMatrix)::RealMatrix

Deviatoric part of 3x3 matrix. That is, dev(A) = A - Tr(A)*I/3.
"""
function dev(A::RealMatrix)::RealMatrix
    return A - 1/3*trace(A)*MAT1
end

"""
    det(A::RealMatrix)::Float64

Determinant of 3x3 matrix.
"""
function det(A::RealMatrix)::Float64
    @inbounds return (
          + A[1]*A[5]*A[9]
          + A[2]*A[6]*A[7]
          + A[3]*A[4]*A[8]
          - A[7]*A[5]*A[3]
          - A[8]*A[6]*A[1]
          - A[9]*A[4]*A[2]
    )
end

"""
    trans(A::RealMatrix)::RealMatrix

Transpose of 3x3 matrix.
"""
function trans(A::RealMatrix)::RealMatrix
    @inbounds return RealMatrix(A[1],A[4],A[7],A[2],A[5],A[8],A[3],A[6],A[9])
end

"""
    cof(A::RealMatrix)::RealMatrix

Cofactor matrix of 3x3 matrix.
"""
function cof(A::RealMatrix)::RealMatrix
    @inbounds return RealMatrix(
        A[5]*A[9] - A[8]*A[6],
        A[7]*A[6] - A[4]*A[9],
        A[4]*A[8] - A[7]*A[5],
        A[8]*A[3] - A[2]*A[9],
        A[1]*A[9] - A[7]*A[3],
        A[7]*A[2] - A[1]*A[8],
        A[2]*A[6] - A[5]*A[3],
        A[4]*A[3] - A[1]*A[6],
        A[1]*A[5] - A[4]*A[2]
    )
end

"""
    inv(A::RealMatrix)::::RealMatrix

Inverse matrix of 3x3 matrix.
"""
function inv(A::RealMatrix)::RealMatrix
    return 1/det(A)*trans(cof(A))
end

"""
    dot(A::RealMatrix, B::RealMatrix)::Float64

Computes the double dot product of two 3x3 matrices.
"""
function dot(A::RealMatrix, B::RealMatrix)::Float64
	@inbounds return (
        + A[1]*B[1]
        + A[2]*B[2]
        + A[3]*B[3]
        + A[4]*B[4]
        + A[5]*B[5]
        + A[6]*B[6]
        + A[7]*B[7]
        + A[8]*B[8]
        + A[9]*B[9]
    )
end

"""
    norm(x::RealMatrix)::Float64

Computes the norm of a matrix.
"""
function norm(x::RealMatrix)::Float64
	return sqrt(dot(x,x))
end

#MATRICES 2x2
#===========#

"""
	FlatMatrix

Static matrix for faster 2d simulations
"""
const FlatMatrix = SArray{Tuple{2,2},Float64,2,4}

Base.:*(A::FlatMatrix, b::RealVector)::RealVector = begin
    @inbounds return RealVector(A[1]*b[1] + A[3]*b[2], A[2]*b[1] + A[4]*b[2], b[3])
end

"""
	FMAT0

Static 2x2 zero matrix.
"""
const FMAT0 = zero(FlatMatrix)

"""
	FMAT1

Static 2x2 identity matrix.
"""
const FMAT1 = FlatMatrix(1., 0.,
                         0., 1.)

"""
    trace(A::FlatMatrix)::Float64

Trace of 2x2 matrix.
"""
function trace(A::FlatMatrix)::Float64
    @inbounds return A[1] + A[4]
end       

"""
    dev(A::FlatMatrix)::FlatMatrix

Deviatoric part of matrix 2x2. That is, dev(A) = A - Tr(A)*I/2.
"""
function dev(A::FlatMatrix)::FlatMatrix
    return A - 0.5*trace(A)*MAT1
end

"""
    det(A::FlatMatrix)::Float64

Determinant of 2x2 matrix.
"""
function det(A::FlatMatrix)::Float64
    @inbounds return A[1]*A[4] - A[2]*A[3]
end

"""
    trans(A::FlatMatrix)::FlatMatrix

Transpose of 2x2 matrix.
"""
function trans(A::FlatMatrix)::FlatMatrix
    @inbounds return FlatMatrix(A[1],A[2],A[3],A[4])
end

"""
    cof(A::FlatMatrix)::FlatMatrix

Cofactor of 2x2 matrix.
"""
function cof(A::FlatMatrix)::FlatMatrix
    @inbounds return FlatMatrix(A[4], -A[3], -A[2], A[1])
end

"""
    inv(A::FlatMatrix)::FlatMatrix

Inverse of 2x2 matrix.
"""
function inv(A::FlatMatrix)::FlatMatrix
    idet = 1/det(A)
    @inbounds return FlatMatrix(idet*A[4], -idet*A[2], -idet*A[3], idet*A[1])
end

"""
    dot(A::FlatMatrix, B::FlatMatrix)::Float64

Computes the double dot product of two 2x2 matrices.
"""
function dot(A::FlatMatrix, B::FlatMatrix)::Float64
	@inbounds return (
        + A[1]*B[1]
        + A[2]*B[2]
        + A[3]*B[3]
        + A[4]*B[4]
    )
end

"""
    norm(x::FlatMatrix)::Float64

Computes the norm of a matrix.
"""
function norm(x::FlatMatrix)::Float64
	return sqrt(dot(x,x))
end