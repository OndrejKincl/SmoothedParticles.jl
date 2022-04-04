using WriteVTK
#using VTKDataIO

"""
    DataStorage

Storage for paraview input/output.
"""
mutable struct DataStorage
    file::WriteVTK.CollectionFile
    path::String
    frame::Int64
end

"""
    new_pvd_file(path::String)::DataStorage

Cteates new `DataStorage` in a given `path`.
"""
function new_pvd_file(path::String)::DataStorage
    if !ispath(path)
        mkpath(path)
        @info("created new path: "*path)
    end
    return DataStorage(paraview_collection(path*"/result"), path, 0)
end

"""
    save_pvd_file(data::DataStorage)

Saves and closes `DataStorage`.
"""
function save_pvd_file(data::DataStorage)
    vtk_save(data.file)
end

function capture_frame(sys::ParticleSystem, data::DataStorage)
    N = length(sys.particles)
    points = zeros(Float64, 3, N)
    for k in 1:N, i in 1:3
        points[i,k] = sys.particles[k].x[i]
    end
    cells = [MeshCell(PolyData.Verts(), [i]) for i in 1:N]
    vtk_file = vtk_grid(data.path*"/frame"*string(data.frame), points, cells)
    return vtk_file
end

"""
    save_frame!(data::DataStorage, sys::ParticleSystem, vars::DataField...)

Inserts one time frame into a `DataStorage` that includes all `vars...` as fields.
"""
function save_frame!(data::DataStorage, sys::ParticleSystem, vars::Symbol...)
    frame = capture_frame(sys, data)
    for var in vars
        field = ParticleField(sys, var)
        Type = attribute_type(get_particle_type(sys), var)
        N = length(sys.particles)
        if Type <: Number
            frame[string(var)] = field
        elseif Type <: AbstractArray
            vals = zeros(size(Type)..., N)
            for k in 1:N
                for i in CartesianIndices(field[k])
                    vals[i,k] = field[k][i]
                end
            end
            frame[string(var)] = vals
        else   
            @error("Cannot export type "*string(Type)*" to VTK.")
        end
    end
    data.file[data.frame] = frame
    data.frame += 1
end

#=
"""
    import_particles!(sys::ParticleSystem, path::String, particle_constructor::Function)

Imports particles from a vtk file in 'path' using 'constructor'.
"""
function import_particles!(sys::ParticleSystem, path::String, particle_constructor::Function)
    input = read_vtk(path)
    (_, N) = size(input.point_coords)
    resize!(sys.particles, N)
    for i in 1:N
        x = [input.point_coords[1,i] for i in 1:3] 
        sys.particles[i] = particle_constructor(RealVector(x))
        for field in fieldnames(sys.particle_type)
            key = string(field)
            #skip particle data not present in file
            if !haskey(input.point_data, key)
                continue
            end
            Type = attribute_type(get_particle_type(sys), field)
            if Type <: Number
                val = input.point_data[key][i]
                setproperty!(sys.particles[i], field, val)
            elseif Type <: RealVector
                val = Type([input.point_data[key][j,i] for j in 1:3])
                setproperty!(sys.particles[i], field, val)
            elseif Type <: RealMatrix
                val = Type([input.point_data[key][j,k] for j in 1:3, k in 1:3])
                setproperty!(sys.particles[i], field, val)
            else
                throw("Cannot import data field with type:" *string(Type))
            end
        end
    end
end
=#