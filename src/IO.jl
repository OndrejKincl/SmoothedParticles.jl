using WriteVTK
using ReadVTK

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


"""
    import_particles!(sys::ParticleSystem, path::String, particle_constructor::Function)

Imports particles from a vtk file in 'path' using 'constructor' and add them to `sys`.
"""
function import_particles!(sys::ParticleSystem, path::String, particle_constructor::Function)
    vtk = VTKFile(path)
    data = get_point_data(vtk)
    points = get_points(vtk)
    N0 = length(sys.particles)
    N = vtk.n_points
    resize!(sys.particles, N0 + N)
    PType = get_particle_type(sys)
    for i in 1:N
        x = [points[j,i] for j in 1:3] 
        sys.particles[N0 + i] = particle_constructor(RealVector(x))
    end
    for fieldname in fieldnames(PType)
        fieldvals = []
        try
            fieldvals = get_data(data[string(fieldname)])
        catch
            continue
        end
        Type = SmoothedParticles.attribute_type(PType, fieldname)
        if Type <: Number
            for i in 1:N
                val = fieldvals[i]
                setproperty!(sys.particles[N0 + i], fieldname, val)
            end
        elseif Type <: RealVector
            for i in 1:N
                val = Type([fieldvals[j,i] for j in 1:3])
                setproperty!(sys.particles[N0 + i], fieldname, val)
            end
        elseif Type <: RealMatrix
            for i in 1:N
                val = Type([fieldvals[j,i] for j in 1:9])
                setproperty!(sys.particles[N0 + i], fieldname, val)
            end
        else
            error("Cannot import data field with type:" *string(Type))
        end
    end
end