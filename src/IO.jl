using WriteVTK
using VTKDataIO

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
    for k in 1:N, i in 1:2
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
        type = attribute_type(sys.particle_type, var)
        N = length(sys.particles)
        if type <: Number
            frame[string(var)] = field
        elseif type == Vec2
            vals = zeros(Float64, 3, N)
            for k in 1:N, i in 1:2
                vals[i,k] = field[k][i]
            end
            frame[string(var)] = vals
        elseif type == Mat2
            vals = zeros(Float64, 3, 3, N)
            for k in 1:N, i in 1:2, j in 1:2
                vals[i,j,k] = field[k][i,j]
            end
        else   
            @error("Cannot export type "*string(type)*" to VTK.")
        end
    end
    data.file[data.frame] = frame
    data.frame += 1
end

function attribute_type(type::DataType, var::Symbol)
    ind = findfirst(s -> s == var, fieldnames(type))
    if typeof(ind) == Nothing
        @error("Variable "*string(var)* 
            " cannot be exported to VTK because it does not exist")
    end
    return fieldtypes(type)[ind]
end

"""
    import_particles!(sys::ParticleSystem, path::String, particle_constructor::Function)

Imports particles from a vtk file in 'path' using 'constructor'.
"""
function read_vtk!(sys::ParticleSystem, path::String, particle_constructor::Function)
    input = read_vtk(path)
    (_, N) = size(input.point_coords)
    resize!(sys.particles, N)
    for i in 1:N
        x1 = input.point_coords[1,i]
        x2 = input.point_coords[2,i]
        sys.particles[i] = particle_constructor(Vec2(x1, x2))
        for field in fieldnames(sys.particle_type)
            key = string(field)
            type = attribute_type(sys.particle_type, field)
            if type <: Number && haskey(input.point_data, key)
                setproperty!(sys.particles[i], field, input.point_data[key][i])
            end
            if type == Vec2 && haskey(input.point_data, key)
                setproperty!(sys.particles[i], field, Vec2(input.point_data[key][1,i], input.point_data[key][2,i]))
            end
        end
    end
end