using WriteVTK

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
        println("created new path: "*path)
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
    for k in 1:N
        points[1,k] = sys.particles[k].x
        points[2,k] = sys.particles[k].y
    end
    cells = [MeshCell(PolyData.Verts(), [i]) for i in 1:N]
    vtk_file = vtk_grid(data.path*"/frame"*string(data.frame), points, cells)
    return vtk_file
end

function save_frame!(data::DataStorage, vtk_file)
    data.file[data.frame] = vtk_file
    data.frame += 1
end

"""
    save_frame!(sys::ParticleSystem, data::DataStorage, vars::DataField...)

Inserts one time frame into a `DataStorage` that includes all `vars...` as fields.
"""
function save_frame!(sys::ParticleSystem, data::DataStorage, vars::DataField...)
    frame = capture_frame(sys, data)
    for var in vars
        frame[var.name] = var
    end
    data.file[data.frame] = frame
    data.frame += 1
end
