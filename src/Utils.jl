module Utils
using Parameters

export InputParameters
export InputFiles
export unpack_selected
export unpack_selected_filenames

"Define struct with constant parameters"
@with_kw struct InputParameters
    voxel_length::Float64 = 7e-6 # Length of voxel
    Nx::Int=300
    Ny::Int=297
    Nz::Int=500
    L::Float64 = Nz*voxel_length; # Length of domain
    dv::Float64 = voxel_length/1 # Minimum velocity increment, voxel length divided by one simulation time step
    len_vector::Int = Nx*Ny*Nz
end

"Define struct with input files"
struct InputFiles
    infile_raw::String  # Input from LBPM, 8-bit binary images savedafter simulations 
    # infile_raw_test::String = "id_t10000.raw" # Input from LBPM, 8-bit binary images savedafter simulations 
    input_files_velocities::String
end

"Extract selected parameters from structs."
unpack_selected(p::InputParameters, fields...) = map(x->getfield(p, x), fields)

"Extract filenames from struct"
unpack_selected_filenames(p::InputFiles, fields...) = map(x->getfield(p, x), fields)

end
