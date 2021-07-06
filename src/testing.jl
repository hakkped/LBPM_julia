
"Used for testing in the REPL, import the database"

using  Gnuplot, BenchmarkTools, LinearAlgebra, Plots, Images, Colors, PyCall, Distributions

const voxel_length = 7e-6 # Length of voxel 
const Nx=300; const Ny=297; const Nz=500
const L = Nz*voxel_length; # Length of domain
const dv = voxel_length/1 # Minimum velocity increment, voxel length divided by one simulation time step

# Read in raw-file

# This part imports the database and generates the array of all velocities and raw-file
# py"""
# import pyvisfile.silo as silo 
# db = silo.SiloFile("00000.silo", create=False, mode=silo.DB_READ)
# f = db.get_quadvar("Velocity_z")
# global A
# A=f.vals[0][0]
# """
# A = py"A" # All velocities in all voxels (z-component)

A_vector= vec(Float64.(A))
A_increments=collect(minimum(A_vector):dv:maximum(A_vector))
