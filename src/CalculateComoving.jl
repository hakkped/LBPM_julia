"""
Calculate the co-moving velocity from real pore space images using LBPM.
"""

module CalcSingleSample
using  Gnuplot, BenchmarkTools, LinearAlgebra, Plots, Images, Colors, PyCall, Distributions, StatsBase
export ap_dv
export aw_dv
export an_dv
export velocity_increments

infile_raw = "mask_water_flooded_water_and_oil.raw" # Input from LBPM, 8-bit binary images savedafter simulations 
infile_raw_test = "id_t10000.raw" # Input from LBPM, 8-bit binary images savedafter simulations 
infile_ξ = "infile_travel.raw" # Input from Visit software, 
infile_velocities = "infile_velocities.raw" # Velocities, to be sorted 
path_infiles = "/mnt/****" #

# Definitions & parameters

# Nx=601; Ny=594; Nz=1311 # Array dimension of input
const voxel_length = 7e-6 # Length of voxel
const Nx=300; const Ny=297; const Nz=500
const L = Nz*voxel_length; # Length of domain
const dv = voxel_length/1 # Minimum velocity increment, voxel length divided by one simulation time step


"Create arrays from input data"
len_vector = Nx*Ny*Nz
data_infile_raw = Vector{UInt8}(undef,len_vector); read!(infile_raw_test, data_infile_raw)

# Read velocities from SILO database using embedded python and the pyvisfile-module
println("Importing velocity from .silo-file...")
py"""
import pyvisfile.silo as silo 
db = silo.SiloFile("00000.silo", create=False, mode=silo.DB_READ)
f = db.get_quadvar("Velocity_z")
global data_velocities
data_velocities=f.vals[0][0]
"""
data_velocities = py"data_velocities" # All velocities in all voxels (z-component)
println("Import succeeded!")
# The data can be plotted using heatmap(reshape(A,(300,297,500))[:,:,500], c=:inferno)

"Reshape data (might not be needed)"
data_infile_raw_rs = reshape(data_infile_raw, Nx,Ny,Nz)

# reinterpret(N0f8, file_ID)
 # file_ID_reshaped = reshape(file_ID[(1*Nx*Ny+1):(2*Nx*Ny)], Nx,Ny)
# println(size(file_ID))
# colorview(Gray, file_ID_reshaped./2)

# plot(file_ID_reshaped)
# plot(Gray.(file_ID_reshaped), show=true)

# Working plotting of 2D images of type Uint8"
# Convert to Int for calculations with the values."
# data_infile_raw_rs_2D=view(data_infile_raw_rs, :,:,100) # These work for plotting 
# plot(colorview(Gray,data_infile_raw_rs_2D./2), show=true)

# img_raw = rawview(file_ID_reshaped)
# println(float.(img_raw))
# plot(Gray.(reinterpret(UInt8,file_ID_reshaped)), show=true)


raw_vector = vec(Int8.(data_infile_raw_rs))   # Convert simulation result raw data to vectors of Int8 and enumerate placementent. [value index]

# Sort velocities in entire array
vel_vector= vec(Float64.(data_velocities)) # Convert array to a vector of Ints
velocity_increments=collect(minimum(vel_vector):dv:maximum(vel_vector))
Nv  = length(velocity_increments)

# Differential area distributions to be integrated over
ap_dv = Vector{Float64}(undef,Nv-1); fill!(ap_dv,0.0)
aw_dv = Vector{Float64}(undef,Nv-1); fill!(aw_dv,0.0)
an_dv = Vector{Float64}(undef,Nv-1); fill!(an_dv,0.0)


# Loop over arrays
vel_set = zeros(Float64,Nx*Ny) # Set of velocities satisfying the condition in a layer.
vel_hist_all = Array{Int64}(undef,(Nv-1,Nz)) # A histogram using a vector with length Nv to generate bins will have Nv-1 bins.
for j in 1:Nz # Loop over layers
println(j)
    # This part was written with a single loop over all ites in mind. It seems like this is not the most convenient way.
   # slice_number = Int.(1+mod(floor(j/(Nx*Ny+0.1)), Nx*Ny)) # Find slice number s.t. the velocities in this layer can be found
   # slice_indices = ((slice_number-1)*Nx*Ny+1):(slice_number*Nx*Ny) # Indices in slice used to index the long vector
   # vel_positions = collect(Int8,((slice_number-1)*Nx*Ny+1):(slice_number*Nx*Ny)) # Generate indices in a layer
   # vel_layer = vel_vector_coord[vel_positions,:] # Collection of velocities in a layer. [velocities indices]
   # vel_layer_sorted = vel_layer[sortperm(vel_layer[:,1]),:] # Sort the velocities in the layer and reshuffle indices accordingly
    # End commented section.

   slice_indices = ((j-1)*Nx*Ny+1):(j*Nx*Ny) # Indices in slice used to index the long vector

   # List of indices for wetting and nonwetting components in slice
   w_indices=findall(x->x==1,raw_vector[slice_indices])
   nw_indices=findall(x->x==2,raw_vector[slice_indices])
   
   "Velocities in slice"
   vel_layer=vel_vector[slice_indices] # Pick out velocities in a single layer
   vel_hist=fit(Histogram, filter(x->x≠0.0,vel_layer), velocity_increments) # Creates histogram of velocities in layer, except for the velocities that are exactly zero (solid).
   vel_set=vel_hist.weights # Get weights/bin values

   vel_layer_w=vel_layer[w_indices]
   vel_hist_w=fit(Histogram, vel_layer_w, velocity_increments)  # Wetting velocity histogram
   vel_set_w=vel_hist_w.weights
    
   vel_layer_nw=vel_layer[nw_indices]
   vel_hist_nw=fit(Histogram, vel_layer_nw, velocity_increments)  # Nonwetting velocity histogram
   vel_set_nw=vel_hist_nw.weights

   # for k in 1:Nx*Ny # Loop on a single layer
    # vel_set = vel_layer_sorted[velocity_increments[i+1] .> vel_layer_sorted[k,1] .> velocity_increments[i], :] # Set of velocities in slice that satisfies condition, [velocity index]
       # if (velocity_increments[i+1] > vel_layer[k] > velocity_increments[i])
           # vel_set = filter(t -> velocity_increments[i] < t < velocity_increments[i+1], vel_layer_sorted[:,1])
           # vel_set[k]=vel_layer[k]
           # Below solution might not work. findall returns indices of the found values.
           # vel_set = vel_layer_sorted[findall(<(velocity_increments[i+1]),findall(>(velocity_increments[i],vel_layer_sorted))),:]
           # ξ_set[k] = ξ_vector[k] # Get set of flows projected onto z in the same locations.
   vel_hist_all[:,j] = vel_set # Store the weights in array

   # Calculate saturations in slice
   solid=countmap(raw_vector[slice_indices])[0] # Number of vozels of solid, nw and w type.
   nw=countmap(raw_vector[slice_indices])[2]; nw_frac=nw/solid
   w=countmap(raw_vector[slice_indices])[1]; w_frac=w/solid
   porosity=(nw+w)/solid # Since all of pore space is filled w. fluid, porosity is give by the sum of w and nw type fluids.

    "Calculate transversal area distribution as a function of v"
    @. ap_dv += (1/L)*vel_set*voxel_length^2
    @. aw_dv += (1/L)*vel_set_w*voxel_length^2 
    @. an_dv += (1/L)*vel_set_nw*voxel_length^2

end

# Calculate A_m = ∫_{-∞}^{∞} dv_z a_m = 0

# Plotting
# sort(collect(Iterators.flatten(vel_layers_all)))

end # Module
