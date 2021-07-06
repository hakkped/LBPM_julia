module CalcSaturationGradients
using PyCall

export FindMaxVelocities
"""
    FindMaxVelocities()

Create list of velocity-files with different saturations read from folder, find max/min velocities, and return velocity increments.
In: Filenames of all velocity-files.
Out: Velocity-increments/vector created from max/min velocities in all samples.

"""
function FindMaxVelocities()
    # Open all files, find max/min velocity and create velocity increments
    dv::Float64=7e-6
    filenames=readdir("velocities/")
    number_samples = length(filenames)
    max_vector=zeros(Float64,number_samples)
    min_vector=zeros(Float64,number_samples)
    py"""
    import pyvisfile.silo as silo
    import numpy as np
    """
    j=1
    for i in filenames
        file =string("velocities/",i)
        py"""
        db = silo.SiloFile($file, create=False, mode=silo.DB_READ)
        f = db.get_quadvar("Velocity_z")
        global max_sample
        global min_sample
        data_velocities=f.vals[0][0]
        max_sample = np.max(data_velocities)    
        min_sample = np.min(data_velocities)
        """
        max_sample=py"max_sample"
        min_sample=py"min_sample"
        max_vector[j]=max_sample
        min_vector[j]=min_sample
        j += 1
    end # i
    max_all=maximum(max_vector)
    min_all=minimum(min_vector)
    velocity_increments=collect(min_all:dv:max_all)
    return velocity_increments
end # FindMaxVelocity

"""
Calculate gradients ∂/∂S_w and ∂/∂S_n and thereby the expression: a_m(v) = S_w*∂/∂S_w(a_w/S_w) + S_n*∂/∂S_w(a_n/S_n)
Input: ap_all,aw_all,an_all,saturations_all
Output: a_m(v)
"""
function CalcGradients(velocity_increments::Vector{Float64},ap_all::Array{Float64,2},aw_all::Array{Float64,2},an_all::Array{Float64,2},saturations_all::Array{Float64,2})
    # Nv = length(velocity_increments)
    # Sw = zeros(Float64)
    # Sn
    # for i in 1:Nv

    #     a_m=
    # end

end # CalcGradients

end # Module
