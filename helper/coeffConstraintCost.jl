# This function evaluates and returns the coefficients for constraints and cost which are used later in the MPC formulation
# Inputs are:
# oldTraj   -> contains information about previous trajectories and Inputs
# mpcCoeff  -> contains information about previous coefficient results
# posInfo   -> contains information about the car's current position along the track
# mpcParams -> contains information about the MPC formulation (e.g. Q, R)
# stateIn   -> current state of the car
# inputIn   -> last input to the system

# structure of oldTrajectory: 1st dimension = state number, 2nd dimension = step number (time equiv.), 3rd dimennsion = lap number

function coeffConstraintCost!(oldTraj::OldTrajectory, mpcCoeff::MpcCoeff, posInfo::PosInfo, mpcParams::MpcParams)
    # this computes the coefficients for the cost and constraints

    # Outputs: 
    # coeffConst
    # coeffCost

    # Read Inputs
    s_start         = posInfo.s_start #currently start alqys 0 , normally start current postion and s till beginn of pred horizon
    s               = posInfo.s 
    s_target        = posInfo.s_target


    # Parameters
    N               = mpcParams.N
    nz              = mpcParams.nz
    R               = mpcParams.R
    Order           = mpcCoeff.order                # interpolation order for cost and constraints

    pLength         = mpcCoeff.pLength              # interpolation length for polynomials

    coeffCost       = zeros(Order+1,2)            # polynomial coefficients for cost, second dimension for number of old trajectories
    coeffConst      = zeros(Order+1,2,3)          # nz-1 beacuse no coeff for s

    # Select the old data
    oldS            = oldTraj.oldTraj[:,1,:]::Array{Float64,2}
    oldeY           = oldTraj.oldTraj[:,2,:]::Array{Float64,2}
    oldePsi         = oldTraj.oldTraj[:,3,:]::Array{Float64,2}
    oldV            = oldTraj.oldTraj[:,4,:]::Array{Float64,2}

    N_points        = size(oldTraj.oldTraj,1)     #  dimension = length #how many points we saved so that we can subract these indeces to get values for second stroed trajectory which beginns at index NPoints+1

    local s_total::Float64        # initialize
    local DistS::Array{Float64}   # initialize
    local idx_s::Array{Int64}     # initialize
    idx_s_target        = 0
    dist_to_s_target    = 0
    qLength             = 0
    local vec_range::Tuple{UnitRange{Int64},UnitRange{Int64}}
    local bS_Vector::Array{Float64}
    local s_forinterpy::Array{Float64}

    # Compute the total s (current position along track)
    #we dont need that anymore 
    s_total =  s % s_target #?? what do we calculate total_s for one round?

    # Compute the index
    DistS = ( s_total - oldS ).^2

    idx_s = findmin(DistS,1)[2]              # contains both indices for the closest distances for both oldS !!
    
    vec_range = (idx_s[1]:idx_s[1]+pLength,idx_s[2]:idx_s[2]+pLength)

    # Create the vectors used for the interpolation
    bS_Vector       = zeros(pLength+1,2)
    for i=1:pLength+1
        bS_Vector[i,1] = oldS[vec_range[1][i]]
        bS_Vector[i,2] = oldS[vec_range[2][i]]
    end
    # bS_Vector       = cat(2, oldS[vec_range[1]],    oldS[vec_range[2]])
    
    # println("************************************** COEFFICIENTS **************************************")
    # println("idx_s[1]  = $(idx_s[1]), idx_s[2] = $(idx_s[2])")
    # println("s_total   = $s_total")
    # println("bS_Vector[1] = $(bS_Vector[:,:,1]')")
    # These matrices (above) contain two vectors each (for both old trajectories), stored in the 3rd dimension
    
    # The states are parametrized with resprect to the curvilinear abscissa,
    # so we select the point used for the interpolation. Need to subtract an
    # offset to be coherent with the MPC formulation
    s_forinterpy   = bS_Vector - s_start 
    if s_total - s_start < 0
        s_forinterpy += s_target
    end
    # println("s_forinterpy[:,1,1]' = $(s_forinterpy[:,1,1]')")
    # Create the Matrices for the interpolation
    MatrixInterp = zeros(pLength+1,Order+1,2)

    for k = 0:Order
        MatrixInterp[:,Order+1-k,:] = s_forinterpy[:,:].^k
    end
    # Compute the constraint coefficients for both old trajectories
    
    coeffConst = zeros(Order+1,2,3)

    for i=1:2
        coeffConst[:,i,1]    = MatrixInterp[:,:,i]\oldeY[vec_range[i]]
        coeffConst[:,i,2]    = MatrixInterp[:,:,i]\oldePsi[vec_range[i]]
        coeffConst[:,i,3]    = MatrixInterp[:,:,i]\oldV[vec_range[i]]
    end

    # Finished with calculating the constraint coefficients
    
    # Now compute the final cost coefficients

    # The Q-function contains for every point in the sampled safe set the minimum cost-to-go-value
    # These values are calculated for both old trajectories
    # The vector bQfunction_Vector contains the cost at each point in the interpolated area to reach the finish line
    # From this vector, polynomial coefficients coeffCost are calculated to approximate this cost

    for i=1:2   
        # in this part we construct a polynomial for the cost associated at each s value. the s value at the curent postion is used to calculate the cost of the old round at this postion
        #we know that each following s has a cost value which is exactl 1 less thne the one before. so we can easiyl do the least squares to get the coeficients for approximation
        iter_to_s_target  = oldTraj.oldCost[i] - (idx_s[i]-N_points*(i-1))  # number of iterations from idx_s to s_target, this has sth todo with the count in the array as we look at values in second row
        bQfunction_Vector = collect(linspace(iter_to_s_target,iter_to_s_target-pLength,pLength+1))    # build a vector that starts at the distance and decreases in equal steps
        coeffCost[:,i]    = MatrixInterp[:,:,i]\bQfunction_Vector           # interpolate this vector with the given s
    end

    mpcCoeff.coeffCost  = coeffCost #this value goes into the variable mpcCoeff in testNode.jl as well variables by reference
    mpcCoeff.coeffConst = coeffConst #this way we dont need to return anything #??speed advantage ? looks weird
    
    nothing
end
