# This function evaluates and returns the coefficients for constraints and cost which are used later in the MPC formulation
# Inputs are:
# oldTraj   -> contains information about previous trajectories and Inputs
# mpcCoeff  -> contains information about previous coefficient results
# posInfo   -> contains information about the car's current position along the track
# mpcParams -> contains information about the MPC formulation (e.g. Q, R)
# stateIn   -> current state of the car
# inputIn   -> last input to the system

# structure of oldTrajectory: 1st dimension = state number, 2nd dimension = step number (time equiv.), 3rd dimennsion = lap number

# z[1] = xDot
# z[2] = yDot
# z[3] = psiDot
# z[4] = ePsi
# z[5] = eY
# z[6] = s

function coeffConstraintCost(oldTraj::OldTrajectory, mpcCoeff::MpcCoeff, posInfo::PosInfo, mpcParams::MpcParams,currentTraj::Array{Float64},currentInput::Array{Float64},lapStatus::LapStatus)
    # this computes the coefficients for the cost and constraints

    # Outputs: 
    # coeffConst
    # coeffCost

    # Read Inputs
    s_start         = posInfo.s_start
    s               = posInfo.s
    s_target        = posInfo.s_target


    # Parameters
    N               = mpcParams.N
    nz              = mpcParams.nz
    R               = mpcParams.R
    Order           = mpcCoeff.order                # interpolation order for cost and constraints
    pLength         = mpcCoeff.pLength              # interpolation length for polynomials

    n_prev          = 10                            # number of points before current s for System ID
    n_ahead         = 20                            # number of points ahead of current s for System ID

    coeffCost       = zeros(Order+1,2)              # polynomial coefficients for cost
    coeffConst      = zeros(Order+1,2,5)            # nz-1 beacuse no coeff for s

    # Select the old data
    oldxDot         = oldTraj.oldTraj[:,1,:]::Array{Float64,3}
    oldyDot         = oldTraj.oldTraj[:,2,:]::Array{Float64,3}
    oldpsiDot       = oldTraj.oldTraj[:,3,:]::Array{Float64,3}
    oldePsi         = oldTraj.oldTraj[:,4,:]::Array{Float64,3}
    oldeY           = oldTraj.oldTraj[:,5,:]::Array{Float64,3}
    oldS            = oldTraj.oldTraj[:,6,:]::Array{Float64,3}
    olda            = oldTraj.oldInput[:,1,:]::Array{Float64,3}
    olddF           = oldTraj.oldInput[:,2,:]::Array{Float64,3}

    N_points        = size(oldTraj.oldTraj,1)     # second dimension = length (=buffersize)

    s_total::Float64        # initialize
    DistS::Array{Float64}   # initialize
    idx_s::Array{Int64}     # initialize
    idx_s_target        = 0
    dist_to_s_target    = 0
    qLength             = 0
    vec_range::Tuple{UnitRange{Int64},UnitRange{Int64}}
    bS_Vector::Array{Float64}
    s_forinterpy::Array{Float64}

    # Compute the total s (current position along track)
    s_total = (s_start + s) % s_target

    # Compute the index
    DistS = ( s_total - oldS ).^2

    idx_s = findmin(DistS,1)[2]              # contains both indices for the closest distances for both oldS !!

    vec_range = (idx_s[1]:idx_s[1]+pLength,idx_s[2]:idx_s[2]+pLength)

    # Create the vectors used for the interpolation
    # bS_vector contains the s-values for later interpolation (total s, s_start is subtracted later)
    bS_Vector       = zeros(pLength+1,2)
    for i=1:pLength+1
        bS_Vector[i,1] = oldS[vec_range[1][i]]
        bS_Vector[i,2] = oldS[vec_range[2][i]]
    end
    if norm(bS_Vector[1,1]-s_total) > 0.3 || norm(bS_Vector[1,2]-s_total) > 0.3
        warn("Couldn't find a close point to current s.")
    end

    # println("************************************** COEFFICIENTS **************************************")
    # println("idx_s[1]  = $(idx_s[1]), idx_s[2] = $(idx_s[2])")
    # println("s_total   = $s_total")
    # println("bS_Vector[1,:] = $(bS_Vector[1,:]')")
    # These matrices (above) contain two vectors each (for both old trajectories), stored in the 3rd dimension
    
    # The states are parametrized with resprect to the curvilinear abscissa,
    # so we select the point used for the interpolation. Need to subtract an
    # offset to be coherent with the MPC formulation
    s_forinterpy   = bS_Vector - s_start
    if s_total - s_start < 0
        s_forinterpy += s_target
    end
    # println("s_forinterpy[:,:]' = $(s_forinterpy[:,:]')")
    # Create the Matrices for the interpolation
    MatrixInterp = zeros(pLength+1,Order+1,2)

    for k = 0:Order
        MatrixInterp[:,Order+1-k,:] = s_forinterpy[:,:].^k
    end
    
    # Compute the coefficients
    coeffConst = zeros(Order+1,2,5)
    for i=1:2
        coeffConst[:,i,1]    = MatrixInterp[:,:,i]\oldxDot[vec_range[i]]
        coeffConst[:,i,2]    = MatrixInterp[:,:,i]\oldyDot[vec_range[i]]
        coeffConst[:,i,3]    = MatrixInterp[:,:,i]\oldpsiDot[vec_range[i]]
        coeffConst[:,i,4]    = MatrixInterp[:,:,i]\oldePsi[vec_range[i]]
        coeffConst[:,i,5]    = MatrixInterp[:,:,i]\oldeY[vec_range[i]]
    end

    # xDotInterp = [s_forinterpy[:,1].^5 s_forinterpy[:,1].^4 s_forinterpy[:,1].^3 s_forinterpy[:,1].^2 s_forinterpy[:,1].^1 s_forinterpy[:,1].^0]*coeffConst[:,1,1]
    # plot(s_forinterpy,oldxDot[vec_range[1]],"*",s_forinterpy,xDotInterp)
    # title("xDotInterp")
    # grid("on")
    # readline()

    # Finished with calculating the constraint coefficients
    
    # Now compute the final cost coefficients

    # The Q-function contains for every point in the sampled safe set the minimum cost-to-go-value
    # These values are calculated for both old trajectories
    # The vector bQfunction_Vector contains the cost at each point in the interpolated area to reach the finish line
    # From this vector, polynomial coefficients coeffCost are calculated to approximate this cost
    for i=1:2   
            dist_to_s_target  = oldTraj.oldCost[i] - (idx_s[i]-N_points*(i-1))  # number of iterations from idx_s to s_target
            bQfunction_Vector = collect(linspace(dist_to_s_target,dist_to_s_target-1,pLength+1))    # build a vector that starts at the distance and
                                                                                                    # decreases in equal steps
            coeffCost[:,i]    = MatrixInterp[:,:,i]\bQfunction_Vector           # interpolate this vector with the given s
    end

    # --------------- SYSTEM IDENTIFICATION --------------- #
    # ----------------------------------------------------- #

    n_ID            = n_prev+n_ahead+1                  # number of system ID points (using only 1 round)
    vec_range_ID    = (idx_s[1]-n_prev:idx_s[1]+n_ahead,idx_s[2]-n_prev:idx_s[2]+n_ahead)  # related index range
    vec_range_ID2   = size(currentTraj,1)-n_prev:size(currentTraj,1)-1      # index range of previous states and inputs

    # psiDot
    y_psi = diff(oldpsiDot[idx_s[1]-n_prev:idx_s[1]+n_ahead+1])
    A_psi = [oldpsiDot[vec_range_ID[1]]./oldxDot[vec_range_ID[1]] oldyDot[vec_range_ID[1]]./oldxDot[vec_range_ID[1]] olddF[vec_range_ID[1]]]
    y_psi = cat(1,y_psi,diff(currentTraj[end-n_prev:end,3]))
    A_psi = cat(1,A_psi,[currentTraj[vec_range_ID2,3]./currentTraj[vec_range_ID2,1] currentTraj[vec_range_ID2,2]./currentTraj[vec_range_ID2,1] currentInput[vec_range_ID2,2]])
 
    # xDot
    y_xDot = diff(oldxDot[idx_s[1]-n_prev:idx_s[1]+n_ahead+1])
    A_xDot = [oldyDot[vec_range_ID[1]] oldpsiDot[vec_range_ID[1]] oldxDot[vec_range_ID[1]].^2 olda[vec_range_ID[1]]]
    y_xDot = cat(1,y_xDot,diff(currentTraj[end-n_prev:end,1]))
    A_xDot = cat(1,A_xDot,[currentTraj[vec_range_ID2,2] currentTraj[vec_range_ID2,3] currentTraj[vec_range_ID2,1] currentInput[vec_range_ID2,1]])

    # yDot
    y_yDot = diff(oldyDot[idx_s[1]-n_prev:idx_s[1]+n_ahead+1])
    A_yDot = [oldyDot[vec_range_ID[1]]./oldxDot[vec_range_ID[1]] oldpsiDot[vec_range_ID[1]].*oldxDot[vec_range_ID[1]] oldpsiDot[vec_range_ID[1]]./oldxDot[vec_range_ID[1]] olddF[vec_range_ID[1]]]
    y_yDot = cat(1,y_yDot,diff(currentTraj[end-n_prev:end,2]))
    A_yDot = cat(1,A_yDot,[currentTraj[vec_range_ID2,2]./currentTraj[vec_range_ID2,1] currentTraj[vec_range_ID2,3].*currentTraj[vec_range_ID2,1] currentTraj[vec_range_ID2,3]./currentTraj[vec_range_ID2,1] currentInput[vec_range_ID2,2]])

    if diff(oldTraj.oldCost) != 0       # if the two saved old trajectories are different (resp. have different costs)
        y_psi  = cat(1,y_psi,diff(oldpsiDot[idx_s[2]-n_prev:idx_s[2]+n_ahead+1]))
        A_psi  = cat(1,A_psi,[oldpsiDot[vec_range_ID[2]]./oldxDot[vec_range_ID[2]] oldyDot[vec_range_ID[2]]./oldxDot[vec_range_ID[2]] olddF[vec_range_ID[2]]])
        y_xDot = cat(1,y_xDot,diff(oldxDot[idx_s[2]-n_prev:idx_s[2]+n_ahead+1]))
        A_xDot = cat(1,A_xDot,[oldyDot[vec_range_ID[2]] oldpsiDot[vec_range_ID[2]] oldxDot[vec_range_ID[2]].^2 olda[vec_range_ID[2]]])
        y_yDot = cat(1,y_yDot,diff(oldyDot[idx_s[2]-n_prev:idx_s[2]+n_ahead+1]))
        A_yDot = cat(1,A_yDot,[oldyDot[vec_range_ID[2]]./oldxDot[vec_range_ID[2]] oldpsiDot[vec_range_ID[2]].*oldxDot[vec_range_ID[2]] oldpsiDot[vec_range_ID[2]]./oldxDot[vec_range_ID[2]] olddF[vec_range_ID[2]]])
    end

    #println("y_yDot:")
    #println(y_yDot)
    #println("y_psi:")
    #println(y_psi)

    # figure(6)
    # subplot(511)
    # plot(A_xDot[:,1],"-o")
    # legend(["yDot"])
    # grid("on")
    # subplot(512)
    # plot(A_xDot[:,2],"-o")
    # legend(["psiDot"])
    # grid("on")
    # subplot(513)
    # plot(A_xDot[:,3],"-o")
    # legend(["xDot.^2"])
    # grid("on")
    # subplot(514)
    # plot(A_xDot[:,4],"-o")
    # legend(["a"])
    # grid("on")
    # subplot(515)
    # plot(y_xDot,"-o")
    # legend(["y_xDot"])
    # grid("on")

    # figure(7)
    # subplot(511)
    # plot(A_yDot[:,1],"-o")
    # legend(["y/x"])
    # grid("on")
    # subplot(512)
    # plot(A_yDot[:,2],"-o")
    # legend(["psi*x"])
    # grid("on")
    # subplot(513)
    # plot(A_yDot[:,3],"-o")
    # legend(["psi/x"])
    # grid("on")
    # subplot(514)
    # plot(A_yDot[:,4],"-o")
    # legend(["d_f"])
    # grid("on")
    # subplot(515)
    # plot(y_yDot,"-o")
    # legend(["y_yDot"])
    # grid("on")
    # readline()

    if any(isnan,y_yDot)            # check if any value in the y_yDot value is NaN
        warn("NaN value detected in coeffConstraintCost! Press to continue...")
        readline()
    end
    mpcCoeff.c_Psi = (A_psi'*A_psi)\A_psi'*y_psi
    mpcCoeff.c_Vx  = (A_xDot'*A_xDot)\A_xDot'*y_xDot
    mpcCoeff.c_Vy  = (A_yDot'*A_yDot)\A_yDot'*y_yDot        # Todo: If all matrix/vector values are really low (<e-22) this might return an error
                                                            # this might happen on long straight parts of the track (no change in psi/y states)
    mpcCoeff.coeffCost  = coeffCost
    mpcCoeff.coeffConst = coeffConst
    nothing
end

# Notes about oldTrajectory:
# oldTrajectory[1:prebuf] = states before finish line (s < 0)
# oldTrajectory[prebuf+1:prebuf+cost] = states between start and finish line (0 <= s < s_target)
# oldTrajectory[prebuf+cost+1:prebuf+cost+postbuf] = states after finish line (s_target <= s)
# once one lap is over, the states are saved (interval prebuf:s_target)
# during the next lap, the states behind the finish line are appended (interval s_target:postbuf)