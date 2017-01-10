using JuMP
using Ipopt
using PyPlot
using JLD

include("barc_lib/classes.jl")
include("barc_lib/LMPC/MPC_models.jl")
include("barc_lib/LMPC/functions.jl")
include("barc_lib/LMPC/coeffConstraintCost.jl")
include("barc_lib/LMPC/solveMpcProblem.jl")
include("barc_lib/simModel.jl")


# IMPORTANT GENERAL DEFINITION:
# x_i+1 = f(x_i,u_i)
# At each timestep i the control inputs u_i and the next state x_i+1 is calculated
# Old trajectory stores x_i and u_i in the same row, even though the u_i is applied in this timestep (and realized in the next step)
# The final state of a lap (at s >= s_target) is used as the initial step in the next lap

function run_sim()
    # DEFINE PARAMETERS
    # Define and initialize variables
    oldTraj                     = OldTrajectory()
    posInfo                     = PosInfo()
    mpcCoeff                    = MpcCoeff()
    lapStatus                   = LapStatus(1,1,false,false,0.3)
    mpcSol                      = MpcSol()
    trackCoeff                  = TrackCoeff()      # info about track (at current position, approximated)
    modelParams                 = ModelParams()
    mpcParams                   = MpcParams()
    mpcParams_pF                = MpcParams()       # for 1st lap (path following)

    buffersize                  = 5000

    InitializeParameters(mpcParams,mpcParams_pF,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,buffersize)
    mdl    = MpcModel(mpcParams,mpcCoeff,modelParams,trackCoeff)
    mdl_pF = MpcModel_pF(mpcParams_pF,modelParams,trackCoeff)

    # Simulation parameters
    dt                          = modelParams.dt::Float64
    t                           = collect(0:dt:100)          # time vector
    zCurr                       = zeros(length(t),6)        # these are the simulated states
    uCurr                       = zeros(length(t),2)        # these are the inputs
    cost                        = zeros(length(t),6)        # these are MPC cost values at each step

    log_z                       = zeros(length(t),6,30)     # log z for 30 laps
    log_u                       = zeros(length(t),2,30)     # log u for 30 laps
    log_xy                      = zeros(length(t),2,30)     # log x-y-data for 30 laps
    log_ParInt                  = zeros(length(t),30)       # log ParInt data

    # Logging parameters
    posInfo.s_target            = 50.49

    z_final         = zeros(6)

    uPrev           = zeros(10,2)

    n_pf            = 3             # number of path-following laps

    s_track = 0.01:.01:50.49
    c_track = zeros(5049)
    c_track[1:300] = 0
    c_track[301:400] = linspace(0,-pi/2,100)
    c_track[401:500] = linspace(-pi/2,0,100)
    c_track[501:900] = 0
    c_track[901:1000] = linspace(0,-pi/2,100)
    c_track[1001:1100] = linspace(-pi/2,0,100)
    c_track[1101:1200] = linspace(0,-pi/4,100)
    c_track[1201:1300] = linspace(-pi/4,0,100)
    c_track[1301:1600] = 0
    c_track[1601:2100] = linspace(0,10*pi/4/10,500)
    c_track[2101:2600] = linspace(10*pi/4/10,0,500)
    c_track[2601:2900] = linspace(0,-pi/3,300)
    c_track[2901:3200] = linspace(-pi/3,0,300)
    c_track[3201:3500] = 0
    c_track[3501:3700] = linspace(0,-2*pi/2/4,200)
    c_track[3701:3900] = linspace(-2*pi/2/4,0,200)
    c_track[3901:4102] = 0
    c_track[4103:4402] = linspace(0,-2*pi/2/6,300)
    c_track[4403:4702] = linspace(-2*pi/2/6,0,300)

    s_track_p, c_track_p = prepareTrack(s_track, c_track)

    no_solution_found = 0

    # Run 10 laps
    for j=1:15
        # Initialize Lap
        lapStatus.currentLap = j

        tt          = zeros(length(t),1)
        zCurr       = zeros(length(t),6)
        uCurr       = zeros(length(t),2)
        zCurr[1,1]  = 0.2

        if j>1                                  # if we are in the second or higher lap
            zCurr[1,:]          = z_final       # use final state as initial state
            zCurr[1,6]          = zCurr[1,6]%posInfo.s_target   # and make sure that it's before the finish line (0 <= s < s_target)
        end

        cost        = zeros(length(t),6)
        finished    = false

        # Start one lap
        # --------------------------------
        i = 1
        while i<length(t) && !finished
            #println("///////////////////////////////// STARTING ONE ITERATION /////////////////////////////////")
            # Define track curvature
            trackCoeff.coeffCurvature = find_curvature(s_track_p,c_track_p,zCurr[i,6],trackCoeff)

            # Calculate coefficients for LMPC (if at least in the 2nd lap)
            posInfo.s   = zCurr[i,6]
            if j > n_pf
                coeffConstraintCost(oldTraj,mpcCoeff,posInfo,mpcParams,lapStatus)
            end

            # Calculate optimal inputs u_i (solve MPC problem)
            if j <= n_pf                    # if we are in the first x laps of path following
                z_pf = zCurr[i,[6,5,4,1]]
                solveMpcProblem_pathFollow(mdl_pF,mpcSol,mpcParams_pF,trackCoeff,posInfo,modelParams,z_pf',uPrev)
            else                            # otherwise: LMPC
                solstat = solveMpcProblem_LMPC(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr[i,:]',uPrev)
                if solstat == false
                    no_solution_found += 1
                else
                    no_solution_found = 0
                end
                if no_solution_found == 5
                    break
                end
                # Plot predictions
                # figure(101)
                # for k=1:5
                #     subplot(5,1,k)
                #     s = mpcSol.z[:,6]
                #     poly1 = [s.^5 s.^4 s.^3 s.^2 s.^1 s.^0]*mpcCoeff.coeffConst[:,1,k]
                #     poly2 = [s.^5 s.^4 s.^3 s.^2 s.^1 s.^0]*mpcCoeff.coeffConst[:,2,k]
                #     plot(mpcSol.z[:,6],mpcSol.z[:,k],"-o")
                #     plot(s,poly1)
                #     plot(s,poly2)
                #     grid("on")
                # end
                # readline()
            end
            cost[i,:]   = mpcSol.cost

            # Simulate the model -> calculate x_i+1 = f(x_i, u_i)
            uCurr[i,:]          = [mpcSol.a_x mpcSol.d_f]
            uPrev               = circshift(uPrev,1)
            uPrev[1,:]          = uCurr[i,:]
            #zCurr[i+1,:]        = simKinModel(zCurr[i,:],uCurr[i,:],modelParams.dt,trackCoeff.coeffCurvature,modelParams)
            zCurr[i+1,:]        = simDynModel(zCurr[i,:],uCurr[i,:],modelParams.dt,modelParams,trackCoeff)

            println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus)")
            # if size(mpcSol.z,2) == 4
            #     println("Prediction error = ",zCurr[i+1,[6,5,4,1]]-mpcSol.z[2,:])
            # else
            #     println("Prediction error = ",zCurr[i+1,:]-mpcSol.z[2,:])
            # end
            # Check if we're crossing the finish line
            if zCurr[i+1,6] >= posInfo.s_target
                oldTraj.idx_end[lapStatus.currentLap] = oldTraj.count[lapStatus.currentLap]
                oldTraj.oldCost[lapStatus.currentLap] = oldTraj.idx_end[lapStatus.currentLap] - oldTraj.idx_start[lapStatus.currentLap]
                println("Reaching finish line at step $(i+1), cost = $(oldTraj.oldCost[lapStatus.currentLap])")
                finished = true
            end

            # Append new states and inputs to old trajectories
            # ===============================================
            oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap],:,lapStatus.currentLap] = zCurr[i,:]
            oldTraj.oldInput[oldTraj.count[lapStatus.currentLap],:,lapStatus.currentLap] = uCurr[i,:]
            oldTraj.count[lapStatus.currentLap] += 1

            # if necessary: append to end of previous lap
            if lapStatus.currentLap > 1 && zCurr[i,6] < 15.0
                oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap-1],:,lapStatus.currentLap-1] = zCurr[i,:]
                oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap-1],6,lapStatus.currentLap-1] += posInfo.s_target
                oldTraj.oldInput[oldTraj.count[lapStatus.currentLap-1],:,lapStatus.currentLap-1] = uCurr[i,:]
                oldTraj.count[lapStatus.currentLap-1] += 1
            end

            #if necessary: append to beginning of next lap
            if zCurr[i,6] > posInfo.s_target - 15.0
                oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap+1],:,lapStatus.currentLap+1] = zCurr[i,:]
                oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap+1],6,lapStatus.currentLap+1] -= posInfo.s_target
                oldTraj.oldInput[oldTraj.count[lapStatus.currentLap+1],:,lapStatus.currentLap+1] = uCurr[i,:]
                oldTraj.count[lapStatus.currentLap+1] += 1
                oldTraj.idx_start[lapStatus.currentLap+1] = oldTraj.count[lapStatus.currentLap+1]
            end

            # Logging
            log_z[i,:,lapStatus.currentLap] = zCurr[i,:]
            log_u[i,:,lapStatus.currentLap] = uCurr[i,:]
            if lapStatus.currentLap > n_pf
                log_ParInt[i,lapStatus.currentLap] = getvalue(mdl.ParInt)
            end
            i = i + 1
        end
        if no_solution_found == 5
            break
        end
        # i = number of steps to *cross* the finish line -> s[i] >= s_target
        z_final             = zCurr[i,:]

        println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i])) s")
        println("Finished Lap Nr. $j with state $(zCurr[i,:])")

        x_xy = transf_s_to_x(s_track,c_track,zCurr[1:i,6],zCurr[1:i,5])
        log_xy[1:i,:,lapStatus.currentLap] = x_xy

        if j>n_pf
            figure(10)
            path_x,xl,xr = s_to_x(s_track,c_track)
            plot(path_x[:,1],path_x[:,2],"b--",xl[:,1],xl[:,2],"b-",xr[:,1],xr[:,2],"b-")
            plot(x_xy[:,1],x_xy[:,2])
            grid("on")
            title("x-y-view")
            axis("equal")

            # figure(4)
            # subplot(211)
            # title("Old Trajectory #1")
            # plot(oldTraj.oldTraj[:,6,1],oldTraj.oldTraj[:,1:5,1])
            # grid("on")
            # legend(["eY","ePsi","v"])
            # subplot(212)
            # title("Old Trajectory #2")
            # plot(oldTraj.oldTraj[:,6,2],oldTraj.oldTraj[:,1:5,2])
            # grid("on")
            # legend(["eY","ePsi","v"])

            # # Print results
            # # --------------------------------
            # figure(2)
            # subplot(211)
            # plot(zCurr[1:i,6],zCurr[1:i,[5,4,1]])
            # title("Real")
            # legend(["eY","ePsi","v"])
            # xlabel("s [m]")
            # grid("on")
            # subplot(212)
            # plot(zCurr[1:i,6],uCurr[1:i,:])
            # legend(["a","d_f"])
            # grid("on")

            # figure(8)
            # plot(zCurr[1:i,6],cost[1:i,1],"r",zCurr[1:i,6],cost[1:i,2],"g",zCurr[1:i,6],cost[1:i,3],"b",zCurr[1:i,6],cost[1:i,4],"y",zCurr[1:i,6],cost[1:i,5],"m",zCurr[1:i,6],cost[1:i,6],"c")
            # grid(1)
            # title("Cost distribution")
            # legend(["z","z_Term","z_Term_const","deriv","control","lane"])
            # println("Press Enter to continue")

            # readline()
        end
    end
    # Save simulation data
    #log_path = "sim_$(Dates.format(now(), "yyyy_mm_dd_HH_MM")).jld"
    log_path = "sim_dyn.jld"
    save(log_path,"t",t,"z",log_z,"u",log_u,"x",log_xy,"totalCost",oldTraj.oldCost)
    println("Saved simulation data.")
end

function prepareTrack(s_track,c_track)
    sz = size(s_track,1)
    s_target = s_track[end]
    s_new = zeros(3*sz)
    c_new = zeros(3*sz)
    s_new[1:sz] = s_track-s_target
    s_new[sz+1:2*sz] = s_track
    s_new[2*sz+1:3*sz] = s_track + s_target
    c_new[1:sz] = c_track
    c_new[sz+1:2*sz] = c_track
    c_new[2*sz+1:3*sz] = c_track
    return s_new, c_new
end
function find_curvature(s_track,c_track,s,trackCoeff)
    sz = size(s_track,1)
    prev = 100                  # 100*0.01 = 1 meter back
    ahea = 200                  # 200*0.01 = 2 meter ahead
    n_tot = prev+ahea+1
    idx_min = indmin((s-s_track).^2)
    idx = idx_min-prev:idx_min+ahea
    intM = zeros(n_tot,trackCoeff.nPolyCurvature+1)
    for i=1:trackCoeff.nPolyCurvature+1
        intM[:,i] = s_track[idx].^(trackCoeff.nPolyCurvature+1-i)
    end
    coeff = intM\c_track[idx]
    return coeff
end

# Sequence of Laps:
# 1st lap:
# Path following, collect data. Actually only the end of the first lap is used for the data of the 2nd lap.
# 2nd lap:
# Path following, append data to first old trajectory and collect further data
# Data of the end of 1st lap is added to data of 2nd lap.
# 3rd lap:
# Start LMPC, use data of previous trajectories for system ID and LMPC