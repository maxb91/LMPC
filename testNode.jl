using JuMP
using Ipopt
using PyPlot

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

    buffersize                  = 700

    InitializeParameters(mpcParams,mpcParams_pF,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,buffersize)
    mdl    = MpcModel(mpcParams,mpcCoeff,modelParams,trackCoeff)
    mdl_pF = MpcModel_pF(mpcParams_pF,modelParams,trackCoeff)

    # Simulation parameters
    dt                          = modelParams.dt::Float64
    t                           = collect(0:dt:40)          # time vector
    zCurr                       = zeros(length(t),8)        # these are the simulated states
    zCurr_meas                  = zeros(length(t),6)        # these are the measured states (added noise)
    uCurr                       = zeros(length(t),2)        # these are the inputs
    cost                        = zeros(length(t),6)        # these are MPC cost values at each step

    # Logging parameters
    coeff_sysID                 = (zeros(length(t),4),zeros(length(t),4),zeros(length(t),3))      # xDot, yDot, psiDot
    step_diff                   = zeros(length(t),6)        # one-step errors

    posInfo.s_target            = 20

    z_final         = zeros(8)
    z_final_meas    = zeros(6)

    z_pf            = zeros(4)

    uPrev           = zeros(10,2)

    n_pf            = 3             # number of path-following laps

    s_track = 0.1:.1:20               # s from 0 to 20 meters
    c_track = zeros(200)
    c_track[1:10] = 0
    c_track[11:30] = linspace(0,0.1,20)
    c_track[31:50] = linspace(0.1,0,20)
    c_track[51:100] = 0
    c_track[101:150] = linspace(0,0.2,50)
    c_track[151:200] = linspace(0.2,0,50)

    s_track, c_track = prepareTrack(s_track, c_track)

    # Run 10 laps
    for j=1:10
        # Initialize Lap
        lapStatus.currentLap = j

        tt          = zeros(length(t),1)
        zCurr       = zeros(length(t),8)
        zCurr_meas  = zeros(length(t),6)
        uCurr       = zeros(length(t),2)
        zCurr[1,1]  = 0.2

        if j>1                                  # if we are in the second or higher lap
            zCurr[1,:]          = z_final       # use final state as initial state
            zCurr_meas[1,:]     = z_final_meas
            zCurr[1,6]          = zCurr[1,6]%posInfo.s_target   # and make sure that it's before the finish line (0 <= s < s_target)
            zCurr_meas[1,6]     = zCurr_meas[1,6]%posInfo.s_target
        end

        cost        = zeros(length(t),6)
        finished    = false

        # Start one lap
        # --------------------------------
        i = 1
        while i<length(t) && !finished
            println("///////////////////////////////// STARTING ONE ITERATION /////////////////////////////////")
            # Define track curvature
            trackCoeff.coeffCurvature = find_curvature(s_track,c_track,zCurr[i,6],trackCoeff)

            # Calculate coefficients for LMPC (if at least in the 2nd lap)
            posInfo.s   = zCurr_meas[i,6]
            if j > n_pf
                coeffConstraintCost(oldTraj,mpcCoeff,posInfo,mpcParams,lapStatus)
                coeff_sysID[1][i,:] = mpcCoeff.c_Vx
                coeff_sysID[2][i,:] = mpcCoeff.c_Vy
                coeff_sysID[3][i,:] = mpcCoeff.c_Psi
            end

            # Calculate optimal inputs u_i (solve MPC problem)
            tic()
            if j <= n_pf                   # if we are in the first x laps of path following
                z_pf = [zCurr[i,6],zCurr[i,5],zCurr[i,4],norm(zCurr[i,1:2])]        # use kinematic model and its states
                solveMpcProblem_pathFollow(mdl_pF,mpcSol,mpcParams_pF,trackCoeff,posInfo,modelParams,z_pf,uPrev)
            else                        # otherwise: use system-ID-model
                solveMpcProblem(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr[i,1:6]',uPrev)
            end
            
            # Ideas: Faster dynamics in car. Add first order damping system in MPC prediction. Can we use the exact same data twice? Overfitting?
            # ------> Use only one lap for fitting!
            tt[i]       = toq()
            cost[i,:]   = mpcSol.cost

            # Simulate the model -> calculate x_i+1 = f(x_i, u_i)
            uCurr[i,:]          = [mpcSol.a_x mpcSol.d_f]
            uPrev               = circshift(uPrev,1)
            uPrev[1,:]          = uCurr[i,:]
            zCurr[i+1,:]        = simDynModel_exact(zCurr[i,:],uCurr[i,:],modelParams.dt,modelParams,trackCoeff)
            zCurr_meas[i+1,:]   = zCurr[i+1,1:6] + randn(1,6)*diagm([0.01,0.01,0.001,0.001,0.001,0.001])

            if j <= n_pf
                step_diff[i,1:4] = (mpcSol.z[2,:]-[zCurr[i+1,6] zCurr[i+1,5] zCurr[i+1,4] zCurr[i+1,1]]).^2
            else
                step_diff[i,:] = (mpcSol.z[2,:]-zCurr[i+1,1:6]).^2
            end

            println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")

            # Check if we're crossing the finish line
            if zCurr_meas[i+1,6] >= posInfo.s_target
                oldTraj.idx_end[lapStatus.currentLap] = oldTraj.count[lapStatus.currentLap]
                oldTraj.oldCost[lapStatus.currentLap] = oldTraj.idx_end[lapStatus.currentLap] - oldTraj.idx_start[lapStatus.currentLap]
                println("Reaching finish line at step $(i+1), cost = $(oldTraj.oldCost[lapStatus.currentLap])")
                finished = true
            end

            # Append new states and inputs to old trajectories
            # ===============================================
            oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap],:,lapStatus.currentLap] = zCurr_meas[i,:]
            oldTraj.oldInput[oldTraj.count[lapStatus.currentLap],:,lapStatus.currentLap] = uCurr[i,:]
            oldTraj.count[lapStatus.currentLap] += 1

            # if necessary: append to end of previous lap
            if lapStatus.currentLap > 1 && zCurr_meas[i,6] < 5.0
                oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap-1],:,lapStatus.currentLap-1] = zCurr_meas[i,:]
                oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap-1],6,lapStatus.currentLap-1] += posInfo.s_target
                oldTraj.oldInput[oldTraj.count[lapStatus.currentLap-1],:,lapStatus.currentLap-1] = uCurr[i,:]
                oldTraj.count[lapStatus.currentLap-1] += 1
            end

            #if necessary: append to beginning of next lap
            if zCurr_meas[i,6] > posInfo.s_target - 5.0
                oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap+1],:,lapStatus.currentLap+1] = zCurr_meas[i,:]
                oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap+1],6,lapStatus.currentLap+1] -= posInfo.s_target
                oldTraj.oldInput[oldTraj.count[lapStatus.currentLap+1],:,lapStatus.currentLap+1] = uCurr[i,:]
                oldTraj.count[lapStatus.currentLap+1] += 1
                oldTraj.idx_start[lapStatus.currentLap+1] = oldTraj.count[lapStatus.currentLap+1]
            end
            #if j>=3
            #    printPrediction(mpcSol)
            #end
            # if j == 3
            #     figure(1)
            #     title("System ID coefficients")
            #     subplot(311)
            #     plot(zCurr[i-1:i,6],coeff_sysID[1][i-1:i,:],"*-")
            #     legend(["1","2","3"])
            #     title("xDot")
            #     grid(1)
            #     subplot(312)
            #     plot(zCurr[i-1:i,6],coeff_sysID[2][i-1:i,:],"*-")
            #     legend(["1","2","3","4"])
            #     title("yDot")
            #     grid(1)
            #     subplot(313)
            #     plot(zCurr[i-1:i,6],coeff_sysID[3][i-1:i,:],"*-")
            #     legend(["1","2","3"])
            #     title("psiDot")
            #     grid(1)
            #     readline()
            # end

            i = i + 1
            lapStatus.currentIt = i
        end

        # i = number of steps to *cross* the finish line -> s[i] >= s_target
        lapStatus.currentIt = i
        z_final             = zCurr[i,:]
        z_final_meas        = zCurr_meas[i,:]

        println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i])) s")
        println("Finished Lap Nr. $j with state $(zCurr[i,:])")

        if j>n_pf
            figure(4)
            subplot(211)
            title("Old Trajectory #1")
            plot(oldTraj.oldTraj[:,6,1],oldTraj.oldTraj[:,1:5,1])
            grid("on")
            legend(["xDot","yDot","psiDot","ePsi","eY"])
            subplot(212)
            title("Old Trajectory #2")
            plot(oldTraj.oldTraj[:,6,2],oldTraj.oldTraj[:,1:5,2])
            grid("on")
            legend(["xDot","yDot","psiDot","ePsi","eY"])

            # Print results
            # --------------------------------
            figure(5)
            plot(zCurr[1:i,6],step_diff[1:i,:])
            grid("on")
            title("One step errors")
            if j>2
                figure(1)
                title("System ID coefficients")
                ax=subplot(311)
                plot(zCurr[1:i,6],coeff_sysID[1][1:i,:])
                legend(["1","2","3"])
                title("xDot")
                grid(1)
                subplot(312,sharex=ax)
                plot(zCurr[1:i,6],coeff_sysID[2][1:i,:])
                legend(["1","2","3","4"])
                title("yDot")
                grid(1)
                subplot(313,sharex=ax)
                plot(zCurr[1:i,6],coeff_sysID[3][1:i,:])
                legend(["1","2","3"])
                title("psiDot")
                grid(1)
            end

            figure(2)
            subplot(211)
            plot(zCurr[1:i,6],zCurr[1:i,1:5],zCurr[1:i,6],zCurr[1:i,7:8])
            title("Real")
            legend(["xDot","yDot","psiDot","ePsi","eY","a","d_f"])
            xlabel("s [m]")
            grid("on")
            subplot(212)
            plot(zCurr[1:i,6],uCurr[1:i,:])
            legend(["a","d_f"])
            grid("on")

            figure(3)
            plot(zCurr_meas[1:i,6],zCurr_meas[1:i,1:5])
            title("Measured")
            legend(["xDot","yDot","psiDot","ePsi","eY","a","d_f"])
            xlabel("s [m]")
            grid("on")

            # figure(2)
            # ax1=subplot(311)
            # plot(t,zCurr[:,6],"y",t,zCurr[:,5],"r",t,zCurr[:,4],"g",t,zCurr[:,1],"b")
            # grid(1)
            # legend(["s","eY","ePsi","v"])
            # title("States")
            # ax2=subplot(312,sharex=ax1)
            # plot(t,uCurr[:,1],"r",t,uCurr[:,2],"g")
            # grid(1)
            # title("Control input")
            # legend(["a","d_f"])
            # ax3=subplot(313,sharex=ax1)
            figure(8)
            plot(zCurr_meas[1:i,6],cost[1:i,1],"r",zCurr_meas[1:i,6],cost[1:i,2],"g",zCurr_meas[1:i,6],cost[1:i,3],"b",zCurr_meas[1:i,6],cost[1:i,4],"y",zCurr_meas[1:i,6],cost[1:i,5],"m",zCurr_meas[1:i,6],cost[1:i,6],"c")
            grid(1)
            title("Cost distribution")
            legend(["z","z_Term","z_Term_const","deriv","control","lane"])
            println("Press Enter to continue")

            readline()
        end
    end
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
    prev = 10
    ahea = 20                   # 30*0.1 = 3 meter ahead
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
function plot_curvature(s_track,c_track,trackCoeff)
    for i=0:.1:20
        s = i-1:.1:i+1
        coeff = find_curvature(s_track,c_track,i,trackCoeff)
        #c = [s.^8 s.^7 s.^6 s.^5 s.^4 s.^3 s.^2 s.^1 s.^0]*coeff
        c = zeros(size(s,1))
        for i=1:trackCoeff.nPolyCurvature+1
            c += s.^(trackCoeff.nPolyCurvature+1-i)*coeff[i]
        end
        plot(s,c)
    end
    plot(s_track,c_track)
    grid("on")
end
# Sequence of Laps:
# 1st lap:
# Path following, collect data. Actually only the end of the first lap is used for the data of the 2nd lap.
# 2nd lap:
# Path following, append data to first old trajectory and collect further data
# Data of the end of 1st lap is added to data of 2nd lap.
# 3rd lap:
# Start LMPC, use data of previous trajectories for system ID and LMPC