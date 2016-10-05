using JuMP
using Ipopt
using PyPlot

include("helper/classes.jl")
include("helper/functions.jl")
include("helper/coeffConstraintCost.jl")
include("helper/solveMpcProblem.jl")
include("helper/simModel.jl")
include("helper/printHelper.jl")


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
    lapStatus                   = LapStatus(1,1)
    mpcSol                      = MpcSol()
    trackCoeff                  = TrackCoeff()      # info about track (at current position, approximated)
    modelParams                 = ModelParams()
    mpcParams                   = MpcParams()
    mpcParams_pF                = MpcParams()       # for 1st lap (path following)
    mdl                         = MpcModel()
    mdl_pF                      = MpcModel_pF()

    buffersize                  = 700

    z_Init    = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]      # xDot needs to be > 0
    z_Init_pF = zeros(4)

    InitializeParameters(mpcParams,mpcParams_pF,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,buffersize)
    InitializeModel(mdl,mpcParams,modelParams,trackCoeff,z_Init)
    InitializeModel_pathFollow(mdl_pF,mpcParams_pF,modelParams,trackCoeff,z_Init_pF)

    zInput = zeros(50,6)
    uInput = zeros(50,2)

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


    posInfo.s_start             = 0
    posInfo.s_target            = 7
    trackCoeff.coeffCurvature   = [0.0,0.0,0.0,0.0,0.0]         # polynomial coefficients for curvature approximation (zeros for straight line)
    
    z_final         = zeros(8)
    u_final         = zeros(2)
    z_final_meas    = zeros(6)

    z_pf            = zeros(4)

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

            oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+1,:,1] = zCurr_meas[1,:]
            oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+1,6,1] += posInfo.s_target
            oldTraj.oldInput[oldTraj.oldCost[1]+oldTraj.prebuf+1,:,1] = uCurr[1,:]
            if j==3
                oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+1,:,2] = zCurr_meas[1,:]
                oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+1,6,2] += posInfo.s_target
                oldTraj.oldInput[oldTraj.oldCost[1]+oldTraj.prebuf+1,:,2] = uCurr[1,:]
            end
        end

        cost        = zeros(length(t),6)
        finished    = false

        # Start one lap
        # --------------------------------
        i = 1
        while i<length(t) && !finished
            # Define track curvature
            if zCurr[i,6] <= 1.0
                trackCoeff.coeffCurvature[5] = 0.0
            elseif zCurr[i,6] <= 3
                trackCoeff.coeffCurvature[5] = 0.2
            elseif zCurr[i,6] <= 5
                trackCoeff.coeffCurvature[5] = -0.2
            else
                trackCoeff.coeffCurvature[5] = 0.0
            end
            println("///////////////////////////////// STARTING ONE ITERATION /////////////////////////////////")

            # Calculate coefficients for LMPC (if at least in the 2nd lap)
            tic()
            posInfo.s   = zCurr_meas[i,6]
            if j > 2
                coeffConstraintCost(oldTraj,mpcCoeff,posInfo,mpcParams,zInput,uInput,lapStatus)
                coeff_sysID[1][i,:] = mpcCoeff.c_Vx
                coeff_sysID[2][i,:] = mpcCoeff.c_Vy
                coeff_sysID[3][i,:] = mpcCoeff.c_Psi
            end
            tt1 = toq()
            println("coeffConstr: $tt1 s")

            # Find last inputs u_i-1
            if i>1                      # last u is needed for smooth MPC solutions (for derivative of u)
                last_u = uCurr[i-1,:]
            else
                last_u = u_final
            end

            # Calculate optimal inputs u_i (solve MPC problem)
            tic()
            if j <= 2                   # if we are in the first two laps
                z_pf = [zCurr_meas[i,6],zCurr_meas[i,5],zCurr_meas[i,4],zCurr_meas[i,1]]        # use kinematic model and its states
                solveMpcProblem_pathFollow(mdl_pF,mpcSol,mpcParams_pF,trackCoeff,posInfo,modelParams,z_pf,last_u')
            else                        # otherwise: use system-ID-model
                solveMpcProblem(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_meas[i,:]',last_u')
            end
            
            # Ideas: Faster dynamics in car. Add first order damping system in MPC prediction. Can we use the exact same data twice? Overfitting?
            # ------> Use only one lap for fitting!
            tt[i]       = toq()
            cost[i,:]   = mpcSol.cost

            # Simulate the model -> calculate x_i+1 = f(x_i, u_i)
            uCurr[i,:]          = [mpcSol.a_x mpcSol.d_f]
            zCurr[i+1,:]        = simDynModel_exact(zCurr[i,:],uCurr[i,:],modelParams.dt,trackCoeff.coeffCurvature,modelParams)
            zCurr_meas[i+1,:]   = zCurr[i+1,1:6] + randn(1,6)*diagm([0.01,0.01,0.001,0.001,0.001,0.001])

            zInput = circshift(zInput,-1)
            uInput = circshift(uInput,-1)
            zInput[end,:] = zCurr_meas[i,:]
            uInput[end,:] = uCurr[i,:]

            if j <= 2
                step_diff[i,1:4] = (mpcSol.z[2,:]-[zCurr[i+1,6] zCurr[i+1,5] zCurr[i+1,4] zCurr[i+1,1]]).^2
            else
                step_diff[i,:] = (mpcSol.z[2,:]-zCurr[i+1,1:6]).^2
            end

            println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")
            println("s_meas[i+1] = $(zCurr_meas[i+1,6])")
            println("z_meas[i+1] = $(zCurr_meas[i+1,:])")

            # Check if we're crossing the finish line
            if zCurr_meas[i+1,6] >= posInfo.s_target
                println("Reaching finish line at step $(i+1)")
                finished = true
            end

            # append new states and inputs to old trajectory
            oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+i,:,1] = zCurr_meas[i,:]
            oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+i,6,1] += posInfo.s_target
            oldTraj.oldInput[oldTraj.oldCost[1]+oldTraj.prebuf+i,:,1] = uCurr[i,:]
            if j==3     # if its the third lap, append to both old trajectories! (since both are the same)
                oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+i,:,2] = zCurr_meas[i,:]
                oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+i,6,2] += posInfo.s_target
                oldTraj.oldInput[oldTraj.oldCost[1]+oldTraj.prebuf+i,:,2] = uCurr[i,:]
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
        u_final             = uCurr[i-1,:]      # only needed for very first iteration in next lap

        println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i])) s")
        println("Finished Lap Nr. $j with state $(zCurr[i,:])")

        

        # Save states in oldTraj:
        # --------------------------------
        saveOldTraj(oldTraj,zCurr_meas,uCurr,lapStatus,posInfo,buffersize,modelParams.dt)

        if j>2
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
            #println("Old Trajectory:")
            #println(oldTraj.oldTraj[:,:,1])
            #println("Old Input:")
            #println(oldTraj.oldInput[:,:,1])
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

# Sequence of Laps:
# 1st lap:
# Path following, collect data. Actually only the end of the first lap is used for the data of the 2nd lap.
# End of lap: Save trajectory
# 2nd lap:
# Path following, append data to first old trajectory and collect further data
# Data of the end of 1st lap is added to data of 2nd lap.
# 3rd lap:
# Start LMPC, use data of previous trajectories for system ID and LMPC