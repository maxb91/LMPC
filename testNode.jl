using JuMP
using Ipopt
using PyPlot

include("helper/classes.jl")
include("helper/functions.jl")
include("helper/coeffConstraintCost.jl")
include("helper/solveMpcProblem.jl")
include("helper/simModel.jl")

# Load Variables and create Model:
#println("Loading and defining variables...")
#include("helper/createModel.jl")

# Initialize model by solving it once
#println("Initial solve...")
#solve(mdl)

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
    mdl_pF                      = MpcModel()

    buffersize                  = 700

    z_Init = zeros(6)

    InitializeParameters(mpcParams,mpcParams_pF,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,buffersize)
    InitializeModel(mdl,mpcParams,modelParams,trackCoeff,z_Init)
    InitializeModel_pathFollow(mdl_pF,mpcParams_pF,modelParams,trackCoeff,z_Init[1:4])

    # Simulation parameters
    dt                          = modelParams.dt::Float64
    t                           = collect(0:dt:40)
    zCurr                       = zeros(length(t),6)
    uCurr                       = zeros(length(t),2)
    cost                        = zeros(length(t),6)

    posInfo.s_start             = 0
    posInfo.s_target            = 6
    trackCoeff.coeffCurvature   = [0.0,0.0,0.0,0.0,0.0]         # polynomial coefficients for curvature approximation (zeros for straight line)
    i = 2
    z_final = zeros(4)
    u_final = zeros(2)
    z_pf    = zeros(2)

    for j=1:2
        lapStatus.currentLap = j

        tt          = zeros(length(t),1)
        zCurr       = zeros(length(t),6)
        zCurr[1,1]      = 0.2
        uCurr       = zeros(length(t),2)
        if j>1                                  # if we are in the second or higher lap
            zCurr[1,:] = z_final
            uCurr[1,:] = u_final
            zCurr[1,1] = 0
        end
        cost        = zeros(length(t),6)
        finished    = false

        # Start one lap
        # --------------------------------
        i = 2
        while i<length(t) && !finished
            # Define track curvature
            if zCurr[i-1,6] <= 1
                trackCoeff.coeffCurvature[5] = 0.0
            elseif zCurr[i-1,6] <= 3
                trackCoeff.coeffCurvature[5] = 0.4
            elseif zCurr[i-1,6] <= 5
                trackCoeff.coeffCurvature[5] = -0.4
            else
                trackCoeff.coeffCurvature[5] = 0.0
            end
            println("===============")

            # Calculate coefficients for LMPC (if at least in the 2nd lap)
            tic()
            posInfo.s   = zCurr[i-1,6]
            if j > 1
                coeffConstraintCost(oldTraj,mpcCoeff,posInfo,mpcParams)
            end
            tt1 = toq()
            println("coeffConstr: $tt1 s")

            # Calculate optimal inputs (solve MPC problem)
            tic()
            if j == 1                   # if we are in the first lap
                z_pf = [zCurr[i-1,6],zCurr[i-1,5],zCurr[i-1,4],zCurr[i-1,1]]        # use kinematic model and its states
                solveMpcProblem_pathFollow(mdl_pF,mpcSol,mpcParams_pF,trackCoeff,posInfo,modelParams,z_pf,uCurr[i-1,:]')
            else                        # otherwise: use system-ID-model
                if zCurr[i-1,1] == 0
                    warn("xDot = 0")
                    zCurr[i-1,1] = 0.05
                end
                solveMpcProblem(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr[i-1,:]',uCurr[i-1,:]')
            end
            
            tt[i]       = toq()
            cost[i,:]   = mpcSol.cost

            # Simulate the model
            uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f]
            zCurr[i,:]  = simDynModel_exact(zCurr[i-1,:],uCurr[i,:],modelParams.dt,trackCoeff.coeffCurvature,modelParams)

            println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")
            println("s = $(zCurr[i,6])")

            # Check if we've crossed the finish line
            if zCurr[i,6] >= posInfo.s_target
                println("Reached finish line at step $i")
                finished = true
            end
            i = i + 1
            lapStatus.currentIt = i
        end
        z_final = zCurr[i-1,:]
        u_final = uCurr[i-1,:]
        println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i-1])) s")
        println("Finished Lap Nr. $j")

        # Save states in oldTraj:
        # --------------------------------
        saveOldTraj(oldTraj,zCurr,uCurr,lapStatus,buffersize,modelParams.dt)

        # Print results
        # --------------------------------

        # figure()
        # plot(zCurr[:,1],zCurr[:,2],"r",zCurr[:,1],zCurr[:,3],"g",zCurr[:,1],zCurr[:,4],"b")
        # grid(1)
        # legend(["eY","ePsi","v"])
        # title("States over s")

        # figure(1)
        # plot(zCurr[:,6],zCurr[:,1:5])
        # legend(["xDot","yDot","psiDot","ePsi","eY"])
        # xlabel("s [m]")
        # grid(1)

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
        # plot(t,cost[:,1],"r",t,cost[:,2],"g",t,cost[:,3],"b",t,cost[:,4],"y",t,cost[:,5],"m",t,cost[:,6],"c")
        # grid(1)
        # title("Cost distribution")
        # legend(["z","z_Term","z_Term_const","deriv","control","lane"])
        # println("Press Enter to continue")
        # readline()
    end

end
