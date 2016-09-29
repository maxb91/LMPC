using JuMP
using Ipopt
using PyPlot

include("helper/classes.jl")
include("helper/functions.jl")
include("helper/coeffConstraintCost.jl")
include("helper/solveMpcProblem.jl")
include("helper/simModel.jl")

#just loads one specified track with distances betwwen points =1 m and returnx x and y coordinatates
include("helper/loadTestMap.jl")
x_track, y_track = loadTestMap()

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
    mdl                         = MpcModel()

    buffersize                  = 700

    z_Init = zeros(4)

    InitializeParameters(mpcParams,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,buffersize)
    InitializeModel(mdl,mpcParams,modelParams,trackCoeff,z_Init)

    # Simulation parameters
    dt                          = modelParams.dt
    t                           = collect(0:dt:40)
    zCurr                       = zeros(length(t),4)
    uCurr                       = zeros(length(t),2)
    cost                        = zeros(length(t),6)

    posInfo.s_start             = 0
    posInfo.s_target            = 6
    trackCoeff.coeffCurvature   = [0.0,0.0,0.0,0.0,0.0]         # polynomial coefficients for curvature approximation (zeros for straight line)
    i = 2
    z_final = zeros(1,4)
    u_final = zeros(1,2)

    for j=1:10
        lapStatus.currentLap = j

        tt          = zeros(length(t),1)
        zCurr       = zeros(length(t)+1,4)
        zCurr[1,4]  = 0.2
        uCurr       = zeros(length(t),2)
        if j>1                                  # if we are in the second or higher lap
            zCurr[1,:] = z_final
            uCurr[1,:] = u_final
            zCurr[1,1] = 0
        end
        cost        = zeros(length(t),6)
        finished    = false

        i = 1
        while i<length(t) && !finished
            if zCurr[i,1] <= 1 # states from old simulation step
                trackCoeff.coeffCurvature[5] = 0.0
            elseif zCurr[i,1] <= 2
                trackCoeff.coeffCurvature[5] = 0.4
            elseif zCurr[i,1] <= 3
                trackCoeff.coeffCurvature[5] = -0.8
            elseif zCurr[i,1] <= 5
                trackCoeff.coeffCurvature[5] = 0.4
            else
                trackCoeff.coeffCurvature[5] = 0.0
            end

            tic()
            posInfo.s   = zCurr[i,1]
            if j > 1
                coeffConstraintCost(oldTraj,mpcCoeff,posInfo,mpcParams)
            end
            tt1 = toq()
            println("coeffConstr: $tt1 s")
            tic()
	    #todo func z curr oaut of xy
            #localizeVehicleCurvAbs(states,x_track,y_track,mpcParams)
            if i > 1
                solveMpcProblem(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr[i,:]',uCurr[i-1,:]')
            else
                solveMpcProblem(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr[i,:]',u_final')
            end
            tt[i]       = toq()
            cost[i,:]   = mpcSol.cost
            uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f]
            #todo sim model with xy cur
            zCurr[i+1,:]  = simModel(zCurr[i,:],uCurr[i,:],modelParams.dt,trackCoeff.coeffCurvature,modelParams)
            println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")
            println("s = $(zCurr[i,1])")
            if zCurr[i+1,1] >= posInfo.s_target
                println("Reached finish line at step $i")
                finished = true
            end
            i = i + 1
            lapStatus.currentIt = i
        end
        i = i-1
	lapStatus.currentIt -= 1
        z_final = zCurr[i,:]
        u_final = uCurr[i,:]
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

        ax1=subplot(311)
        plot(t,zCurr[1:end-1,1],"y",t,zCurr[1:end-1,2],"r",t,zCurr[1:end-1,3],"g",t,zCurr[1:end-1,4],"b")
        grid(1)
        legend(["s","eY","ePsi","v"])
        title("States")
        ax2=subplot(312,sharex=ax1)
        plot(t,uCurr[:,1],"r",t,uCurr[:,2],"g")
        grid(1)
        title("Control input")
        legend(["a","d_f"])
        ax3=subplot(313,sharex=ax1)
        plot(t,cost[:,1],"r",t,cost[:,2],"g",t,cost[:,3],"b",t,cost[:,4],"y",t,cost[:,5],"m",t,cost[:,6],"c")
        grid(1)
        title("Cost distribution")
        legend(["z","z_Term","z_Term_const","deriv","control","lane"])
        println("Press Enter to continue")
        readline()
    end

end
