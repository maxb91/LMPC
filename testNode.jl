using JuMP
using Ipopt
using PyPlot

include("helper/status.jl")
include("helper/coeffConstraintCost.jl")
include("helper/solveMpcProblem.jl")
include("helper/ComputeCostLap.jl")
include("simModel.jl")

# Load Variables and create Model:
println("Loading and defining variables...")
include("createModel.jl")

# Initialize model by solving it once
println("Initial solve...")
solve(mdl)

# Run Model
println("Run 1: ")
# Find Coefficients

function run_sim()

    # DEFINE PARAMETERS
    oldTraj         = OldTrajectory()
    lapStatus       = LapStatus(1,1)
    mpcCoeff        = MpcCoeff()
    posInfo         = PosInfo()
    mpcSol          = MpcSol()

    buffersize                  = 700
    oldTraj.oldTraj             = zeros(buffersize,4,2)
    oldTraj.oldInput            = zeros(buffersize,2,2)

    posInfo.s_start             = 0
    posInfo.s_target            = 2
    mpcCoeff.order              = 5
    mpcCoeff.pLength            = 4*mpcParams.N        # small values here may lead to numerical problems since the functions are only approximated in a short horizon


    mpcParams.QderivZ       = 0.0*[1 1 1 1]     # cost matrix for derivative cost of states
    mpcParams.QderivU       = 0.1*[1 1]         # cost matrix for derivative cost of inputs
    mpcParams.R             = 0.0*[1 1]        # cost matrix for control inputs
    mpcParams.Q             = [0.0 10.0 10.0 1.0]     # put weights on ey, epsi and v

    global mdl, trackCoeff

    # Simulate System
    t           = collect(0:dt:40)
    zCurr       = zeros(length(t),4)
    uCurr       = zeros(length(t),2)
    cost        = zeros(length(t),6)
    posInfo.s_target            = 6

    trackCoeff.coeffCurvature   = [0.0,0.0,0.0,0.0,0.0]         # polynomial coefficients for curvature approximation (zeros for straight line)
    i = 2
    z_final = zeros(1,4)
    u_final = zeros(1,2)

    for j=1:10
        lapStatus.currentLap = j

        tt          = zeros(length(t),1)
        zCurr       = zeros(length(t),4)
        zCurr[1,4]  = 0.2
        uCurr       = zeros(length(t),2)
        if j>1                                  # if we are in the second or higher lap
            zCurr[1,:] = z_final
            uCurr[1,:] = u_final
            zCurr[1,1] = 0
        end
        cost        = zeros(length(t),6)
        finished    = false

        i = 2
        while i<length(t) && !finished
            if zCurr[i-1,1] <= 1
                trackCoeff.coeffCurvature[5] = 0.0
            elseif zCurr[i-1,1] <= 2
                trackCoeff.coeffCurvature[5] = 0.4
            elseif zCurr[i-1,1] <= 3
                trackCoeff.coeffCurvature[5] = -0.8
            elseif zCurr[i-1,1] <= 5
                trackCoeff.coeffCurvature[5] = 0.4
            else
                trackCoeff.coeffCurvature[5] = 0.0
            end

            tic()
            posInfo.s   = zCurr[i-1,1]
            mpcCoeff    = coeffConstraintCost(oldTraj,lapStatus,mpcCoeff,posInfo,mpcParams)
            tt1 = toc()
            println("coeffConstr: $tt1")
            tic()
            mpcSol      = solveMpcProblem(mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr[i-1,:]',uCurr[i-1,:]')
            tt[i]       = toc()
            cost[i,:]   = mpcSol.cost
            uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f]
            zCurr[i,:]  = simModel(zCurr[i-1,:],uCurr[i-1,:],modelParams.dt,trackCoeff.coeffCurvature,modelParams)
            println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")
            println("s = $(zCurr[i,1])")
            if zCurr[i,1] >= posInfo.s_target
                println("Reached finish line at step $i")
                finished = true
            end
            i = i + 1
        end
        z_final = zCurr[i-1,:]
        u_final = uCurr[i-1,:]
        println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i-1])) s")
        println("Finished Lap Nr. $j")


        # Save states in oldTraj:
        # --------------------------------
        zCurr_export = cat(1,zCurr[1:i-1,:], [zCurr[i-1,1]+collect(1:buffersize-i+1)*dt*zCurr[i-1,4] ones(buffersize-i+1,1)*zCurr[i-1,2:4]])
        uCurr_export = cat(1,uCurr[1:i-1,:], zeros(buffersize-i+1,2))
        costLap = computeCostLap(zCurr,posInfo.s_target)
        println("costLap = $costLap")

        if lapStatus.currentLap == 1
            oldTraj.oldTraj[:,:,1]  = zCurr_export
            oldTraj.oldInput[:,:,1] = uCurr_export
            oldTraj.oldTraj[:,:,2]  = zCurr_export
            oldTraj.oldInput[:,:,2] = uCurr_export
            oldTraj.oldCost = [costLap,costLap]
        else
            if oldTraj.oldCost[1] < oldTraj.oldCost[2]      # if the first traj is better than the second
                oldTraj.oldTraj[:,:,2]  = zCurr_export     # ...write the new traj in the second
                oldTraj.oldInput[:,:,2] = uCurr_export
                oldTraj.oldCost[2] = costLap
            else
                oldTraj.oldTraj[:,:,1]  = zCurr_export     # if the second traj is better than the first
                oldTraj.oldInput[:,:,1] = uCurr_export     # ...write the new traj in the first
                oldTraj.oldCost[1] = costLap
            end
        end

        # Print results
        # ax1=subplot(311)
        # plot(t,zCurr[:,1],"y",t,zCurr[:,2],"r",t,zCurr[:,3],"g",t,zCurr[:,4],"b")
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