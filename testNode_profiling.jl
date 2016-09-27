using JuMP
using Ipopt
using PyPlot
using Polynomials

include("helper/status.jl")
include("helper/coeffConstraintCost.jl")
include("helper/solveMpcProblem.jl")
include("helper/ComputeCostLap.jl")
include("helper/simModel.jl")

# Load Variables and create Model:
println("Loading and defining variables...")
include("helper/createModel.jl")

# Initialize model by solving it once
println("Initial solve...")
solve(mdl)

# Run Model
println("Run 1: ")
# Find Coefficients

QderivZ     = 0*[1 1 1 1]       # cost matrix for derivative cost of states
QderivU     = 0.1*[1 1]         # cost matrix for derivative cost of inputs
mpcParams.R = 0.1*[1 1]         # cost matrix for control inputs

# Simulate System
zCurr       = zeros(length(t),4)
uCurr       = zeros(length(t),2)

trackCoeff.coeffCurvature   = [0,0,0,0,0]         # polynomial coefficients for curvature approximation (zeros for straight line)

for j=1:3
    lapStatus.currentLap = j
    println(j)

    t           = collect(1:dt:20)
    tt          = zeros(length(t),1)
    zCurr       = zeros(length(t),4)
    zCurr[1,:]  = [0 0 0 0]
    uCurr       = zeros(length(t),2)
    #cost        = zeros(length(t),5)

    for i=2:length(t)
        #tic()
        posInfo.s   = zCurr[i-1,1]
        mpcCoeff    = coeffConstraintCost(oldTraj,lapStatus,mpcCoeff,posInfo,mpcParams)
        mpcSol      = solveMpcProblem(mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr[i-1,:]',uCurr[i-1,:]')
        #tt[i]       = toc()
        #cost[i,:]   = mpcSol.cost
        uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f]#mpcSol.u[:,1]'
        zCurr[i,:]  = simModel(zCurr[i-1,:],uCurr[i,:],modelParams.dt,trackCoeff.coeffCurvature,modelParams)
      #  println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")
    end
    #println("=================\nFinished Solving. Avg. time = $(mean(tt)) s")
    #println("Finished Lap Nr. $j")


    # Save states in oldTraj:
    # --------------------------------
    zCurr_export = cat(1,zCurr, [zCurr[end,1]+collect(1:buffersize-length(t))*dt*zCurr[end,4] ones(buffersize-length(t),1)*zCurr[end,2:4]])
    uCurr_export = cat(1,uCurr, zeros(buffersize-length(t),2))
    costLap = computeCostLap(zCurr,posInfo.s_target)
    #println("costLap = $costLap")

    if lapStatus.currentLap == 1
        oldTraj.oldTraj[:,:,1]  = zCurr_export'
        oldTraj.oldInput[:,:,1] = uCurr_export'
        oldTraj.oldTraj[:,:,2]  = zCurr_export'
        oldTraj.oldInput[:,:,2] = uCurr_export'
        oldTraj.oldCost = [costLap,costLap]
    else
        if oldTraj.oldCost[1] < oldTraj.oldCost[2]      # if the first traj is better than the second
            oldTraj.oldTraj[:,:,2]  = zCurr_export'     # ...write the new traj in the second
            oldTraj.oldInput[:,:,2] = uCurr_export'
            oldTraj.oldCost[2] = costLap
        else
            oldTraj.oldTraj[:,:,1]  = zCurr_export'     # if the second traj is better than the first
            oldTraj.oldInput[:,:,1] = uCurr_export'     # ...write the new traj in the first
            oldTraj.oldCost[1] = costLap
        end
    end
end