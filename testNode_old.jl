using JuMP
using Ipopt
using PyPlot
using Polynomials

include("helper/status.jl")
include("helper/coeffConstraintCost.jl")
include("helper/solveMpcProblem.jl")
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
mpcCoeff    = coeffConstraintCost(oldTraj,lapStatus,mpcCoeff,posInfo,mpcParams,stateIn,inputIn)

QderivZ = 0*[1 1 1 1]       # cost matrix for derivative cost of states
QderivU = 0.1*[1 1]         # cost matrix for derivative cost of inputs
mpcParams.R = 0.1*[1, 1]        # cost matrix for control inputs

# Simulate System

t           = collect(1:dt:10)
tt          = zeros(length(t),1)
zCurr       = zeros(length(t),4)
zCurr[1,:]  = [0 0.3 0.3 0]
uCurr       = zeros(length(t),2)
cost        = zeros(length(t),5)

for i=2:length(t)
    tic()
    mpcSol      = solveMpcProblem(mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr[i-1,:]',uCurr[i-1,:]')
    tt[i] = toc();
    cost[i,:]   = mpcSol.cost
    uCurr[i,:]  = mpcSol.u[:,1]'
    zCurr[i,:]  = simModel(zCurr[i-1,:],uCurr[i,:],modelParams.dt,trackCoeff.coeffCurvature,modelParams)
    println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")
end
println("=================\nFinished Solving. Avg. time = $(mean(tt)) s")


# Print results
subplot(311)
plot(t,zCurr)
legend(["s","eY","ePsi","v"])
grid()
title("States")
subplot(312)
plot(t,uCurr)
grid()
title("Control input")
legend(["a","d_f"])
subplot(313)
plot(t,cost)
grid()
title("Cost distribution")
legend(["z","z_Term","z_Term_const","deriv","control"])