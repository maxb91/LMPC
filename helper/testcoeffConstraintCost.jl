
include("status.jl")
include("coeffConstraintCost.jl")
# Define Variables
oldTraj         = OldTrajectory()
lapStatus       = LapStatus(2,1)
mpcCoeff        = MpcCoeff()
posInfo         = PosInfo()
mpcParams       = MpcParams()
stateIn         = zeros(7,1)        # stateIn: xdot, ydot, psidot, epsi, y, s, a_f
inputIn         = zeros(2,1)

# ===============================
# Initialize Parameters
# ===============================
t = 1:500
oldTraj.oldTraj         = rand(4,500,2)
oldTraj.oldTraj[1,:,1]  = t*0.9
oldTraj.oldTraj[1,:,2]  = t*0.8
oldTraj.oldTraj[2,:,1]  = cos(t*0.25)
oldTraj.oldTraj[3,:,1]  = sin(t*0.25)
oldTraj.oldTraj[4,:,1]  = cos(t*0.25)
oldTraj.oldInput        = rand(2,500,2)

mpcParams.N             = 50
mpcParams.nz            = 4
mpcParams.OrderCostCons = 4
mpcParams.R             = [100, 0];
mpcParams.Q             = [0 1 1 1]
mpcParams.vPathFollowing = 0.2

mpcCoeff.currentTime    = 0
mpcCoeff.currentTraj    = zeros(6,100)
mpcCoeff.currentInput   = zeros(2,100)
mpcCoeff.coeffCost      = 0
mpcCoeff.coeffConst     = 0


mpcCoeff = coeffConstraintCost(oldTraj, lapStatus, mpcCoeff, posInfo, mpcParams, stateIn, inputIn)

println("Finished.")