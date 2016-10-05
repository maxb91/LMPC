# Important information about oldTraj:
# ======================================
# oldTraj.oldTraj contains all state information of one lap. It is structured like follows:
# The first <prebuf> values are the end of one lap before the next lap starts (necessary for inter-lap-system-ID)
# The last value of <prebuf> ends with the last value *before* the finish line (s < s_target)
# The s-values within <prebuf> need to be below zero (-> before the finish line). Otherwise they can be mistaken as the end of the current (not previous) lap.
# After <prebuf> follow the recorded states of the trajectory (as many as there were during one lap)
# The number of recorded states during one trajectory (0 <= s < s_target) is equal to <oldCost>.
# After the recorded trajectory, the rest of the vector (until <buffersize>) is filled up with constant values

function saveOldTraj(oldTraj::OldTrajectory,zCurr::Array{Float64},uCurr::Array{Float64},lapStatus::LapStatus,posInfo::PosInfo,buffersize::Int64,dt::Float64)
                
                i               = lapStatus.currentIt-1         # i = number of points for 0 <= s < s_target (= cost of this lap)
                prebuf          = oldTraj.prebuf                # so many points of the end of the previous old traj will be attached to the beginning
                zCurr_export    = zeros(buffersize,6)
                uCurr_export    = zeros(buffersize,2)

                zCurr_export    = cat(1,oldTraj.oldTraj[oldTraj.oldCost[1]+1:oldTraj.oldCost[1]+prebuf,:,1],
                                        zCurr[1:i,:], NaN*ones(buffersize-i-prebuf,6))#[ones(buffersize-i-prebuf,1)*zCurr[i,1:5] zCurr[i,6]+collect(1:buffersize-i-prebuf)*dt*zCurr[i,1]])
                uCurr_export    = cat(1,oldTraj.oldInput[oldTraj.oldCost[1]+1:oldTraj.oldCost[1]+prebuf,:,1],
                                        uCurr[1:i,:], NaN*ones(buffersize-i-prebuf,2))#zeros(buffersize-i-prebuf,2))

                zCurr_export[1:prebuf,6] -= posInfo.s_target       # make the prebuf-values below zero
                costLap                   = i                      # the cost of the current lap is the time it took to reach the finish line

                # Save all data in oldTrajectory:
                if lapStatus.currentLap <= 2                        # if it's the first or second lap
                    oldTraj.oldTraj[:,:,1]  = zCurr_export          # ... just save everything
                    oldTraj.oldInput[:,:,1] = uCurr_export
                    oldTraj.oldTraj[:,:,2]  = zCurr_export
                    oldTraj.oldInput[:,:,2] = uCurr_export
                    oldTraj.oldCost = [costLap,costLap]
                else                                                # idea: always copy the new trajectory in the first array!
                    if oldTraj.oldCost[1] < oldTraj.oldCost[2]      # if the first old traj is better than the second
                        oldTraj.oldTraj[:,:,2]  = oldTraj.oldTraj[:,:,1]    # ... copy the first in the second
                        oldTraj.oldInput[:,:,2] = oldTraj.oldInput[:,:,1]   # ... same for the input
                        oldTraj.oldCost[2] = oldTraj.oldCost[1]
                    end
                    oldTraj.oldTraj[:,:,1]  = zCurr_export                 # ... and write the new traj in the first
                    oldTraj.oldInput[:,:,1] = uCurr_export
                    oldTraj.oldCost[1] = costLap
                end
end

function InitializeParameters(mpcParams::MpcParams,mpcParams_pF::MpcParams,trackCoeff::TrackCoeff,modelParams::ModelParams,
                                posInfo::PosInfo,oldTraj::OldTrajectory,mpcCoeff::MpcCoeff,lapStatus::LapStatus,buffersize::Int64)
    mpcParams.N                 = 5
    mpcParams.Q_term            = 0.1*[0.01,1.0,1.0,1.0,1.0]     # weights for terminal constraints (LMPC, for xDot,yDot,psiDot,ePsi,eY)
    mpcParams.R                 = 0*[1.0,1.0]                   # put weights on a and d_f
    mpcParams.QderivZ           = 0.0*[0,0,0.1,0,0,0]           # cost matrix for derivative cost of states
    mpcParams.QderivU           = 0.1*[1,10]                    # cost matrix for derivative cost of inputs
    mpcParams.Q_term_cost       = 0.01                          # scaling of Q-function

    mpcParams_pF.N              = 5
    mpcParams_pF.Q              = [0.0,10.0,0.1,1.0]
    mpcParams_pF.R              = 0*[1.0,1.0]               # put weights on a and d_f
    mpcParams_pF.QderivZ        = 0.0*[0,0,0.1,0]           # cost matrix for derivative cost of states
    mpcParams_pF.QderivU        = 0.1*[1,1]                 # cost matrix for derivative cost of inputs
    mpcParams_pF.vPathFollowing = 0.5                       # reference speed for first lap of path following

    trackCoeff.nPolyCurvature   = 4                         # 4th order polynomial for curvature approximation
    trackCoeff.coeffCurvature   = zeros(trackCoeff.nPolyCurvature+1)         # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.width            = 0.6                       # width of the track (0.5m)

    modelParams.u_lb            = ones(mpcParams.N,1) * [-0.2  -pi/6]                           # lower bounds on steering
    modelParams.u_ub            = ones(mpcParams.N,1) * [1.2   pi/6]                            # upper bounds
    modelParams.z_lb            = ones(mpcParams.N+1,1)*[-Inf -Inf -Inf -0.5]                   # lower bounds on states
    modelParams.z_ub            = ones(mpcParams.N+1,1)*[ Inf  Inf  Inf  2.0]                   # upper bounds
    modelParams.l_A             = 0.125
    modelParams.l_B             = 0.125
    modelParams.dt              = 0.1
    modelParams.m               = 1.98
    modelParams.I_z             = 0.24
    modelParams.c_f             = 0.63                 # friction coefficient: xDot = - c_f*xDot² (aerodynamic+tire)

    posInfo.s_start             = 0.0
    posInfo.s_target            = 5.0

    oldTraj.oldTraj             = NaN*ones(buffersize,6,2)
    oldTraj.oldInput            = NaN*ones(buffersize,2,2)
    oldTraj.oldCost             = ones(Int64,2)                   # dummies for initialization
    oldTraj.prebuf              = 30
    oldTraj.postbuf             = 30

    mpcCoeff.order              = 5
    mpcCoeff.coeffCost          = zeros(mpcCoeff.order+1,2)
    mpcCoeff.coeffConst         = zeros(mpcCoeff.order+1,2,5)
    mpcCoeff.pLength            = 4*mpcParams.N        # small values here may lead to numerical problems since the functions are only approximated in a short horizon
    mpcCoeff.c_Vx               = zeros(4)
    mpcCoeff.c_Vy               = zeros(4)
    mpcCoeff.c_Psi              = zeros(3)

    lapStatus.currentLap        = 1         # initialize lap number
    lapStatus.currentIt         = 0         # current iteration in lap
end