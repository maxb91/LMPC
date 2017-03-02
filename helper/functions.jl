function saveOldTraj(oldTraj,zCurr::Array{Float64}, zCurr_x::Array{Float64},uCurr::Array{Float64},lapStatus::classes.LapStatus,buffersize::Int64,dt::Float64, 
    load_safeset::Bool, costs::Array{Float64},  lambda_log::Array{Float64},z_pred_log::Array{Float64},u_pred_log::Array{Float64},ssInfOn_log::Array{Float64},
     mpcSol::classes.MpcSol, obstacle, distance2obst::Array{Float64}, curvature_curr::Array{Float64}, copyInfo::Array{Float64,2})
    #!!::classes.OldTrajectory
    i               = lapStatus.currentIt           # current iteration number, just to make notation shorter
    costLap         = lapStatus.currentIt
    zCurr_export    = zeros(buffersize,4)
    uCurr_export    = zeros(buffersize,2)
    
    println(size(ones(buffersize-i+1,1)))

    #we just take i-1 beacuse the last s value is already the begining ogf the new round and we jsut coun from 0 <= s < s_target
    zCurr_export    = cat(1,zCurr[1:i-1,:], [zCurr[i-1,1]+collect(1:buffersize-i+1)*dt*zCurr[i-1,4] ones(buffersize-i+1,1)*zCurr[i-1,2:4]']) # extrapolate values for after the finish line so that the old trjectory has feasible vlaues for the interpolation in the next round 
    uCurr_export    = cat(1,uCurr[1:i-1,:], zeros(buffersize-i+1,2)) # we need the i-1 as we do not want to keep the last vlaues which is s >= s_target
                   # the cost of the current lap is the time it took to reach the finish line
    
    #the cost at each time step is the weighted current iterations to the target plus the costs of the obstacle, inputs and drivatives, not soft constraints
    inpCost = zeros(buffersize)
    derivStateCost = zeros(buffersize)
    currCostObst = zeros(buffersize)
    derivInpCost = zeros(buffersize)
    cost2target = zeros(buffersize)
    
    #save the terminal cost
    for j = 1:buffersize
        cost2target[j] = mpcParams.Q_cost*(costLap-j+1)
    end

    # Save all data in oldTrajectory:
    if lapStatus.currentLap == 1 && load_safeset == false    # if it's the first lap
        for k= 1:oldTraj.n_oldTraj                 
            oldTraj.oldTraj[:,:,k]  = zCurr_export          # ... just save everything
            oldTraj.oldInput[:,:,k] = uCurr_export
            oldTraj.oldTrajXY[1:size(zCurr_x)[1],:,k] = zCurr_x
            oldTraj.oldNIter[k] = costLap
            oldTraj.costs[:,:,k] = costs
            oldTraj.lambda_sol[:,:,k] = lambda_log
            oldTraj.z_pred_sol[:,:,:,k] = z_pred_log
            oldTraj.u_pred_sol[:,:,:,k] = u_pred_log
            oldTraj.ssInfOn_sol[:,:,k]= ssInfOn_log
            oldTraj.eps[:,:,k] = mpcSol.eps
            oldTraj.cost2Target[:,k] = cost2target
            oldTraj.distance2obst[:,k,:] = distance2obst
            oldTraj.curvature[:,k] = curvature_curr
            oldTraj.copyInfo[:,:,k] = copyInfo

            obstacle.s_obstacle[:,k] = obstacle.s_obstacle[:,1]# in the else part we dont shift s_obscle because we do that in the beginning of a new round
            obstacle.sy_obstacle[:,k] = obstacle.sy_obstacle[:,1]
            obstacle.v[:,k] = obstacle.v[:,1]
        end
    else
        for k = oldTraj.n_oldTraj-3:-1:1
            oldTraj.oldTraj[:,:,k+1]  = oldTraj.oldTraj[:,:,k]    # ... copy the first in the second
            oldTraj.oldInput[:,:,k+1] = oldTraj.oldInput[:,:,k]   # ... same for the input
            oldTraj.oldTrajXY[:,:,k+1]  = oldTraj.oldTrajXY[:,:,k]   
            oldTraj.oldNIter[k+1] = oldTraj.oldNIter[k]
            oldTraj.costs[:,:,k+1] = oldTraj.costs[:,:,k]
            oldTraj.lambda_sol[:,:,k+1] = oldTraj.lambda_sol[:,:,k]
            oldTraj.z_pred_sol[:,:,:,k+1] = oldTraj.z_pred_sol[:,:,:,k]
            oldTraj.u_pred_sol[:,:,:,k+1] = oldTraj.u_pred_sol[:,:,:,k]
            oldTraj.ssInfOn_sol[:,:,k+1]= oldTraj.ssInfOn_sol[:,:,k]
            oldTraj.eps[:,:,k+1] = oldTraj.eps[:,:,k]
            oldTraj.cost2Target[:,k+1] = oldTraj.cost2Target[:,k]
            oldTraj.distance2obst[:,k+1,:] = oldTraj.distance2obst[:,k,:]
            oldTraj.curvature[:,k+1] = oldTraj.curvature[:,k]
            oldTraj.copyInfo[:,:,k+1] = oldTraj.copyInfo[:,:,k]
        end
        if costLap > 40
            oldTraj.oldTraj[:,:,1]  = zCurr_export                 # ... and write the new traj in the first
            oldTraj.oldInput[:,:,1] = uCurr_export
            oldTraj.oldTrajXY[1:size(zCurr_x)[1],:,1] = zCurr_x
            oldTraj.oldNIter[1] = costLap
            
            oldTraj.costs[:,:,1] = costs
            oldTraj.lambda_sol[:,:,1] = lambda_log
            oldTraj.z_pred_sol[:,:,:,1] = z_pred_log
            oldTraj.u_pred_sol[:,:,:,1] = u_pred_log
            oldTraj.ssInfOn_sol[:,:,1]= ssInfOn_log
            oldTraj.eps[:,:,1] = mpcSol.eps
            oldTraj.cost2Target[:,1] = cost2target
            oldTraj.distance2obst[:,1,:] = distance2obst
            oldTraj.curvature[:,1] = curvature_curr
            oldTraj.copyInfo[:,:,1] = copyInfo
        else # if cost are unrealisticaly low dont copy round
            warn("trajectory apparently erroneous, copy path folowing instead")
            oldTraj.oldTraj[:,:,1]  = oldTraj.oldTraj[:,:,18]    # ... copy the first in the second
            oldTraj.oldInput[:,:,1] = oldTraj.oldInput[:,:,18]   # ... same for the input
            oldTraj.oldTrajXY[:,:,1]  = oldTraj.oldTrajXY[:,:,18]   
            oldTraj.oldNIter[1] = oldTraj.oldNIter[18]
            oldTraj.costs[:,:,1] = oldTraj.costs[:,:,18]
            oldTraj.lambda_sol[:,:,1] = oldTraj.lambda_sol[:,:,18]
            oldTraj.z_pred_sol[:,:,:,1] = oldTraj.z_pred_sol[:,:,:,18]
            oldTraj.u_pred_sol[:,:,:,1] = oldTraj.u_pred_sol[:,:,:,18]
            oldTraj.ssInfOn_sol[:,:,1]= oldTraj.ssInfOn_sol[:,:,18]
            oldTraj.eps[:,:,1] = oldTraj.eps[:,:,18]
            oldTraj.cost2Target[:,1] = oldTraj.cost2Target[:,18]
            oldTraj.distance2obst[:,1,:] = oldTraj.distance2obst[:,18,:]
            oldTraj.curvature[:,1] = oldTraj.curvature[:,18]
            oldTraj.copyInfo[:,:,1] = oldTraj.copyInfo[:,:,18]
        end

    end
    
end



function InitializeParameters(mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,modelParams::classes.ModelParams,
                                oldTraj::classes.OldTrajectory,mpcCoeff::classes.MpcCoeff,mpcSol::classes.MpcSol,lapStatus::classes.LapStatus,obstacle::classes.Obstacle,buffersize::Int64)
    mpcParams.N                 = 10                        #lenght of prediction horizon
    mpcParams.nz                = 4                         #number of States
    mpcParams.Q                 = [0.0,10.0,0.1,10.0]  #0 10 0 1    # put weights on ey, epsi and v, just for first round of PathFollowing
    mpcParams.Q_term            = 100*[10.0,2.0,1.0]           # weights for terminal constraints (LMPC, for e_y, e_psi, and v)
    mpcParams.Q_cost            = 0.7                           #factor for terminal cost
    mpcParams.Q_obstacle        = 0.3 #
    mpcParams.Q_obstacleNumer   = 0.004#0.04#0.0025#0.0019
    mpcParams.Q_lane            = 2000.0
    mpcParams.Q_velocity        = 1000.0
    mpcParams.R                 = 0.0*[1.0,1.0]             # put weights on a and d_f
    mpcParams.QderivZ           = 0.0*[0,0.0,0.1,0.1]             # cost matrix for derivative cost of states
    mpcParams.QderivU           = 0.1*[1.0,10]               # cost matrix for derivative cost of inputs
    mpcParams.vPathFollowing    = 0.6                 # reference speed for first lap of path following

    trackCoeff.nPolyCurvature   = 4                       # 4th order polynomial for curvature approximation
    trackCoeff.coeffCurvature   = zeros(trackCoeff.nPolyCurvature+1)         # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.width            = 0.6                   # width of the track (0.6m)
    trackCoeff.ds               = 1//10#4//100#1//10 # is defined as a rational number so we can use it to calculate indices in matrix. with float becomes error

    modelParams.u_lb            = ones(mpcParams.N,1) * [-1.0  -pi/6]  #-0.6 for braking                  # lower bou9nds on steering
    modelParams.u_ub            = ones(mpcParams.N,1) * [ 2.5   pi/6]       #1.2           # upper bounds
    modelParams.z_lb            = ones(mpcParams.N+1,1) * [-Inf -Inf -Inf  0.1]                    # lower bounds on states
    modelParams.z_ub            = ones(mpcParams.N+1,1) * [ Inf  Inf  Inf  Inf]               # upper bounds
    modelParams.l_A             = 0.125
    modelParams.l_B             = 0.125 #0.125
    modelParams.v_max           = 2.0
    modelParams.max_alpha       = 12
    modelParams.mass            = 1.98 # kg
    modelParams.mu              = 0.85
    modelParams.g               = 9.81 # m/s^2
    modelParams.I_z             = 0.03 # kg * m^2
    modelParams.B               = 6.0#1.0
    modelParams.C               = 1.6#1.25
    modelParams.dt              = 0.1#0.1

    oldTraj.n_oldTraj     = 20 #number of old Trajectories for safe set
    oldTraj.oldTraj             = zeros(buffersize,4,oldTraj.n_oldTraj)
    oldTraj.oldTrajXY           = zeros(buffersize,6,oldTraj.n_oldTraj)
    oldTraj.distance2obst       = zeros(buffersize,oldTraj.n_oldTraj, obstacle.n_obstacle)
    oldTraj.curvature           = zeros(buffersize,oldTraj.n_oldTraj)
    oldTraj.oldInput            = zeros(buffersize,2,oldTraj.n_oldTraj)
    oldTraj.oldNIter             = 1000*ones(Int64,oldTraj.n_oldTraj)    # dummies for initialization
    oldTraj.costs = zeros(8,buffersize,oldTraj.n_oldTraj)
    oldTraj.lambda_sol = zeros(oldTraj.n_oldTraj,buffersize,oldTraj.n_oldTraj)
    oldTraj.z_pred_sol =zeros(mpcParams.N+1,4,buffersize,oldTraj.n_oldTraj)
    oldTraj.u_pred_sol = zeros(mpcParams.N,2,buffersize,oldTraj.n_oldTraj)
    oldTraj.ssInfOn_sol = zeros(oldTraj.n_oldTraj,buffersize,oldTraj.n_oldTraj)
    oldTraj.eps = zeros(3,buffersize,oldTraj.n_oldTraj)
    oldTraj.cost2Target = zeros(buffersize,oldTraj.n_oldTraj)
    oldTraj.copyInfo = zeros(buffersize,4, oldTraj.n_oldTraj)

    mpcSol.u  = zeros(mpcParams.N,2)
    mpcSol.z  = zeros(mpcParams.N+1,4) 
    mpcSol.lambda = zeros(oldTraj.n_oldTraj)
    mpcSol.lambda[1] = 1
    mpcSol.ssInfOn = ones(oldTraj.n_oldTraj)
    mpcSol.eps = zeros(3,buffersize)

    mpcCoeff.order              = 5
    mpcCoeff.coeffCost          = zeros(buffersize,mpcCoeff.order+1,oldTraj.n_oldTraj)
    mpcCoeff.coeffConst         = zeros(buffersize,mpcCoeff.order+1,oldTraj.n_oldTraj,3) # nz-1 because no coeff for s
    mpcCoeff.pLength            = 4*mpcParams.N        # small values here may lead to numerical problems since the functions are only approximated in a short horizon

    lapStatus.currentLap        = 1         # initialize lap number
    lapStatus.currentIt         = 0         # current iteration in lap #?? why not 1?

    obstacle.s_obstacle = zeros(buffersize,oldTraj.n_oldTraj,obstacle.n_obstacle)
    obstacle.sy_obstacle = zeros(buffersize,oldTraj.n_oldTraj,obstacle.n_obstacle)
    obstacle.v = zeros(buffersize, oldTraj.n_oldTraj,obstacle.n_obstacle)
    obstacle.rs = 0.5
    obstacle.ry = 0.21
    obstacle.xy_vector   = zeros(buffersize,2,oldTraj.n_oldTraj,obstacle.n_obstacle)
    obstacle.axis_y_up   = zeros(buffersize,2,oldTraj.n_oldTraj,obstacle.n_obstacle)
    obstacle.axis_y_down = zeros(buffersize,2,oldTraj.n_oldTraj,obstacle.n_obstacle)
    obstacle.axis_s_up   = zeros(buffersize,2,oldTraj.n_oldTraj,obstacle.n_obstacle)
    obstacle.axis_s_down = zeros(buffersize,2,oldTraj.n_oldTraj,obstacle.n_obstacle)
end
