function saveOldTraj(oldTraj::classes.OldTrajectory,zCurr::Array{Float64}, zCurr_x::Array{Float64},uCurr::Array{Float64},lapStatus::classes.LapStatus,buffersize::Int64,dt::Float64)
    
    i               = lapStatus.currentIt           # current iteration number, just to make notation shorter
    zCurr_export    = zeros(buffersize,4)
    uCurr_export    = zeros(buffersize,2)
    
    println(size(ones(buffersize-i+1,1)))

    #we just take i-1 beacuse the last s value is already the begining ogf the new round and we jsut coun from 0 <= s < s_target
    #why bufferisze because length t is not in ros 
    zCurr_export    = cat(1,zCurr[1:i-1,:], [zCurr[i-1,1]+collect(1:buffersize-i+1)*dt*zCurr[i-1,4] ones(buffersize-i+1,1)*zCurr[i-1,2:4]']) # extrapolate values for after the finish line so that the old trjectory has feasible vlaues for the interpolation in the next round 
    uCurr_export    = cat(1,uCurr[1:i-1,:], zeros(buffersize-i+1,2)) # we need the i-1 as we do not want to keep the last vlaues which is s >= s_target
    costLap         = lapStatus.currentIt               # the cost of the current lap is the time it took to reach the finish line
    
    # Save all data in oldTrajectory:
    if lapStatus.currentLap == 1     
        for i= 1:5                   # if it's the first lap
        oldTraj.oldTraj[:,:,i]  = zCurr_export          # ... just save everything
        oldTraj.oldInput[:,:,i] = uCurr_export
        oldTraj.oldTrajXY[1:size(zCurr_x)[1],:,i] = zCurr_x
        oldTraj.oldCost[i] = costLap
        end
    else
        for i = 4:-1:1
            oldTraj.oldTraj[:,:,i+1]  = oldTraj.oldTraj[:,:,i]    # ... copy the first in the second
            oldTraj.oldInput[:,:,i+1] = oldTraj.oldInput[:,:,i]   # ... same for the input
            oldTraj.oldTrajXY[:,:,i+1]  = oldTraj.oldTrajXY[:,:,i]   
            oldTraj.oldCost[i+1] = oldTraj.oldCost[i]
        end
        oldTraj.oldTraj[:,:,1]  = zCurr_export                 # ... and write the new traj in the first
        oldTraj.oldInput[:,:,1] = uCurr_export
        oldTraj.oldTrajXY[1:size(zCurr_x)[1],:,1] = zCurr_x
        oldTraj.oldCost[1] = costLap
    end
    
end

function InitializeModel(m::classes.MpcModel,mpcParams::classes.MpcParams,modelParams::classes.ModelParams,trackCoeff::classes.TrackCoeff,z_Init::Array{Float64,1}, obstacle::classes.Obstacle)

    dt        = modelParams.dt
    L_a       = modelParams.l_A
    L_b       = modelParams.l_B
    c0        = modelParams.c0
    u_lb      = modelParams.u_lb
    u_ub      = modelParams.u_ub
    z_lb      = modelParams.z_lb
    z_ub      = modelParams.z_ub
    s_obst    = obstacle.s_obstacle
    sy_obst    = obstacle.sy_obstacle

    N         = mpcParams.N

    n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation
    
    m.mdl = Model(solver = IpoptSolver(print_level=0))#,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))

    @variable( m.mdl, m.z_Ol[1:(N+1),1:4])      # z = s, ey, epsi, v
    @variable( m.mdl, m.u_Ol[1:N,1:2])          # overwrtie dim of in classes.jl?
    @variable( m.mdl, 0 <= m.lambda[1:5] <= 1)

    #!!
    #@variable( m.mdl, m.t[1:N+1])

    for i=1:2       # I don't know why but somehow the short method returns errors sometimes
        for j=1:N
            setlowerbound(m.u_Ol[j,i], modelParams.u_lb[j,i])
            setupperbound(m.u_Ol[j,i], modelParams.u_ub[j,i])
        end
    end
    for i=1:4
        for j=1:N+1
            setlowerbound(m.z_Ol[j,i], modelParams.z_lb[j,i])
            setupperbound(m.z_Ol[j,i], modelParams.z_ub[j,i])
        end
    end

    @NLparameter(m.mdl, m.z0[i=1:4] == z_Init[i])
    @NLconstraint(m.mdl, [i=1:4], m.z_Ol[1,i] == m.z0[i])

    @constraint(m.mdl, m.lambda[1]+m.lambda[2]+m.lambda[3]+m.lambda[4]+m.lambda[5]== 1)

   
    #!! object avoidance constraint
    # set s and y obst as nlparamter
    #@NLconstraint(m.mdl, [i=1:N+1], ((m.z_Ol[i,1]-s_obst)/0.3)^2+((m.z_Ol[i,2]-sy_obst)/0.2)^2 == 1+m.t[i]^2)
   #1/m.t

    @NLparameter(m.mdl, m.coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])


    #@NLexpression(m.mdl, m.c[i = 1:N],    m.coeff[1]*m.z_Ol[1,i]^4+m.coeff[2]*m.z_Ol[1,i]^3+m.coeff[3]*m.z_Ol[1,i]^2+m.coeff[4]*m.z_Ol[1,i]+m.coeff[5])
    @NLexpression(m.mdl, m.c[i = 1:N],    sum{m.coeff[j]*m.z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + m.coeff[n_poly_curv+1]) 
    @NLexpression(m.mdl, m.bta[i = 1:N],  atan( L_b / (L_a + L_b) * tan(m.u_Ol[i,2]) ) )
    @NLexpression(m.mdl, m.dsdt[i = 1:N], m.z_Ol[i,4]*cos(m.z_Ol[i,3]+m.bta[i])/(1-m.z_Ol[i,2]*m.c[i]))
    # System dynamics
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,1]  == m.z_Ol[i,1] + dt*m.dsdt[i]  )
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,2]  == m.z_Ol[i,2] + dt*m.z_Ol[i,4]*sin(m.z_Ol[i,3]+m.bta[i])  )                     # ey
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,3]  == m.z_Ol[i,3] + dt*(m.z_Ol[i,4]/L_b*sin(m.bta[i])-m.dsdt[i]*m.c[i])  )            # epsi
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,4]  == m.z_Ol[i,4] + dt*(m.u_Ol[i,1]))#- 0.23*abs(m.z_Ol[i,4]) * m.z_Ol[i,4]))#0.63  # v
    # @NLconstraint(m.mdl, [i=1:N+1], m.z_Ol[i,5]  == s_obst - m.z_Ol[i,1] )
    # @NLconstraint(m.mdl, [i=1:N+1], m.z_Ol[i,6]  == sy_obst - m.z_Ol[i,2] )
    # for i=1:N
       
    #     #@NLconstraint(m.mdl, m.z_Ol[i+1,1]  == m.z_Ol[i,1] + dt*m.dsdt[i]  )                                             # s
    #     @NLconstraint(m.mdl, m.z_Ol[i+1,2]  == m.z_Ol[i,2] + dt*m.z_Ol[i,4]*sin(m.z_Ol[i,3]+m.bta[i])  )                     # ey
    #     @NLconstraint(m.mdl, m.z_Ol[i+1,3]  == m.z_Ol[i,3] + dt*(m.z_Ol[i,4]/L_a*sin(m.bta[i])-m.dsdt[i]*m.c[i])  )            # epsi
    #     @NLconstraint(m.mdl, m.z_Ol[i+1,4]  == m.z_Ol[i,4] + dt*(m.u_Ol[i,1] - 0.03*abs(m.z_Ol[i,4]) * m.z_Ol[i,4]))#0.63  # v
    # end

end

function InitializeParameters(mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,modelParams::classes.ModelParams,
                                posInfo::classes.PosInfo,oldTraj::classes.OldTrajectory,mpcCoeff::classes.MpcCoeff,lapStatus::classes.LapStatus,obstacle::classes.Obstacle,buffersize::Int64, n_rounds::Int64)
    mpcParams.N                 = 10                        #lenght of prediction horizon
    mpcParams.nz                = 4                         #number of States
    mpcParams.Q                 = [0.0,10.0,0.0,10.0]  #0 10 0 1    # put weights on ey, epsi and v, just for first round of PathFollowing
    mpcParams.Q_term            = 100*[1.0,1.0,0.1]           # weights for terminal constraints (LMPC, for e_y, e_psi, and v)
    mpcParams.Q_cost            = 1
    mpcParams.R                 = 0*[1.0,1.0]             # put weights on a and d_f
    mpcParams.QderivZ           = 1.0*[0,0.0,0.1,0.1]             # cost matrix for derivative cost of states
    mpcParams.QderivU           = 0.1*[1,10]               # cost matrix for derivative cost of inputs
    mpcParams.vPathFollowing    = 0.6                 # reference speed for first lap of path following

    trackCoeff.nPolyCurvature   = 4                       # 4th order polynomial for curvature approximation
    trackCoeff.coeffCurvature   = zeros(trackCoeff.nPolyCurvature+1)         # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.width            = 0.6                   # width of the track (0.6m)
    trackCoeff.ds               = 1//10 # is defined as a rational number so we can use it to calculate indices in matrix. with float becomes error

    modelParams.u_lb            = ones(mpcParams.N,1) * [-1.0  -pi/6]                    # lower bounds on steering
    modelParams.u_ub            = ones(mpcParams.N,1) * [ 5.0   pi/6]       #1.2           # upper bounds
    modelParams.z_lb            = ones(mpcParams.N+1,1) * [-Inf -Inf -Inf -0.1]                    # lower bounds on states
    modelParams.z_ub            = ones(mpcParams.N+1,1) * [ Inf  Inf  Inf  2.0]                 # upper bounds
    modelParams.l_A             = 0.125
    modelParams.l_B             = 0.125 #0.125

    modelParams.dt              = 0.1

    posInfo.s_start             = 0.0
    posInfo.s_target            = 5.0

    oldTraj.oldTraj             = zeros(buffersize,4,5)
    # oldTraj.oldTraj[:,1,1]      = 0.1:0.1:buffersize/10 # are these necessarxy?
    # oldTraj.oldTraj[:,1,2]      = 1:buffersize
    oldTraj.oldTrajXY           = zeros(buffersize,4,5)
    oldTraj.oldInput            = zeros(buffersize,2,5)
    oldTraj.oldCost             = 1000*ones(Int64,5)                   # dummies for initialization

    mpcCoeff.order              = 5
    mpcCoeff.coeffCost          = zeros(mpcCoeff.order+1,5)
    mpcCoeff.coeffConst         = zeros(mpcCoeff.order+1,5,3) # nz-1 because no coeff for s
    mpcCoeff.pLength            = 4*mpcParams.N        # small values here may lead to numerical problems since the functions are only approximated in a short horizon

    lapStatus.currentLap        = 1         # initialize lap number
    lapStatus.currentIt         = 0         # current iteration in lap #?? why not 1?

    obstacle.s_obstacle = zeros(buffersize,n_rounds)
    obstacle.sy_obstacle = zeros(buffersize,n_rounds)
    obstacle.rs = 0
    obstacle.ry = 0
    obstacle.index = zeros(buffersize)##is not used at the moment can be deleted in final version
    obstacle.xy_vector = zeros(buffersize,2,n_rounds)
    obstacle.axis_y_up = zeros(buffersize,2,n_rounds)
    obstacle.axis_y_down = zeros(buffersize,2,n_rounds)
    obstacle.axis_s_up = zeros(buffersize,2,n_rounds)
    obstacle.axis_s_down = zeros(buffersize,2,n_rounds)
end