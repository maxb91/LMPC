function saveOldTraj(oldTraj::OldTrajectory,zCurr::Array{Float64}, zCurr_x::Array{Float64},uCurr::Array{Float64},lapStatus::LapStatus,buffersize::Int64,dt::Float64)
    
    i               = lapStatus.currentIt           # current iteration number, just to make notation shorter
    zCurr_export    = zeros(buffersize,4)
    uCurr_export    = zeros(buffersize,2)
    
    println(size(ones(buffersize-i+1,1)))
    println(zCurr[i-1,2:4])

    #we just take i-1 beacuse the last s value is already the begining ogf the new round and we jsut coun from 0 <= s < s_target
    #why bufferisze because length t is not in ros 
    #why do we need the extrapolated values to use them for the mpc if we are close to the finish line in the foloowing round
    zCurr_export    = cat(1,zCurr[1:i-1,:], [zCurr[i-1,1]+collect(1:buffersize-i+1)*dt*zCurr[i-1,4] ones(buffersize-i+1,1)*zCurr[i-1,2:4]']) #?? werte nach zeillinie wofür? doppelt beschriebB?
    uCurr_export    = cat(1,uCurr[1:i-1,:], zeros(buffersize-i+1,2)) #?? do we need the i-1 we alredy counted it down outisde the loop should just be i ?
    costLap         = lapStatus.currentIt               # the cost of the current lap is the time it took to reach the finish line
    
    # Save all data in oldTrajectory:
    if lapStatus.currentLap == 1                        # if it's the first lap
        oldTraj.oldTraj[:,:,1]  = zCurr_export          # ... just save everything
        oldTraj.oldInput[:,:,1] = uCurr_export
        oldTraj.oldTrajXY[1:size(zCurr_x)[1],:,1] = zCurr_x
        oldTraj.oldTraj[:,:,2]  = zCurr_export
        oldTraj.oldInput[:,:,2] = uCurr_export
        oldTraj.oldTrajXY[1:size(zCurr_x)[1],:,2] = zCurr_x
        oldTraj.oldCost = [costLap,costLap]
    else                                                # idea: always copy the new trajectory in the first array!
        if oldTraj.oldCost[1] < oldTraj.oldCost[2]      # if the first old traj is better than the second
            oldTraj.oldTraj[:,:,2]  = oldTraj.oldTraj[:,:,1]    # ... copy the first in the second
            oldTraj.oldInput[:,:,2] = oldTraj.oldInput[:,:,1]   # ... same for the input
            oldTraj.oldTrajXY[:,:,2]  = oldTraj.oldTrajXY[:,:,1]   
            oldTraj.oldCost[2] = oldTraj.oldCost[1]
        end
        oldTraj.oldTraj[:,:,1]  = zCurr_export                 # ... and write the new traj in the first
        oldTraj.oldInput[:,:,1] = uCurr_export
        oldTraj.oldTrajXY[1:size(zCurr_x)[1],:,1] = zCurr_x
        oldTraj.oldCost[1] = costLap
    end
    
end

function InitializeModel(m::MpcModel,mpcParams::MpcParams,modelParams::ModelParams,trackCoeff::TrackCoeff,z_Init::Array{Float64,1})

    dt   = modelParams.dt
    L_a  = modelParams.l_A
    L_b  = modelParams.l_B
    c0   = modelParams.c0
    u_lb = modelParams.u_lb
    u_ub = modelParams.u_ub
    z_lb = modelParams.z_lb
    z_ub = modelParams.z_ub

    N    = mpcParams.N

    n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation
    
    m.mdl = Model(solver = IpoptSolver(print_level=3,max_cpu_time=0.5))#,linear_solver="ma57",print_user_options="yes"))

    @variable( m.mdl, m.z_Ol[1:(N+1),1:4])      # z = s, ey, epsi, v
    @variable( m.mdl, m.u_Ol[1:N,1:2])          #?? different dim then in classes.jl?
    @variable( m.mdl, 0 <= m.ParInt[1:1] <= 1)

    for i=1:2       # I don't know why but somehow the short method returns errors sometimes #?? short method
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
    #@variable( m.mdl, 1 >= m.ParInt >= 0 ) #?? better like above with [1:1]?

    @NLparameter(m.mdl, m.z0[i=1:4] == z_Init[i])
    @NLconstraint(m.mdl, [i=1:4], m.z_Ol[1,i] == m.z0[i])

    @NLparameter(m.mdl, m.coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])
    #@NLparameter(m.mdl, m.s_startC == 0.1)

    #@NLexpression(m.mdl, m.c[i = 1:N],    m.coeff[1]*m.z_Ol[1,i]^4+m.coeff[2]*m.z_Ol[1,i]^3+m.coeff[3]*m.z_Ol[1,i]^2+m.coeff[4]*m.z_Ol[1,i]+m.coeff[5])
    @NLexpression(m.mdl, m.c[i = 1:N],    sum{m.coeff[j]*m.z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + m.coeff[n_poly_curv+1]) #??  last value x^0  poly to appxo cuvature syntax? what is this,approx of states?
    @NLexpression(m.mdl, m.bta[i = 1:N],  atan( L_a / (L_a + L_b) * ( m.u_Ol[i,2] ) ) )
    @NLexpression(m.mdl, m.dsdt[i = 1:N], m.z_Ol[i,4]*cos(m.z_Ol[i,3]+m.bta[i])/(1-m.z_Ol[i,2]*m.c[i]))
    # System dynamics
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,1]  == m.z_Ol[i,1] + dt*m.dsdt[i]  )
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,2]  == m.z_Ol[i,2] + dt*m.z_Ol[i,4]*sin(m.z_Ol[i,3]+m.bta[i])  )                     # ey
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,3]  == m.z_Ol[i,3] + dt*(m.z_Ol[i,4]/L_a*sin(m.bta[i])-m.dsdt[i]*m.c[i])  )            # epsi
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,4]  == m.z_Ol[i,4] + dt*(m.u_Ol[i,1]))# - 0.63*abs(m.z_Ol[i,4]) * m.z_Ol[i,4]))#0.63  # v
    # for i=1:N
       
    #     #@NLconstraint(m.mdl, m.z_Ol[i+1,1]  == m.z_Ol[i,1] + dt*m.dsdt[i]  )                                             # s
    #     @NLconstraint(m.mdl, m.z_Ol[i+1,2]  == m.z_Ol[i,2] + dt*m.z_Ol[i,4]*sin(m.z_Ol[i,3]+m.bta[i])  )                     # ey
    #     @NLconstraint(m.mdl, m.z_Ol[i+1,3]  == m.z_Ol[i,3] + dt*(m.z_Ol[i,4]/L_a*sin(m.bta[i])-m.dsdt[i]*m.c[i])  )            # epsi
    #     @NLconstraint(m.mdl, m.z_Ol[i+1,4]  == m.z_Ol[i,4] + dt*(m.u_Ol[i,1] - 0.03*abs(m.z_Ol[i,4]) * m.z_Ol[i,4]))#0.63  # v
    # end

end

function InitializeParameters(mpcParams::MpcParams,trackCoeff::TrackCoeff,modelParams::ModelParams,
                                posInfo::PosInfo,oldTraj::OldTrajectory,mpcCoeff::MpcCoeff,lapStatus::LapStatus,buffersize::Int64)
    mpcParams.N                 = 10                        #lenght of prediction horizon
    mpcParams.nz                = 4                         #number of States
    mpcParams.Q                 = [0.0,10.0,1.0,1.0]  #0 10 0 1    # put weights on ey, epsi and v, just for first round of PathFollowing
    mpcParams.Q_term            = 100*[1.0,1.0,0.1]           # weights for terminal constraints (LMPC, for e_y, e_psi, and v)
    mpcParams.Q_cost            = 0.04
    mpcParams.R                 = 0*[1.0,1.0]             # put weights on a and d_f
    mpcParams.QderivZ           = 1.0*[0,0.0,0.1,0.1]             # cost matrix for derivative cost of states
    mpcParams.QderivU           = 0.1*[1,10]               # cost matrix for derivative cost of inputs
    mpcParams.vPathFollowing    = 0.6#!!was 0.6                   # reference speed for first lap of path following

    trackCoeff.nPolyCurvature   = 4                       # 4th order polynomial for curvature approximation
    trackCoeff.coeffCurvature   = zeros(trackCoeff.nPolyCurvature+1)         # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.width            = 0.6                   # width of the track (0.6m)
    trackCoeff.ds               =1//10 # is defined as a rational number so we can use it to calculate indices in matrix. with float becomes error

    modelParams.u_lb            = ones(mpcParams.N,1) * [-1.0  -pi/6]                    # lower bounds on steering
    modelParams.u_ub            = ones(mpcParams.N,1) * [1.0  pi/6]       #1.2           # upper bounds
    modelParams.z_lb            = ones(mpcParams.N+1,1)*[-Inf -Inf -Inf -0.1]                    # lower bounds on states
    #changeMetersforBARC
    modelParams.z_ub            = ones(mpcParams.N+1,1)*[ Inf  Inf  Inf  2.0]       #!! was 2             # upper bounds
    modelParams.l_A             = 0.125
    modelParams.l_B             = 0.125 #0.125

    modelParams.dt              = 0.1

    posInfo.s_start             = 0.0
    posInfo.s_target            = 5.0

    oldTraj.oldTraj             = zeros(buffersize,4,2)
    # oldTraj.oldTraj[:,1,1]      = 0.1:0.1:buffersize/10 # are these necessarxy?
    # oldTraj.oldTraj[:,1,2]      = 1:buffersize
    oldTraj.oldTrajXY           = zeros(buffersize,4,2)
    oldTraj.oldInput            = zeros(buffersize,2,2)
    oldTraj.oldCost             = 100*ones(Int64,2)                   # dummies for initialization

    mpcCoeff.order              = 5
    mpcCoeff.coeffCost          = zeros(mpcCoeff.order+1,2)
    mpcCoeff.coeffConst         = zeros(mpcCoeff.order+1,2,3) # nz-1 because no coeff for s
    mpcCoeff.pLength            = 4*mpcParams.N        # small values here may lead to numerical problems since the functions are only approximated in a short horizon

    lapStatus.currentLap        = 1         # initialize lap number
    lapStatus.currentIt         = 0         # current iteration in lap #?? why not 1?
end