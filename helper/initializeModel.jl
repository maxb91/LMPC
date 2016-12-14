function InitializeModel(m::classes.MpcModel,mpcParams::classes.MpcParams,modelParams::classes.ModelParams,trackCoeff::classes.TrackCoeff,z_Init::Array{Float64,1}, obstacle::classes.Obstacle, oldTraj::classes.OldTrajectory)

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
    rs = obstacle.rs
    ry = obstacle.ry
    Q               = mpcParams.Q #Cost of states just for path following
    Q_term          = mpcParams.Q_term
    Q_cost          = mpcParams.Q_cost
    Q_obstacle      = mpcParams.Q_obstacle
    Q_obstacleNumer = mpcParams.Q_obstacleNumer
    Q_lane          = mpcParams.Q_lane
    R               = mpcParams.R # cost for control is always used but curently 0
    coeffTermCost   = mpcCoeff.coeffCost::Array{Float64,2}
    coeffTermConst  = mpcCoeff.coeffConst::Array{Float64,3}
    order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
    s_start         = posInfo.s_start
    s_target        = posInfo.s_target
    QderivZ         = mpcParams.QderivZ::Array{Float64,1}
    QderivU         = mpcParams.QderivU::Array{Float64,1}
    v_ref           = mpcParams.vPathFollowing

    N           = mpcParams.N
    ey_max      = trackCoeff.width/2
    n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation
    n_oldTraj   = oldTraj.n_oldTraj
    # Create function-specific parameters
    local z_Ref::Array{Float64,2}
    z_Ref           = cat(2,s_target*ones(N+1,1),zeros(N+1,2),v_ref*ones(N+1,1))     # Reference trajectory: path following -> stay on line and keep constant velocity
    u_Ref           = zeros(N,2)
    
    #########################################################
    m.mdl = Model(solver = IpoptSolver(print_level=0))#mu_strategy=adaptive,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))
    #########################################################


    @variable( m.mdl, modelParams.z_lb[i,j]<= m.z_Ol[i =1:(N+1),j = 1:4]<= modelParams.z_ub[i,j])      # z = s, ey, epsi, v
    @variable( m.mdl,  modelParams.u_lb[i,j] <= m.u_Ol[i=1:N,j=1:2]<= modelParams.u_ub[i,j])          # overwrtie dim of in classes.jl?
    @variable( m.mdl, 0 <= m.lambda[1:n_oldTraj] <= 1)
    @variable( m.mdl, m.eps[1:2] >=0)

    # for i=1:4
    #     for j=1:N+1
    #         setlowerbound(m.z_Ol[j,i], modelParams.z_lb[j,i])
    #         setupperbound(m.z_Ol[j,i], modelParams.z_ub[j,i])
    #     end
    # end
    @NLparameter(m.mdl,m.ssInfOn[1:n_oldTraj]== 1)
    #@NLparameter(m.mdl,m.ssAddTraj[i=1:1] == 1500)
    @NLparameter(m.mdl, m.z0[i=1:4] == z_Init[i])
    @NLconstraint(m.mdl, [i=1:4], m.z_Ol[1,i] == m.z0[i])

    #@constraint(m.mdl, m.lambda[1]+m.lambda[2]+m.lambda[3]+m.lambda[4]+m.lambda[5]== 1)
    @NLconstraint(m.mdl, sum{m.lambda[j],j=1:n_oldTraj}== 1)
    @NLconstraint(m.mdl,[i = 1:(N+1)], m.z_Ol[i,2] <=  trackCoeff.width/2 + m.eps[1] )
    @NLconstraint(m.mdl,[i = 1:(N+1)], m.z_Ol[i,2] >= -trackCoeff.width/2 - m.eps[2] )
    # object avoidance constraint now done with slack cost instead slack constraint
    # set s and y obst as nlparamter
    #@NLconstraint(m.mdl, [i=1:N+1], ((m.z_Ol[i,1]-m.sCoord_obst[1])/rs )^2 + ( (m.z_Ol[i,2]-m.sCoord_obst[2])/ry )^2 >= 1 - m.eps[3])
   #1/m.t

    @NLparameter(m.mdl, m.coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])
    @NLparameter(m.mdl, m.uCurr[i=1:2] == 0)
    @NLparameter(m.mdl, m.sCoord_obst[j=1:(N+1),i=1:2] == 100)
    @NLparameter(m.mdl, m.coeffTermConst[i=1:order+1,k=1:n_oldTraj,j=1:3] == coeffTermConst[i,k,j])
    @NLparameter(m.mdl, m.coeffTermCost[i=1:order+1,k=1:n_oldTraj] == coeffTermCost[i,k])

    @NLexpression(m.mdl, m.c[i = 1:N],  m.coeff[1]*m.z_Ol[i,1]^4+m.coeff[2]*m.z_Ol[i,1]^3+m.coeff[3]*m.z_Ol[i,1]^2+m.coeff[4]*m.z_Ol[i,1]+m.coeff[5])
    #@NLexpression(m.mdl, m.c[i = 1:N],    sum{m.coeff[j]*m.z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + m.coeff[n_poly_curv+1]) 
    #@NLexpression(m.mdl, m.c[i = 1:N],0)
    @NLexpression(m.mdl, m.bta[i = 1:N],  atan( L_b / (L_a + L_b) * tan(m.u_Ol[i,2]) ) )
    #@NLexpression(m.mdl, m.bta[i = 1:N],0)
    @NLexpression(m.mdl, m.dsdt[i = 1:N], m.z_Ol[i,4]*cos(m.z_Ol[i,3]+m.bta[i])/(1-m.z_Ol[i,2]*m.c[i]))
    # System dynamics
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,1]  == m.z_Ol[i,1] + dt*m.dsdt[i]  )
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,2]  == m.z_Ol[i,2] + dt*m.z_Ol[i,4]*sin(m.z_Ol[i,3]+m.bta[i])  )                     # ey
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,3]  == m.z_Ol[i,3] + dt*(m.z_Ol[i,4]/L_b*sin(m.bta[i])-m.dsdt[i]*m.c[i])  )            # epsi
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,4]  == m.z_Ol[i,4] + dt*(m.u_Ol[i,1]))#- 0.23*abs(m.z_Ol[i,4]) * m.z_Ol[i,4]))#0.63  # v
    # @NLconstraint(m.mdl, [i=1:N+1], m.z_Ol[i,5]  == m.sCoord_obst[i,1] - m.z_Ol[i,1]  )
    # @NLconstraint(m.mdl, [i=1:N+1], m.z_Ol[i,6]  == sy_obst - m.z_Ol[i,2] )
 


    #define expressions for cost
    # Derivative cost
    # ---------------------------------
    @NLexpression(m.mdl, derivCost, sum{QderivZ[j]*(sum{(m.z_Ol[i,j]-m.z_Ol[i+1,j])^2,i=1:N}),j=1:4} +  #(m.z0[j]-m.z_Ol[1,j])^2
                                      sum{QderivU[j]*((m.uCurr[j]-m.u_Ol[1,j])^2+sum{(m.u_Ol[i,j]-m.u_Ol[i+1,j])^2,i=1:N-1}),j=1:2})
    m.derivCost= derivCost
    # Lane cost
    # ---------------------------------
    #@NLexpression(m.mdl, laneCost, Q_lane*sum{m.z_Ol[i,2]^2*((0.5+0.5*tanh(35*(m.z_Ol[i,2]-ey_max-0.09))) + (0.5-0.5*tanh(35*(m.z_Ol[i,2]+ey_max+0.09)))),i=1:N+1})
    @NLexpression(m.mdl, laneCost, Q_lane*sum{10.0*m.eps[i]+7.0*m.eps[i]^2,i=1:2})
    m.laneCost = laneCost
    # Control Input cost
    # ---------------------------------
    @NLexpression(m.mdl, controlCost, 0.5*sum{R[j]*sum{(m.u_Ol[i,j]-u_Ref[i,j])^2,i=1:N},j=1:2})
    m.controlCost = controlCost

    # Terminal constraints (soft), starting from 2nd lap
    #constraints force trajectory to end up on ss
    # ---------------------------------  
   #@NLexpression(m.mdl, constZTerm, sum{Q_term[j]* sum{m.ssInfOn[k]*m.lambda[k]*(sum{m.coeffTermConst[i,k,j]*m.z_Ol[N+1,1]^(order+1-i),i=1:order+1}-m.z_Ol[N+1,j+1]),k=1:n_oldTraj}^2,j=1:3})                           
    
    @NLexpression(m.mdl, constZTerm, sum{Q_term[j]*(sum{m.ssInfOn[k]*m.lambda[k]*sum{m.coeffTermConst[i,k,j]*m.z_Ol[N+1,1]^(order+1-i),i=1:order+1},k=1:n_oldTraj}-m.z_Ol[N+1,j+1])^2,j=1:3})                           
    m.constZTerm = constZTerm           
   #basic idea:
    #@NLexpression(m.mdl, constZTerm, sum{Q_term[j]*(sum{coeffTermConst[i,1,j]*m.z_Ol[N+1,1]^(order+1-i),i=1:order+1}-m.z_Ol[N+1,j+1])^2,j=1:3})


    # Terminal cost
    # ---------------------------------
    # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
    @NLexpression(m.mdl, costZTerm,  sum{m.ssInfOn[k]*m.lambda[k]*sum{m.coeffTermCost[i,k]*m.z_Ol[N+1,1]^(order+1-i),i=1:order+1},k=1:n_oldTraj})
    m.costZTerm = costZTerm
    #basic idea    
    #@NLexpression(m.mdl, costZTerm, Q_cost*sum{coeffTermCost[i,1]*m.z_Ol[N+1,1]^(order+1-i),i=1:order+1})

    # State cost
    # ---------------------------------
    # if we're in the first lap, just do path following
    @NLexpression(m.mdl, costPath, 0.5*sum{Q[i]*sum{(m.z_Ol[j,i]-z_Ref[j,i])^2,j=1:N+1},i=1:4})    # Follow trajectory
    m.costPath = costPath
    #put cost on z (actually should put cost only on z before finishing the lap)
    #@NLexpression(m.mdl, costZ_h, 0)          # zero state cost after crossing the finish line
    #@NLexpression(m.mdl, costZ, 1 + (costZ_h-1) * (0.5+0.5*tanh(50*(m.z_Ol[1,N+1]+s_start-s_target))))
    #@NLexpression(m.mdl, costZ, Q_cost*1)
    @NLexpression(m.mdl, costZ, Q_cost*sum{1,i=1:N+1})

    m.costZ = costZ
    ## Cost to avoid obstacle. increases when car is near obstacle currently implemented as : a *1/(0.1+cost)
    @NLexpression(m.mdl, costObstacle, 
        #sum{Q_obstacleNumer*1/(0.01+(Q_obstacle* (( (m.z_Ol[i,1]-m.sCoord_obst[i,1])/rs )^2 + ( (m.z_Ol[i,2]-m.sCoord_obst[i,2])/ry) ^2 - 1))^4),i=1:N+1})
        sum{Q_obstacleNumer*1/(0.01+(Q_obstacle* (( (m.z_Ol[i,1]-m.sCoord_obst[i,1])/rs )^2 + ( (m.z_Ol[i,2]-m.sCoord_obst[i,2])/ry) ^2 - 1))^4)+
        Q_obstacleNumer*3/(0.01+0.6*(((m.z_Ol[i,1]-m.sCoord_obst[i,1])/rs)^2+((m.z_Ol[i,2]-m.sCoord_obst[i,2])/ry)^2)),i=1:N+1})
        #3*Q_obstacleNumer/((((m.z_Ol[i,1]-m.sCoord_obst[i,1])/rs)^2+((m.z_Ol[i,2]-m.sCoord_obst[i,2])/ry)^2)),i=1:N+1})
    #@NLexpression(m.mdl, costObstacle,    (10*m.eps[3]+5.0*m.eps[3]^2))
    #@NLexpression(m.mdl, costObstacle,    0)       
    m.costObstacle = costObstacle
    #@NLobjective(m.mdl, Min, m.costPath)
end