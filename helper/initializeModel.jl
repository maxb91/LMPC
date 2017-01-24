function InitializeModel(m::classes.MpcModel,mpcParams::classes.MpcParams,modelParams::classes.ModelParams,trackCoeff::classes.TrackCoeff,z_Init::Array{Float64,1}, 
    obstacle::classes.Obstacle, oldTraj::classes.OldTrajectory)

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
    Q_velocity      = mpcParams.Q_velocity
    R               = mpcParams.R # cost for control is always used but curently 0
    coeffTermCost   = mpcCoeff.coeffCost[1,:,:]::Array{Float64,2}
    coeffTermConst  = mpcCoeff.coeffConst[1,:,:,:]::Array{Float64,3}
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
    z_Ref           = cat(2,v_ref*ones(N+1,1),zeros(N+1,3),0.0*ones(N+1,1),s_target*ones(N+1,1))     # Reference trajectory: path following -> stay on line and keep constant velocity
    u_Ref           = zeros(N,2)
    
    #########################################################
    m.mdl = Model(solver = IpoptSolver(print_level=0))#, max_cpu_time=10.0))#mu_strategy=adaptive,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))
    #########################################################
    # c_Vx            = [-0.012521482551127934,-0.12341315079450611,0.24925430976232502] 
    # c_Vy            = [-0.4489520316881863,-0.003816534068778571,0.11170845000402227,-0.16451185081929146]
    # c_Psi           = [-0.6796565974307567,-0.2094159184787298, 2.84751043369531]

    # c_Vx            = [-0.07180500739657202,-0.09745914885807667,0.20428541637830308] 
    # c_Vy            = [-0.5604464371949326,-0.0026778881078367424,0.10449300305620522,-0.16279116816645695]
    # c_Psi           = [-0.4828479437301086,-0.25196103153406685, 2.191800531472544]


    m.z_Ol = @variable( m.mdl, [i =1:(N+1),j = 1:6], lowerbound = modelParams.z_lb[i,j], upperbound = modelParams.z_ub[i,j])      # z = s, ey, epsi, v
    m.u_Ol = @variable( m.mdl, [i=1:N,j=1:2], lowerbound = modelParams.u_lb[i,j], upperbound = modelParams.u_ub[i,j])          # overwrtie dim of in classes.jl?
    m.lambda = @variable( m.mdl, [1:n_oldTraj],lowerbound = 0, upperbound = 1)
    m.eps = @variable( m.mdl,[1:3], lowerbound = 0)

    # for i=1:4
    #     for j=1:N+1
    #         setlowerbound(m.z_Ol[j,i], modelParams.z_lb[j,i])
    #         setupperbound(m.z_Ol[j,i], modelParams.z_ub[j,i])
    #     end
    # end
    @NLparameter(m.mdl,m.ssInfOn[1:n_oldTraj]== 1)
    @NLparameter(m.mdl, m.z0[i=1:6] == z_Init[i])
    @NLconstraint(m.mdl, [i=1:6], m.z_Ol[1,i] == m.z0[i])

    #@constraint(m.mdl, m.lambda[1]+m.lambda[2]+m.lambda[3]+m.lambda[4]+m.lambda[5]== 1)
    @NLconstraint(m.mdl, sum(m.lambda[j] for j=1:n_oldTraj)== 1)
    @NLconstraint(m.mdl,[i = 1:(N+1)], m.z_Ol[i,5] <=  trackCoeff.width/2 + m.eps[1] )
    @NLconstraint(m.mdl,[i = 1:(N+1)], m.z_Ol[i,5] >= -trackCoeff.width/2 - m.eps[2] )

    v_max = 2.0
    @NLconstraint(m.mdl,[i = 1:(N+1)], m.z_Ol[i,1] <= v_max +m.eps[3] )
    # object avoidance constraint now done with slack cost instead slack constraint
    # set s and y obst as nlparamter
    #@NLconstraint(m.mdl, [i=1:N+1], ((m.z_Ol[i,1]-m.sCoord_obst[1])/rs )^2 + ( (m.z_Ol[i,2]-m.sCoord_obst[2])/ry )^2 >= 1 - m.eps[3])
   #1/m.t

    @NLparameter(m.mdl, m.coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])
    @NLparameter(m.mdl, m.uCurr[i=1:2] == 0)
    @NLparameter(m.mdl, m.sCoord_obst[j=1:(N+1),i=1:2] == 100)
    @NLparameter(m.mdl, m.coeffTermConst[i=1:order+1,k=1:n_oldTraj,j=1:5] == coeffTermConst[i,k,j])
    @NLparameter(m.mdl, m.coeffTermCost[i=1:order+1,k=1:n_oldTraj] == coeffTermCost[i,k])

    @NLexpression(m.mdl, m.c[i = 1:N],    m.coeff[1]*m.z_Ol[i,6]^4+m.coeff[2]*m.z_Ol[i,6]^3+m.coeff[3]*m.z_Ol[i,6]^2+m.coeff[4]*m.z_Ol[i,6]+m.coeff[5])
    #@NLexpression(m.mdl, m.c[i = 1:N],    sum{m.coeff[j]*m.z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + m.coeff[n_poly_curv+1]) 
    #@NLexpression(m.mdl, m.c[i = 1:N],0)
    # @NLexpression(m.mdl, m.bta[i = 1:N],  atan( L_b / (L_a + L_b) * tan(m.u_Ol[i,2]) ) )
    #@NLexpression(m.mdl, m.bta[i = 1:N],0)
    # @NLexpression(m.mdl, m.dsdt[i = 1:N], m.z_Ol[i,4]*cos(m.z_Ol[i,3]+m.bta[i])/(1-m.z_Ol[i,2]*m.c[i]))
    @NLexpression(m.mdl, m.dsdt[i = 1:N], (m.z_Ol[i,1]*cos(m.z_Ol[i,4]) - m.z_Ol[i,2]*sin(m.z_Ol[i,4]))/(1-m.z_Ol[i,5]*m.c[i]))
    # @NLexpression(m.mdl, m.dsdt[i = 1:N], (m.z_Ol[i,1]*1 - m.z_Ol[i,2]*(m.z_Ol[i,4]))/(1-m.z_Ol[i,5]*m.c[i]))
    # System dynamics
    
    # @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,1]  == m.z_Ol[i,1] + c_Vx[1]*m.z_Ol[i,2]*m.z_Ol[i,3] + c_Vx[2]*m.z_Ol[i,1] + c_Vx[3]*m.u_Ol[i,1]) #v_x
    # @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,2]  == m.z_Ol[i,2] + c_Vy[2]*m.z_Ol[i,1]*m.z_Ol[i,3] + c_Vy[1]*m.z_Ol[i,2]/m.z_Ol[i,1] + c_Vy[3]*m.z_Ol[i,3]/m.z_Ol[i,1] + c_Vy[4]*m.u_Ol[i,2]) #v_y
    # @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,3]  == m.z_Ol[i,3] + c_Psi[1]*m.z_Ol[i,3]/m.z_Ol[i,1] + c_Psi[2]*m.z_Ol[i,2]/m.z_Ol[i,1] +c_Psi[3]*m.u_Ol[i,2] )                     # psi_dot
    
    mass = 1.98 # kg
    mu  = 0.85
    g = 9.81 # m/s^2
    I_z = 0.03 # kg * m^2
    B = 1.0#1.0
    C = 1.25#1.25
    FMax = mu*mass*g / 2.0 
    l_A             = 0.125
    l_B             = 0.125   
    # if F_xr > FMax      
    #     F_xr = FMax    
    # elseif F_xr < -FMax  
    #     F_xr = -FMax 
    # end

    # determine slip angles  
    # m.z_Ol is always<=0.1 see boundaries    
    # if v_x < 0.1        
    #     alpha_f = 0.0        
    #     alpha_r = 0.0      
    # else        
    # @NLexpression(m.mdl, alpha_f[i = 1:N], atan( (m.z_Ol[i,2]+l_A*m.z_Ol[i,3]) / m.z_Ol[i,1] ) - m.u_Ol[i,2])
    # @NLexpression(m.mdl, alpha_r[i = 1:N], atan( (m.z_Ol[i,2]-l_B*m.z_Ol[i,3]) / m.z_Ol[i,1]) )   
    # @NLexpression(m.mdl, alpha_f[i = 1:N],  (m.z_Ol[i,2]+l_A*m.z_Ol[i,3]) / m.z_Ol[i,1]  - m.u_Ol[i,2])
    # @NLexpression(m.mdl, alpha_r[i = 1:N],  (m.z_Ol[i,2]-l_B*m.z_Ol[i,3]) / m.z_Ol[i,1] )   
    # end
    # if exact_sim_i==1 && max(abs(alpha_f),abs(alpha_r))>30/180*pi
    #     warn("Large slip angles: alpha_f = $(alpha_f*180/pi)°, alpha_r = $(alpha_r*180/pi)° , x =$x, y = $y")
    # end
    
    # @NLexpression(m.mdl, F_yf[i = 1:N], -FMax * sin(C*atan(B*atan( (m.z_Ol[i,2]+l_A*m.z_Ol[i,3]) / m.z_Ol[i,1] ) - m.u_Ol[i,2])))
    # @NLexpression(m.mdl, F_yr[i = 1:N], -FMax * sin(C*atan(B*atan( (m.z_Ol[i,2]-l_B*m.z_Ol[i,3]) / m.z_Ol[i,1]))))

    @NLexpression(m.mdl, F_yf[i = 1:N], -FMax * sin(C*atan(B*( (m.z_Ol[i,2]+l_A*m.z_Ol[i,3]) / m.z_Ol[i,1] ) - m.u_Ol[i,2])))
    @NLexpression(m.mdl, F_yr[i = 1:N], -FMax * sin(C*atan(B*( (m.z_Ol[i,2]-l_B*m.z_Ol[i,3]) / m.z_Ol[i,1]))))
    # @NLexpression(m.mdl, F_yf[i = 1:N], -FMax * 1.02*atan( (m.z_Ol[i,2]+l_A*m.z_Ol[i,3]) / m.z_Ol[i,1] ) - m.u_Ol[i,2] )
    # @NLexpression(m.mdl, F_yr[i = 1:N], -FMax * 1.02*atan( (m.z_Ol[i,2]-l_B*m.z_Ol[i,3]) / m.z_Ol[i,1]) )

    max_alpha = 15

    @NLconstraint(m.mdl, [i=1:N], atan( (m.z_Ol[i,2]-l_B*m.z_Ol[i,3]) / m.z_Ol[i,1]) <=  max_alpha/180*pi)
    @NLconstraint(m.mdl, [i=1:N], atan( (m.z_Ol[i,2]-l_B*m.z_Ol[i,3]) / m.z_Ol[i,1]) >= -max_alpha/180*pi)
    @NLconstraint(m.mdl, [i=1:N], atan( (m.z_Ol[i,2]+l_A*m.z_Ol[i,3]) / m.z_Ol[i,1] ) - m.u_Ol[i,2] <=  max_alpha/180*pi)
    @NLconstraint(m.mdl, [i=1:N], atan( (m.z_Ol[i,2]+l_A*m.z_Ol[i,3]) / m.z_Ol[i,1] ) - m.u_Ol[i,2] >= -max_alpha/180*pi)

    # if F_yr > sqrt(FxMax^2 - F_xr^2)        
    #     F_yr = sqrt(FxMax^2 - F_xr^2)    
    # elseif  F_yr < -sqrt(FxMax^2 - F_xr^2)  
    #     F_yr = -sqrt(FxMax^2 - F_xr^2) 
    # end


    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,1]  == m.z_Ol[i,1] + dt*(m.z_Ol[i,3]*m.z_Ol[i,2]+m.u_Ol[i,1]-1/mass*(F_yf[i]*sin(m.u_Ol[i,2])))) #v_x
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,2]  == m.z_Ol[i,2] + dt*(-m.z_Ol[i,3]*m.z_Ol[i,1]+1/mass*(F_yf[i]*cos(m.u_Ol[i,2])+F_yr[i]))) #v_y
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,3]  == m.z_Ol[i,3] + dt*(1/I_z*(l_A*F_yf[i]*cos(m.u_Ol[i,2])-l_B*F_yr[i])))    # psi_dot             
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,4]  == m.z_Ol[i,4] + dt*(m.z_Ol[i,3]-m.dsdt[i]*m.c[i]) )            # epsi
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,5]  == m.z_Ol[i,5] + dt*(m.z_Ol[i,1]*sin(m.z_Ol[i,4])+m.z_Ol[i,2]*cos(m.z_Ol[i,4])))#e_y
    # @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,5]  == m.z_Ol[i,5] + dt*(m.z_Ol[i,1]*m.z_Ol[i,4]+m.z_Ol[i,2]))#e_y
    @NLconstraint(m.mdl, [i=1:N], m.z_Ol[i+1,6]  == m.z_Ol[i,6] + dt*m.dsdt[i]) #s



    #define expressions for cost
    # Derivative cost
    # ---------------------------------
    @NLexpression(m.mdl, derivCost, sum(QderivZ[j] * sum((m.z_Ol[i,j] - m.z_Ol[i + 1,j]) ^ 2 for i = 1:N) for j = 1:6) +  #(m.z0[j]-m.z_Ol[1,j])^2
                                      sum(QderivU[j]*((m.uCurr[j]-m.u_Ol[1,j])^2+sum((m.u_Ol[i,j]-m.u_Ol[i+1,j])^2 for i=1:N-1)) for j=1:2))
    m.derivCost= derivCost
    # Lane cost
    # ---------------------------------
    #@NLexpression(m.mdl, laneCost, Q_lane*sum{m.z_Ol[i,2]^2*((0.5+0.5*tanh(35*(m.z_Ol[i,2]-ey_max-0.09))) + (0.5-0.5*tanh(35*(m.z_Ol[i,2]+ey_max+0.09)))),i=1:N+1})
    @NLexpression(m.mdl, laneCost, Q_lane*sum(10.0*m.eps[i]+7.0*m.eps[i]^2 for i=1:2))
    m.laneCost = laneCost


    # soft constraint for max velocity
    @NLexpression(m.mdl, velocityCost, Q_velocity*(10.0*m.eps[3]+7.0*m.eps[3]^2)  )
    m.velocityCost = velocityCost
    # Control Input cost
    # ---------------------------------
    @NLexpression(m.mdl, controlCost, 0.5*sum(R[j]*sum((m.u_Ol[i,j]-u_Ref[i,j])^2 for i=1:N) for j=1:2))
    m.controlCost = controlCost

    # Terminal constraints (soft), starting from 2nd lap
    #constraints force trajectory to end up on ss
    # ---------------------------------  
   #@NLexpression(m.mdl, constZTerm, sum{Q_term[j]* sum{m.ssInfOn[k]*m.lambda[k]*(sum{m.coeffTermConst[i,k,j]*m.z_Ol[N+1,1]^(order+1-i),i=1:order+1}-m.z_Ol[N+1,j+1]),k=1:n_oldTraj}^2,j=1:3})                           
    
    @NLexpression(m.mdl, constZTerm, sum(Q_term[j]* (sum(m.ssInfOn[k]*m.lambda[k]*sum(m.coeffTermConst[i,k,j]*m.z_Ol[N+1,6]^(order+1-i) for i=1:order+1) for k=1:n_oldTraj)-m.z_Ol[N+1,j])^2 for j=1:5))                           
    m.constZTerm = constZTerm           
   #basic idea:
    #@NLexpression(m.mdl, constZTerm, sum{Q_term[j]*(sum{coeffTermConst[i,1,j]*m.z_Ol[N+1,1]^(order+1-i),i=1:order+1}-m.z_Ol[N+1,j+1])^2,j=1:3})


    # Terminal cost
    # ---------------------------------
    # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
    @NLexpression(m.mdl, costZTerm,  sum(m.ssInfOn[k]*m.lambda[k]*sum(m.coeffTermCost[i,k]*m.z_Ol[N+1,6]^(order+1-i) for i=1:order+1) for k=1:n_oldTraj))
    m.costZTerm = costZTerm
    #basic idea    
    #@NLexpression(m.mdl, costZTerm, Q_cost*sum{coeffTermCost[i,1]*m.z_Ol[N+1,1]^(order+1-i),i=1:order+1})

    # State cost
    # ---------------------------------
    # if we're in the first lap, just do path following
    @NLexpression(m.mdl, costPath, 0.5*sum(Q[i]*sum((m.z_Ol[j,i]-z_Ref[j,i])^2for j=1:N+1) for i=1:6))    # Follow trajectory
    m.costPath = costPath
    #put cost on z (actually should put cost only on z before finishing the lap)
    #@NLexpression(m.mdl, costZ_h, 0)          # zero state cost after crossing the finish line
    #@NLexpression(m.mdl, costZ, 1 + (costZ_h-1) * (0.5+0.5*tanh(50*(m.z_Ol[1,N+1]+s_start-s_target))))
    #@NLexpression(m.mdl, costZ, Q_cost*1)
    @NLexpression(m.mdl, costZ, Q_cost*sum(1 for i=1:N+1))

    m.costZ = costZ
    ## Cost to avoid obstacle. increases when car is near obstacle currently implemented as : a *1/(0.1+cost)
    @NLexpression(m.mdl, costObstacle, 
        # sum{-Q_obstacleNumer*log(0.1+ (( (m.z_Ol[i,1]-m.sCoord_obst[i,1])/rs )^2 + ( (m.z_Ol[i,2]-m.sCoord_obst[i,2])/ry) ^2 - 1)^4)+
        #   -2*Q_obstacleNumer*log(( (m.z_Ol[i,1]-m.sCoord_obst[i,1])/rs )^2 + ( (m.z_Ol[i,2]-m.sCoord_obst[i,2])/ry) ^2),i=1:N+1})

        # sum{Q_obstacleNumer*1/(0.01+(Q_obstacle* (( (m.z_Ol[i,1]-m.sCoord_obst[i,1])/rs )^2 + ( (m.z_Ol[i,2]-m.sCoord_obst[i,2])/ry) ^2 - 1))^4),i=1:N+1})

        sum(Q_obstacleNumer*1/(0.01+(Q_obstacle* (( (m.z_Ol[i,6]-m.sCoord_obst[i,1])/rs )^2 + ( (m.z_Ol[i,5]-m.sCoord_obst[i,2])/ry) ^2 - 1))^4)+
        Q_obstacleNumer*3/(0.01+0.6*(((m.z_Ol[i,6]-m.sCoord_obst[i,1])/rs)^2+((m.z_Ol[i,5]-m.sCoord_obst[i,2])/ry)^2)) for i=1:N+1))



    #3*Q_obstacleNumer/((((m.z_Ol[i,1]-m.sCoord_obst[i,1])/rs)^2+((m.z_Ol[i,2]-m.sCoord_obst[i,2])/ry)^2)),i=1:N+1})
    #@NLexpression(m.mdl, costObstacle,    (10*m.eps[3]+5.0*m.eps[3]^2))
    #@NLexpression(m.mdl, costObstacle,    0)       
    m.costObstacle = costObstacle
    #@NLobjective(m.mdl, Min, m.costPath)
end
