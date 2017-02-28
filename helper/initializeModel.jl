type initPathFollowingModel
    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}
    uCurr::Array{JuMP.NonlinearParameter,1}
    sCoord_obst::Array{JuMP.NonlinearParameter,2}

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    eps::JuMP.Variable

    dsdt::Array{JuMP.NonlinearExpression,1}
    c::Array{JuMP.NonlinearExpression,1}

    costPath::JuMP.NonlinearExpression
    derivCost::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression

    #soft constraints
    velocityCost::JuMP.NonlinearExpression
    costObstacle::JuMP.NonlinearExpression

    function initPathFollowingModel(mpcParams::classes.MpcParams,modelParams::classes.ModelParams,trackCoeff::classes.TrackCoeff,mpcCoeff::classes.MpcCoeff, 
    obstacle::classes.Obstacle, n_oldTraj::Int64)

        m = new()

        dt        = modelParams.dt
        l_A       = modelParams.l_A
        l_B       = modelParams.l_B
        u_lb      = modelParams.u_lb
        u_ub      = modelParams.u_ub
        z_lb      = modelParams.z_lb
        z_ub      = modelParams.z_ub
        s_obst    = obstacle.s_obstacle
        sy_obst   = obstacle.sy_obstacle
        rs        = obstacle.rs
        ry        = obstacle.ry

        Q               = mpcParams.Q #Cost of states just for path following
        Q_obstacle      = mpcParams.Q_obstacle
        Q_obstacleNumer = mpcParams.Q_obstacleNumer
        QderivZ         = mpcParams.QderivZ::Array{Float64,1}
        QderivU         = mpcParams.QderivU::Array{Float64,1}
        Q_velocity      = mpcParams.Q_velocity
        R               = mpcParams.R # cost for control is always used but curently 0
        order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
        
        v_ref           = mpcParams.vPathFollowing

        N           = mpcParams.N
        ey_max      = trackCoeff.width/2
        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation 
        # Create function-specific parameters
        local z_Ref::Array{Float64,2}
        s_target = 100 # the weight for s in pathfollowing is set to zero, so this value is not used
        z_Ref           = cat(2,v_ref*ones(N+1,1),zeros(N+1,3),0.0*ones(N+1,1),s_target*ones(N+1,1))     # Reference trajectory: path following -> stay on line and keep constant velocity
        u_Ref           = zeros(N,2)
        z_Init          = zeros(6)
        z_Init[1]       = 0.6

        v_max = modelParams.v_max
        max_alpha =modelParams.max_alpha
        mass = modelParams.mass
        mu  = modelParams.mu
        g = modelParams.g
        I_z = modelParams.I_z
        B = modelParams.B
        C = modelParams.C
        FMax = mu*mass*g / 2.0 
 
        #########################################################
        mdl = Model(solver = IpoptSolver(print_level=0))#, max_cpu_time=10.0))#mu_strategy=adaptive,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))
        #########################################################
        # c_Vx            = [-0.012521482551127934,-0.12341315079450611,0.24925430976232502] 
        # c_Vy            = [-0.4489520316881863,-0.003816534068778571,0.11170845000402227,-0.16451185081929146]
        # c_Psi           = [-0.6796565974307567,-0.2094159184787298, 2.84751043369531]

        # c_Vx            = [-0.07180500739657202,-0.09745914885807667,0.20428541637830308] 
        # c_Vy            = [-0.5604464371949326,-0.0026778881078367424,0.10449300305620522,-0.16279116816645695]
        # c_Psi           = [-0.4828479437301086,-0.25196103153406685, 2.191800531472544]


        z_Ol   = @variable( mdl, [i =1:(N+1),j = 1:6], lowerbound = z_lb[i,j], upperbound = z_ub[i,j])      # z = s, ey, epsi, v
        u_Ol   = @variable( mdl, [i=1:N,j=1:2], lowerbound = u_lb[i,j], upperbound = u_ub[i,j])          # overwrtie dim of in classes.jl?
        eps    = @variable( mdl, lowerbound = 0)

     
        @NLparameter(mdl, z0[i=1:6] == z_Init[i])
        @NLconstraint(mdl, [i=1:6], z_Ol[1,i] == z0[i])
        @NLconstraint(mdl,[i = 1:(N+1)], z_Ol[i,1] <= v_max + eps )


        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])
        @NLparameter(mdl, uCurr[i=1:2] == 0)
        @NLparameter(mdl, sCoord_obst[j=1:(N+1),i=1:2] == 100)

        @NLexpression(mdl, c[i = 1:N],    coeff[1]*z_Ol[i,6]^4+coeff[2]*z_Ol[i,6]^3+coeff[3]*z_Ol[i,6]^2+coeff[4]*z_Ol[i,6]+coeff[5])
        #@NLexpression(mdl, c[i = 1:N],    sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1])
        
        @NLexpression(mdl, dsdt[i = 1:N], (z_Ol[i,1]*cos(z_Ol[i,4]) - z_Ol[i,2]*sin(z_Ol[i,4]))/(1-z_Ol[i,5]*c[i]))
        # @NLexpression(mdl, dsdt[i = 1:N], (z_Ol[i,1]*1 - z_Ol[i,2]*(z_Ol[i,4]))/(1-z_Ol[i,5]*c[i]))
       
     

        # System dynamics
        # ---------------------------------
        
        # @NLexpression(mdl, F_yf[i = 1:N], -FMax * sin(C*atan(B*atan( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2])))
        # @NLexpression(mdl, F_yr[i = 1:N], -FMax * sin(C*atan(B*atan( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]))))
        # @NLexpression(mdl, F_yf[i = 1:N], -FMax * sin(C*atan(B*( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2])))
        # @NLexpression(mdl, F_yr[i = 1:N], -FMax * sin(C*atan(B*( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]))))
        @NLexpression(mdl, F_yf[i = 1:N], -2*FMax * (((z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2]))
        @NLexpression(mdl, F_yr[i = 1:N], -2*FMax *  ((z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) )

        

        # @NLconstraint(mdl, [i=1:N], atan( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) <=  max_alpha/180*pi)
        # @NLconstraint(mdl, [i=1:N], atan( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) >= -max_alpha/180*pi)
        # @NLconstraint(mdl, [i=1:N], atan( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2] <=  max_alpha/180*pi)
        # @NLconstraint(mdl, [i=1:N], atan( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2] >= -max_alpha/180*pi)
        @NLconstraint(mdl, [i=1:N], ( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) <=  max_alpha/180*pi)
        @NLconstraint(mdl, [i=1:N], ( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) >= -max_alpha/180*pi)
        @NLconstraint(mdl, [i=1:N], ( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2] <=  max_alpha/180*pi)
        @NLconstraint(mdl, [i=1:N], ( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2] >= -max_alpha/180*pi)



        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + c_Vx[1]*z_Ol[i,2]*z_Ol[i,3] + c_Vx[2]*z_Ol[i,1] + c_Vx[3]*u_Ol[i,1]) #v_x
        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + c_Vy[2]*z_Ol[i,1]*z_Ol[i,3] + c_Vy[1]*z_Ol[i,2]/z_Ol[i,1] + c_Vy[3]*z_Ol[i,3]/z_Ol[i,1] + c_Vy[4]*u_Ol[i,2]) #v_y
        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + c_Psi[1]*z_Ol[i,3]/z_Ol[i,1] + c_Psi[2]*z_Ol[i,2]/z_Ol[i,1] +c_Psi[3]*u_Ol[i,2] )# psi_dot
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + dt*(z_Ol[i,3]*z_Ol[i,2]+u_Ol[i,1]-1/mass*(F_yf[i]*sin(u_Ol[i,2])))) #v_x
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + dt*(-z_Ol[i,3]*z_Ol[i,1]+1/mass*(F_yf[i]*cos(u_Ol[i,2])+F_yr[i]))) #v_y
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + dt*(1/I_z*(l_A*F_yf[i]*cos(u_Ol[i,2])-l_B*F_yr[i])))    # psi_dot              
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,4]  == z_Ol[i,4] + dt*(z_Ol[i,3]-dsdt[i]*c[i]) )            # epsi
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,5]  == z_Ol[i,5] + dt*(z_Ol[i,1]*sin(z_Ol[i,4])+z_Ol[i,2]*cos(z_Ol[i,4])))#e_y
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,6]  == z_Ol[i,6] + dt*dsdt[i]) #s



        # Define expressions for cost
        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum(QderivZ[j] * sum((z_Ol[i,j] - z_Ol[i + 1,j]) ^ 2 for i = 1:N) for j = 1:6) +
                                          sum(QderivU[j]*((uCurr[j]-u_Ol[1,j])^2+sum((u_Ol[i,j]-u_Ol[i+1,j])^2 for i=1:N-1)) for j=1:2))
        


        # soft constraint for max velocity
        @NLexpression(mdl, velocityCost, Q_velocity*(10.0*eps+7.0*eps^2)  )

        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*sum(R[j]*sum((u_Ol[i,j]-u_Ref[i,j])^2 for i=1:N) for j=1:2))
   
        # Pathfollowing
        # ---------------------------------
        @NLexpression(mdl, costPath, 0.5*sum(Q[i]*sum((z_Ol[j,i]-z_Ref[j,i])^2for j=1:N+1) for i=1:6))    # Follow trajectory

        # Cost to avoid obstacle
        # ---------------------------------
        @NLexpression(mdl, costObstacle, sum(Q_obstacleNumer*1/(0.01+(Q_obstacle* (( (z_Ol[i,6]-sCoord_obst[i,1])/rs )^2 + ( (z_Ol[i,5]-sCoord_obst[i,2])/ry) ^2 - 1))^4)+
        Q_obstacleNumer*3/(0.01+0.6*(((z_Ol[i,6]-sCoord_obst[i,1])/rs)^2+((z_Ol[i,5]-sCoord_obst[i,2])/ry)^2)) for i=1:N+1))


        #Objective formulation
        @NLobjective(mdl, Min, costPath + derivCost + controlCost + velocityCost + costObstacle)

        m.mdl = mdl
        m.coeff = coeff # curvature coefficients
        m.uCurr = uCurr #last applied input
        m.sCoord_obst = sCoord_obst
        m.u_Ol = u_Ol
        m.z_Ol = z_Ol
        m.eps = eps
        m.z0  = z0
        m.derivCost= derivCost
        m.velocityCost = velocityCost
        m.controlCost = controlCost 
        m.costObstacle = costObstacle
        m.costPath = costPath

        return m
    end
end



type initLearningModel
    mdl::JuMP.Model

    ssInfOn::Array{JuMP.NonlinearParameter,1}
    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}
    uCurr::Array{JuMP.NonlinearParameter,1}
    sCoord_obst::Array{JuMP.NonlinearParameter,2}
    coeffTermConst::Array{JuMP.NonlinearParameter,3}
    coeffTermCost::Array{JuMP.NonlinearParameter,2}
    #s_startC::JuMP.NonlinearParameter

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    lambda::Array{JuMP.Variable,1}
    eps::Array{JuMP.Variable,1}

    dsdt::Array{JuMP.NonlinearExpression,1}
    c::Array{JuMP.NonlinearExpression,1}

    costZ::JuMP.NonlinearExpression
    costZTerm::JuMP.NonlinearExpression
    derivCost::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression

    #soft constraints
    constZTerm::JuMP.NonlinearExpression
    laneCost::JuMP.NonlinearExpression
    velocityCost::JuMP.NonlinearExpression
    costObstacle::JuMP.NonlinearExpression

    function initLearningModel(mpcParams::classes.MpcParams,modelParams::classes.ModelParams,trackCoeff::classes.TrackCoeff,mpcCoeff::classes.MpcCoeff, 
    obstacle::classes.Obstacle, n_oldTraj::Int64)

        m = new()

        dt        = modelParams.dt
        l_A       = modelParams.l_A
        l_B       = modelParams.l_B
        u_lb      = modelParams.u_lb
        u_ub      = modelParams.u_ub
        z_lb      = modelParams.z_lb
        z_ub      = modelParams.z_ub
        s_obst    = obstacle.s_obstacle
        sy_obst    = obstacle.sy_obstacle
        rs = obstacle.rs
        ry = obstacle.ry
        # Q               = mpcParams.Q #Cost of states just for path following
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
        QderivZ         = mpcParams.QderivZ::Array{Float64,1}
        QderivU         = mpcParams.QderivU::Array{Float64,1}
        v_ref           = mpcParams.vPathFollowing

        N           = mpcParams.N
        ey_max      = trackCoeff.width/2
        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation 
        # Create function-specific parameters

        u_Ref           = zeros(N,2)

        z_Init       = zeros(N+1,6)
        z_Init[:,1]  = v_ref*ones(N+1)

        v_max = modelParams.v_max
        max_alpha =modelParams.max_alpha
        mass = modelParams.mass
        mu  = modelParams.mu
        g = modelParams.g
        I_z = modelParams.I_z
        B = modelParams.B
        C = modelParams.C
        FMax = mu*mass*g / 2.0 

        #########################################################
        mdl = Model(solver = IpoptSolver(print_level=0))#, max_cpu_time=10.0))#mu_strategy=adaptive,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))
        #########################################################
        # c_Vx            = [-0.012521482551127934,-0.12341315079450611,0.24925430976232502] 
        # c_Vy            = [-0.4489520316881863,-0.003816534068778571,0.11170845000402227,-0.16451185081929146]
        # c_Psi           = [-0.6796565974307567,-0.2094159184787298, 2.84751043369531]

        # c_Vx            = [-0.07180500739657202,-0.09745914885807667,0.20428541637830308] 
        # c_Vy            = [-0.5604464371949326,-0.0026778881078367424,0.10449300305620522,-0.16279116816645695]
        # c_Psi           = [-0.4828479437301086,-0.25196103153406685, 2.191800531472544]


        z_Ol   = @variable( mdl, [i =1:(N+1),j = 1:6], lowerbound = z_lb[i,j], upperbound = z_ub[i,j], start =z_Init[i,j])      # z = s, ey, epsi, v
        u_Ol   = @variable( mdl, [i=1:N,j=1:2], lowerbound = u_lb[i,j], upperbound = u_ub[i,j], start = 0)          # overwrtie dim of in classes.jl?
        lambda = @variable( mdl, [1:n_oldTraj],lowerbound = 0, upperbound = 1, start = 1/n_oldTraj)
        eps    = @variable( mdl,[1:3], lowerbound = 0, start = 0)

     
        @NLparameter(mdl,ssInfOn[1:n_oldTraj]== 1)
        @NLparameter(mdl, z0[i=1:6] == z_Init[1,i])
        @NLconstraint(mdl, [i=1:6], z_Ol[1,i] == z0[i])

        #@constraint(mdl, lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]== 1)
        @NLconstraint(mdl, sum(lambda[j] for j=1:n_oldTraj)== 1)
        @NLconstraint(mdl,[i = 1:(N+1)], z_Ol[i,5] <=  trackCoeff.width/2 + eps[1] )
        @NLconstraint(mdl,[i = 1:(N+1)], z_Ol[i,5] >= -trackCoeff.width/2 - eps[2] )

        
        @NLconstraint(mdl,[i = 1:(N+1)], z_Ol[i,1] <= v_max + eps[3] )

        # object avoidance constraint now done with slack cost instead slack constraint
        #@NLconstraint(mdl, [i=1:N+1], ((z_Ol[i,1]-sCoord_obst[1])/rs )^2 + ( (z_Ol[i,2]-sCoord_obst[2])/ry )^2 >= 1 - eps[3])


        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])
        @NLparameter(mdl, uCurr[i=1:2] == 0)
        @NLparameter(mdl, sCoord_obst[j=1:(N+1),i=1:2] == 100)
        @NLparameter(mdl, coeffTermConst[i=1:order+1,k=1:n_oldTraj,j=1:5] == coeffTermConst[i,k,j])
        @NLparameter(mdl, coeffTermCost[i=1:order+1,k=1:n_oldTraj] == coeffTermCost[i,k])

        @NLexpression(mdl, c[i = 1:N],    coeff[1]*z_Ol[i,6]^4+coeff[2]*z_Ol[i,6]^3+coeff[3]*z_Ol[i,6]^2+coeff[4]*z_Ol[i,6]+coeff[5])
        #@NLexpression(mdl, c[i = 1:N],    sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1]) 

        @NLexpression(mdl, dsdt[i = 1:N], (z_Ol[i,1]*cos(z_Ol[i,4]) - z_Ol[i,2]*sin(z_Ol[i,4]))/(1-z_Ol[i,5]*c[i]))
        # @NLexpression(mdl, dsdt[i = 1:N], (z_Ol[i,1]*1 - z_Ol[i,2]*(z_Ol[i,4]))/(1-z_Ol[i,5]*c[i]))
        
        
        # System dynamics
        # ---------------------------------
        # @NLexpression(mdl, F_yf[i = 1:N], -FMax * sin(C*atan(B*atan( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2])))
        # @NLexpression(mdl, F_yr[i = 1:N], -FMax * sin(C*atan(B*atan( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]))))
        # @NLexpression(mdl, F_yf[i = 1:N], -FMax * sin(C*atan(B*( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2])))
        # @NLexpression(mdl, F_yr[i = 1:N], -FMax * sin(C*atan(B*( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]))))
        @NLexpression(mdl, F_yf[i = 1:N], -2*FMax * (((z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2]))
        @NLexpression(mdl, F_yr[i = 1:N], -2*FMax *  ((z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) )


        # @NLconstraint(mdl, [i=1:N], atan( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) <=  max_alpha/180*pi)
        # @NLconstraint(mdl, [i=1:N], atan( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) >= -max_alpha/180*pi)
        # @NLconstraint(mdl, [i=1:N], atan( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2] <=  max_alpha/180*pi)
        # @NLconstraint(mdl, [i=1:N], atan( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2] >= -max_alpha/180*pi)
        @NLconstraint(mdl, [i=1:N], ( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) <=  max_alpha/180*pi)
        @NLconstraint(mdl, [i=1:N], ( (z_Ol[i,2]-l_B*z_Ol[i,3]) / z_Ol[i,1]) >= -max_alpha/180*pi)
        @NLconstraint(mdl, [i=1:N], ( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2] <=  max_alpha/180*pi)
        @NLconstraint(mdl, [i=1:N], ( (z_Ol[i,2]+l_A*z_Ol[i,3]) / z_Ol[i,1] ) - u_Ol[i,2] >= -max_alpha/180*pi)



        
        
        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + dt*(z_Ol[i,3]*z_Ol[i,2]+u_Ol[i,1]-1/mass*(F_yf[i]*u_Ol[i,2]))) #v_x
        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + dt*(-z_Ol[i,3]*z_Ol[i,1]+1/mass*(F_yf[i]*1+F_yr[i]))) #v_y
        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + dt*(1/I_z*(l_A*F_yf[i]*1-l_B*F_yr[i])))    # psi_dot  
        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,5]  == z_Ol[i,5] + dt*(z_Ol[i,1]*z_Ol[i,4]+z_Ol[i,2]))#e_y

        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + c_Vx[1]*z_Ol[i,2]*z_Ol[i,3] + c_Vx[2]*z_Ol[i,1] + c_Vx[3]*u_Ol[i,1]) #v_x
        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + c_Vy[2]*z_Ol[i,1]*z_Ol[i,3] + c_Vy[1]*z_Ol[i,2]/z_Ol[i,1] + c_Vy[3]*z_Ol[i,3]/z_Ol[i,1] + c_Vy[4]*u_Ol[i,2]) #v_y
        # @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + c_Psi[1]*z_Ol[i,3]/z_Ol[i,1] + c_Psi[2]*z_Ol[i,2]/z_Ol[i,1] +c_Psi[3]*u_Ol[i,2] )# psi_dot
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + dt*(z_Ol[i,3]*z_Ol[i,2]+u_Ol[i,1]-1/mass*(F_yf[i]*sin(u_Ol[i,2])))) #v_x
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + dt*(-z_Ol[i,3]*z_Ol[i,1]+1/mass*(F_yf[i]*cos(u_Ol[i,2])+F_yr[i]))) #v_y
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + dt*(1/I_z*(l_A*F_yf[i]*cos(u_Ol[i,2])-l_B*F_yr[i])))    # psi_dot        
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,4]  == z_Ol[i,4] + dt*(z_Ol[i,3]-dsdt[i]*c[i]) )            # epsi
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,5]  == z_Ol[i,5] + dt*(z_Ol[i,1]*sin(z_Ol[i,4])+z_Ol[i,2]*cos(z_Ol[i,4])))#e_y
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,6]  == z_Ol[i,6] + dt*dsdt[i]) #s



        #define expressions for cost
        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum(QderivZ[j] * sum((z_Ol[i,j] - z_Ol[i + 1,j]) ^ 2 for i = 1:N) for j = 1:6) +
                                      sum(QderivU[j]*((uCurr[j]-u_Ol[1,j])^2+sum((u_Ol[i,j]-u_Ol[i+1,j])^2 for i=1:N-1)) for j=1:2))
        
        # Lane cost
        # ---------------------------------
        @NLexpression(mdl, laneCost, Q_lane*sum(10.0*eps[i]+7.0*eps[i]^2 for i=1:2))


        # soft constraint for max velocity
        @NLexpression(mdl, velocityCost, Q_velocity*(10.0*eps[3]+7.0*eps[3]^2)  )

        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*sum(R[j]*sum((u_Ol[i,j]-u_Ref[i,j])^2 for i=1:N) for j=1:2))


        # Terminal constraints (soft)
        # ---------------------------------
        @NLexpression(mdl, constZTerm, sum(Q_term[j]* (sum(ssInfOn[k]*lambda[k]*sum(coeffTermConst[i,k,j]*z_Ol[N+1,6]^(order+1-i) for i=1:order+1) for k=1:n_oldTraj)-z_Ol[N+1,j])^2 for j=1:5))
        #basic idea:
        #@NLexpression(mdl, constZTerm, sum{Q_term[j]*(sum{coeffTermConst[i,1,j]*z_Ol[N+1,1]^(order+1-i),i=1:order+1}-z_Ol[N+1,j+1])^2,j=1:3})


        # Terminal cost
        # ---------------------------------
        @NLexpression(mdl, costZTerm,  sum(ssInfOn[k]*lambda[k]*sum(coeffTermCost[i,k]*z_Ol[N+1,6]^(order+1-i) for i=1:order+1) for k=1:n_oldTraj))
        #basic idea    
        #@NLexpression(mdl, costZTerm, Q_cost*sum{coeffTermCost[i,1]*z_Ol[N+1,1]^(order+1-i),i=1:order+1})


        #put cost on z (actually should put cost only on z before finishing the lap)
        #@NLexpression(mdl, costZ, 1 + (costZ_h-1) * (0.5+0.5*tanh(50*(z_Ol[1,N+1]+s_start-s_target))))
        @NLexpression(mdl, costZ, Q_cost*sum(1 for i=1:N+1))

        ## Cost to avoid obstacle. increases when car is near obstacle currently implemented as : a *1/(0.1+cost)
        # ---------------------------------
        @NLexpression(mdl, costObstacle, sum(Q_obstacleNumer*1/(0.01+(Q_obstacle* (( (z_Ol[i,6]-sCoord_obst[i,1])/rs )^2 + ( (z_Ol[i,5]-sCoord_obst[i,2])/ry) ^2 - 1))^4)+
                                             Q_obstacleNumer*3/(0.01+0.6*(((z_Ol[i,6]-sCoord_obst[i,1])/rs)^2+((z_Ol[i,5]-sCoord_obst[i,2])/ry)^2)) for i=1:N+1))

        #Setting objective formulation
        # ---------------------------------
        @NLobjective(mdl, Min, costZ + costZTerm + constZTerm + derivCost + controlCost + laneCost+ velocityCost + costObstacle)





        m.mdl = mdl
        m.coeff = coeff # curvature coefficients
        m.coeffTermCost = coeffTermCost #coefficients of terminal Cost
        m.coeffTermConst = coeffTermConst   # coefficients for terminal constraints
        m.uCurr = uCurr #last applied input
        m.sCoord_obst = sCoord_obst
        m.u_Ol = u_Ol
        m.z_Ol = z_Ol
        m.lambda = lambda
        m.ssInfOn = ssInfOn
        m.eps = eps
        m.z0  = z0
        m.derivCost= derivCost
        m.laneCost = laneCost
        m.velocityCost = velocityCost
        m.controlCost = controlCost
        m.constZTerm = constZTerm 
        m.costZTerm = costZTerm
        m.costZ = costZ       
        m.costObstacle = costObstacle

        return m
    end
end
