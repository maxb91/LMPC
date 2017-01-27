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

        dt         = modelParams.dt
        L_a        = modelParams.l_A
        L_b        = modelParams.l_B
        u_lb       = modelParams.u_lb
        u_ub       = modelParams.u_ub
        z_lb       = modelParams.z_lb
        z_ub       = modelParams.z_ub
        s_obst     = obstacle.s_obstacle
        sy_obst    = obstacle.sy_obstacle
        @show rs = obstacle.rs
        ry = obstacle.ry
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
        s_target = 100 # the weight for s in pathfollowing is set to zero, so this value is not used
        local z_Ref::Array{Float64,2}
        z_Ref           = cat(2,s_target*ones(N+1,1),zeros(N+1,2),v_ref*ones(N+1,1))     # Reference trajectory: path following -> stay on line and keep constant velocity
        u_Ref           = zeros(N,2)
	    z_Init= zeros(N+1,4)
        z_Init[:,4] = 0.6*ones(N+1)
	    v_max = modelParams.v_max
        max_alpha =modelParams.max_alpha
        #########################################################
        mdl = Model(solver = IpoptSolver(print_level=0))#mu_strategy=adaptive,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))
        #########################################################
        

        z_Ol = @variable( mdl, [i =1:(N+1),j = 1:4], lowerbound = modelParams.z_lb[i,j], upperbound = modelParams.z_ub[i,j], start = z_Init[i,j])      # z = s, ey, epsi, v
        u_Ol = @variable( mdl, [i=1:N,j=1:2], lowerbound = modelParams.u_lb[i,j], upperbound = modelParams.u_ub[i,j], start = 0)          # overwrtie dim of in classes.jl?
        eps  = @variable( mdl, lowerbound = 0, start = 0)

        @NLparameter(mdl, z0[i=1:4] == z_Init[1,i])
        @NLconstraint(mdl, [i=1:4], z_Ol[1,i] == z0[i])

        #@constraint(mdl, lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]== 1)

        
        @NLconstraint(mdl,[i = 1:(N+1)], z_Ol[i,4] <= v_max +eps )

        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])
        @NLparameter(mdl, uCurr[i=1:2] == 0)
        @NLparameter(mdl, sCoord_obst[j=1:(N+1),i=1:2] == 100)
  

        @NLexpression(mdl, c[i = 1:N],    coeff[1]*z_Ol[i,1]^4+coeff[2]*z_Ol[i,1]^3+coeff[3]*z_Ol[i,1]^2+coeff[4]*z_Ol[i,1]+coeff[5])
        #@NLexpression(mdl, c[i = 1:N],    sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1]) 
        @NLexpression(mdl, bta[i = 1:N],  atan( L_b / (L_a + L_b) * tan(u_Ol[i,2]) ) )

        @NLexpression(mdl, dsdt[i = 1:N], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))
        # System dynamics
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + dt*dsdt[i]  )
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + dt*z_Ol[i,4]*sin(z_Ol[i,3]+bta[i])  )                     # ey
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + dt*(z_Ol[i,4]/L_b*sin(bta[i])-dsdt[i]*c[i])  )            # epsi
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,4]  == z_Ol[i,4] + dt*(u_Ol[i,1]))#- 0.23*abs(z_Ol[i,4]) * z_Ol[i,4]))#0.63  # v
     


        #define expressions for cost
        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum(QderivZ[j] * sum((z_Ol[i,j] - z_Ol[i + 1,j]) ^ 2 for i = 1:N) for j = 1:4) +
                                          sum(QderivU[j]*((uCurr[j]-u_Ol[1,j])^2+sum((u_Ol[i,j]-u_Ol[i+1,j])^2 for i=1:N-1)) for j=1:2))
        derivCost= derivCost
       

        # soft constraint for max velocity
        @NLexpression(mdl, velocityCost, Q_velocity*(10.0*eps+7.0*eps^2)  )
        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*sum(R[j]*sum((u_Ol[i,j]-u_Ref[i,j])^2 for i=1:N) for j=1:2))

         # ---------------------------------
        # if we're in the first lap, just do path following
        @NLexpression(mdl, costPath, 0.5*sum(Q[i]*sum((z_Ol[j,i]-z_Ref[j,i])^2for j=1:N+1) for i=1:4))    # Follow trajectory

        ## Cost to avoid obstacle. increases when car is near obstacle currently implemented as : a *1/(0.1+cost)
        @NLexpression(mdl, costObstacle, sum(Q_obstacleNumer*1/(0.01+(Q_obstacle* (( (z_Ol[i,1]-sCoord_obst[i,1])/rs )^2 + ( (z_Ol[i,2]-sCoord_obst[i,2])/ry) ^2 - 1))^4)+
                                             Q_obstacleNumer*3/(0.01+0.6*(((z_Ol[i,1]-sCoord_obst[i,1])/rs)^2+((z_Ol[i,2]-sCoord_obst[i,2])/ry)^2)) for i=1:N+1))

        #objective formulation, minimize the sum of all parts of the objective
        @NLobjective(mdl, Min, costPath + derivCost + controlCost + velocityCost + costObstacle)
        m.mdl = mdl

        # Update model values
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
        obstacle::classes.Obstacle, n_oldTraj)

        m = new()

        dt        = modelParams.dt
        L_a       = modelParams.l_A
        L_b       = modelParams.l_B
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
        coeffCost   = mpcCoeff.coeffCost[1,:,:]::Array{Float64,2}
        coeffConst  = mpcCoeff.coeffConst[1,:,:,:]::Array{Float64,3}
        order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
        QderivZ         = mpcParams.QderivZ::Array{Float64,1}
        QderivU         = mpcParams.QderivU::Array{Float64,1}
        v_ref           = mpcParams.vPathFollowing

        N           = mpcParams.N
        ey_max      = trackCoeff.width/2
        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation
        # Create function-specific parameters
        u_Ref           = zeros(N,2)

        z_Init       = zeros(N+1,4)
        z_Init[:,4]  = v_ref*ones(N+1)

        v_max = modelParams.v_max
        max_alpha =modelParams.max_alpha
        
        #########################################################
        mdl = Model(solver = IpoptSolver(print_level=0))#mu_strategy=adaptive,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))
        #########################################################

        z_Ol = @variable( mdl, [i =1:(N+1),j = 1:4], lowerbound = z_lb[i,j], upperbound = z_ub[i,j], start = z_Init[i,j])      # z = s, ey, epsi, v
        u_Ol = @variable( mdl, [i=1:N,j=1:2], lowerbound = u_lb[i,j], upperbound = u_ub[i,j], start =0)          # overwrtie dim of in classes.jl?
        lambda = @variable( mdl, [1:n_oldTraj],lowerbound = 0, upperbound = 1, start = 0)
        eps = @variable( mdl,[1:3], lowerbound = 0, start = 0)

        @NLparameter(mdl,ssInfOn[1:n_oldTraj]== 1)
        @NLparameter(mdl, z0[i=1:4] == z_Init[1,i])
        @NLconstraint(mdl, [i=1:4], z_Ol[1,i] == z0[i])

        #@constraint(mdl, lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]== 1)
        @NLconstraint(mdl, sum(lambda[j] for j=1:n_oldTraj)== 1)
        @NLconstraint(mdl,[i = 1:(N+1)], z_Ol[i,2] <=  trackCoeff.width/2 + eps[1] )
        @NLconstraint(mdl,[i = 1:(N+1)], z_Ol[i,2] >= -trackCoeff.width/2 - eps[2] )

        @NLconstraint(mdl,[i = 1:(N+1)], z_Ol[i,4] <= v_max +eps[3] )

        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])
        @NLparameter(mdl, uCurr[i=1:2] == 0)
        @NLparameter(mdl, sCoord_obst[j=1:(N+1),i=1:2] == 100)
        @NLparameter(mdl, coeffTermConst[i=1:order+1,k=1:n_oldTraj,j=1:3] == coeffConst[i,k,j])
        @NLparameter(mdl, coeffTermCost[i=1:order+1,k=1:n_oldTraj] == coeffCost[i,k])

        @NLexpression(mdl, c[i = 1:N],    coeff[1]*z_Ol[i,1]^4+coeff[2]*z_Ol[i,1]^3+coeff[3]*z_Ol[i,1]^2+coeff[4]*z_Ol[i,1]+coeff[5])
        #@NLexpression(mdl, c[i = 1:N],    sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1]) 

        @NLexpression(mdl, bta[i = 1:N],  atan( L_b / (L_a + L_b) * tan(u_Ol[i,2]) ) )
        @NLexpression(mdl, dsdt[i = 1:N], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))
        # System dynamics
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + dt*dsdt[i]  )
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + dt*z_Ol[i,4]*sin(z_Ol[i,3]+bta[i])  )                     # ey
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + dt*(z_Ol[i,4]/L_b*sin(bta[i])-dsdt[i]*c[i])  )            # epsi
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,4]  == z_Ol[i,4] + dt*(u_Ol[i,1]))#- 0.23*abs(z_Ol[i,4]) * z_Ol[i,4]))#0.63  # v
     
        # @NLconstraint(mdl,[i =1:(N+1)],( (z_Ol[i,1]-sCoord_obst[i,1])/rs )^2 + ( (z_Ol[i,2]-sCoord_obst[i,2])/ry) ^2 - 1>=0)

        #define expressions for cost
        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum(QderivZ[j] * sum((z_Ol[i,j] - z_Ol[i + 1,j]) ^ 2 for i = 1:N) for j = 1:4) +  #(z0[j]-z_Ol[1,j])^2
                                          sum(QderivU[j]*((uCurr[j]-u_Ol[1,j])^2+sum((u_Ol[i,j]-u_Ol[i+1,j])^2 for i=1:N-1)) for j=1:2))

        # Lane cost
        # ---------------------------------
        #@NLexpression(mdl, laneCost, Q_lane*sum{z_Ol[i,2]^2*((0.5+0.5*tanh(35*(z_Ol[i,2]-ey_max-0.09))) + (0.5-0.5*tanh(35*(z_Ol[i,2]+ey_max+0.09)))),i=1:N+1})
        @NLexpression(mdl, laneCost, Q_lane*sum(10.0*eps[i]+7.0*eps[i]^2 for i=1:2))



        # soft constraint for max velocity
        @NLexpression(mdl, velocityCost, Q_velocity*(10.0*eps[3]+7.0*eps[3]^2)  )

        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*sum(R[j]*sum((u_Ol[i,j]-u_Ref[i,j])^2 for i=1:N) for j=1:2))


        # Terminal constraints (soft), starting from 2nd lap

        # ---------------------------------                   
        @NLexpression(mdl, constZTerm, sum(Q_term[j]* (sum(ssInfOn[k]*lambda[k]*sum(coeffTermConst[i,k,j]*z_Ol[N+1,1]^(order+1-i) for i=1:order+1) for k=1:n_oldTraj)-z_Ol[N+1,j+1])^2 for j=1:3))
       #basic idea:
        #@NLexpression(mdl, constZTerm, sum{Q_term[j]*(sum{coeffTermConst[i,1,j]*z_Ol[N+1,1]^(order+1-i),i=1:order+1}-z_Ol[N+1,j+1])^2,j=1:3})

        # Terminal cost
        # ---------------------------------
        # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
        @NLexpression(mdl, costZTerm,  sum(ssInfOn[k]*lambda[k]*sum(coeffTermCost[i,k]*z_Ol[N+1,1]^(order+1-i) for i=1:order+1) for k=1:n_oldTraj))
        #basic idea    
        #@NLexpression(mdl, costZTerm, Q_cost*sum{coeffTermCost[i,1]*z_Ol[N+1,1]^(order+1-i),i=1:order+1})

        # State cost
        # ---------------------------------
        @NLexpression(mdl, costZ, Q_cost*sum(1 for i=1:N+1))

        ## Cost to avoid obstacle. increases when car is near obstacle currently implemented as : a *1/(0.1+cost)
        # @NLexpression(mdl, costObstacle, sum(((N+1.2-0.2*i)/(N+1))*(Q_obstacleNumer*1/(0.01+(Q_obstacle* (( (z_Ol[i,1]-sCoord_obst[i,1])/rs )^2 + ( (z_Ol[i,2]-sCoord_obst[i,2])/ry) ^2 - 1))^4)+
        #                                      Q_obstacleNumer*3/(0.01+0.6*(((z_Ol[i,1]-sCoord_obst[i,1])/rs)^2+((z_Ol[i,2]-sCoord_obst[i,2])/ry)^2))) for i=1:N+1))
        @NLexpression(mdl, costObstacle, sum(Q_obstacleNumer*1/(0.01+(Q_obstacle* (( (z_Ol[i,1]-sCoord_obst[i,1])/rs )^2 + ( (z_Ol[i,2]-sCoord_obst[i,2])/ry) ^2 - 1))^4)+
                                             Q_obstacleNumer*3/(0.01+0.6*(((z_Ol[i,1]-sCoord_obst[i,1])/rs)^2+((z_Ol[i,2]-sCoord_obst[i,2])/ry)^2)) for i=1:N+1))

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
