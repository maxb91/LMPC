# Variable definitions
# mdl.z_Ol[i,j] = z_OpenLoop, open loop prediction of the state, i = state, j = step

# States:
# i = 1 -> s
# i = 2 -> ey
# i = 3 -> epsi
# i = 4 -> v

function solveMpcProblem!(mdl::classes.MpcModel,mpcSol::classes.MpcSol,mpcCoeff::classes.MpcCoeff,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,lapStatus::classes.LapStatus,posInfo::classes.PosInfo,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64}, obstacle::classes.Obstacle, iter::Int64)
    #tic()
    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    N               = mpcParams.N
    Q               = mpcParams.Q #Cost of states just for path following
    Q_term          = mpcParams.Q_term
    Q_cost          = mpcParams.Q_cost
    R               = mpcParams.R # cost for control is always used but curently 0
    coeffTermCost   = mpcCoeff.coeffCost::Array{Float64,2}
    coeffTermConst  = mpcCoeff.coeffConst::Array{Float64,3}
    order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
    s_start         = posInfo.s_start
    s_target        = posInfo.s_target
    ey_max          = trackCoeff.width/2

    s_obst     = obstacle.s_obstacle[iter]
    sy_obst    = obstacle.sy_obstacle[iter]
    rs         = obstacle.rs
    ry         = obstacle.ry


    QderivZ         = mpcParams.QderivZ::Array{Float64,1}
    QderivU         = mpcParams.QderivU::Array{Float64,1}

    v_ref           = mpcParams.vPathFollowing

   local sol_u::Array{Float64,2} 
   local sol_z::Array{Float64,2}

    # println("************************************** MPC SOLVER **************************************")
    # println("zCurr    = $(zCurr')")
    # println("s_start  = $s_start")
    # println("s_target = $s_target")
    # println("s_total  = $((zCurr[1]+s_start)%s_target)")

    # Create function-specific parameters
    local z_Ref::Array{Float64,2}
    z_Ref           = cat(2,s_target*ones(N+1,1),zeros(N+1,2),v_ref*ones(N+1,1))     # Reference trajectory: path following -> stay on line and keep constant velocity
    u_Ref           = zeros(N,2)

    # Update current initial condition
    setvalue(mdl.z0,zCurr')
    #setvalue(mdl.z0[1],zCurr[1])
    #setvalue(mdl.s_startC, s_start)

    # Update curvature
    setvalue(mdl.coeff,coeffCurvature)
    #println("z0 = $(getvalue(mdl.z0))")
    #println("coeffCurvature = $(getvalue(mdl.coeff))")

    @NLexpression(mdl.mdl, costZ,       0)
    @NLexpression(mdl.mdl, costZTerm,   0)
    @NLexpression(mdl.mdl, constZTerm,  0)
    @NLexpression(mdl.mdl, derivCost,   0)
    @NLexpression(mdl.mdl, laneCost,    0)
    @NLexpression(mdl.mdl, costObstacle,    0)

    #n_poly_curv = trackCoeff.nPolyCurvature 
        #@NLexpression(mdl.mdl, mdl.c[i = 1:N],    sum{mdl.coeff[j]*(mdl.z_Ol[i,1]-s_start)^(n_poly_curv-j+1),j=1:n_poly_curv} + mdl.coeff[n_poly_curv+1])
    # Derivative cost
    # ---------------------------------
    @NLexpression(mdl.mdl, derivCost, sum{QderivZ[j]*((zCurr[j]-mdl.z_Ol[1,j])^2+sum{(mdl.z_Ol[i,j]-mdl.z_Ol[i+1,j])^2,i=1:N}),j=1:4} +
                                      sum{QderivU[j]*((uCurr[j]-mdl.u_Ol[1,j])^2+sum{(mdl.u_Ol[i,j]-mdl.u_Ol[i+1,j])^2,i=1:N-1}),j=1:2})

    # Lane cost
    # ---------------------------------
    @NLexpression(mdl.mdl, laneCost, 1*sum{mdl.z_Ol[i,2]^2*((0.5+0.5*tanh(30*(mdl.z_Ol[i,2]-ey_max-0.1))) + (0.5-0.5*tanh(30*(mdl.z_Ol[i,2]+ey_max+0.1)))),i=1:N+1})

    # Control Input cost
    # ---------------------------------
    @NLexpression(mdl.mdl, controlCost, 0.5*sum{R[j]*sum{(mdl.u_Ol[i,j]-u_Ref[i,j])^2,i=1:N},j=1:2})



    # Terminal constraints (soft), starting from 2nd lap
    #constraints force trajectory to end up on ss
    # ---------------------------------
    if lapStatus.currentLap > 2    # if at least in the 3rd lap, as of the third round we have two old trajectories between which we can interpolate
        @NLexpression(mdl.mdl, constZTerm, (sum{Q_term[j]*( mdl.lambda[1]*sum{coeffTermConst[i,1,j]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}+
                                                            mdl.lambda[2]*sum{coeffTermConst[i,2,j]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}+
                                                            mdl.lambda[3]*sum{coeffTermConst[i,3,j]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}+
                                                            mdl.lambda[4]*sum{coeffTermConst[i,4,j]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}+
                                                            mdl.lambda[5]*sum{coeffTermConst[i,5,j]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}-
                                                            mdl.z_Ol[N+1,j+1])^2,j=1:3}))
    elseif lapStatus.currentLap == 2        # if in the 2nd lap
        @NLexpression(mdl.mdl, constZTerm, sum{Q_term[j]*(sum{coeffTermConst[i,1,j]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}-mdl.z_Ol[N+1,j+1])^2,j=1:3})
    end

    # Terminal cost
    # ---------------------------------
    # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
    if lapStatus.currentLap > 2     # if at least in the 3rd lap
        @NLexpression(mdl.mdl, costZTerm, Q_cost*(  mdl.lambda[1]*sum{coeffTermCost[i,1]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}+
                                                    mdl.lambda[2]*sum{coeffTermCost[i,2]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}+
                                                    mdl.lambda[3]*sum{coeffTermCost[i,3]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}+
                                                    mdl.lambda[4]*sum{coeffTermCost[i,4]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}+
                                                    mdl.lambda[5]*sum{coeffTermCost[i,5]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1}))
    elseif lapStatus.currentLap == 2         # if we're in the second second lap
        @NLexpression(mdl.mdl, costZTerm, Q_cost*sum{coeffTermCost[i,1]*mdl.z_Ol[N+1,1]^(order+1-i),i=1:order+1})
    end

    # State cost
    # ---------------------------------
    if lapStatus.currentLap <= 1      # if we're in the first lap, just do path following
        @NLexpression(mdl.mdl, costZ, 0.5*sum{Q[i]*sum{(mdl.z_Ol[j,i]-z_Ref[j,i])^2,j=2:N+1},i=1:4})    # Follow trajectory

    else        # if we're in another lap, put cost on z (actually should put cost only on z before finishing the lap)
        #@NLexpression(mdl.mdl, costZ_h, 0)          # zero state cost after crossing the finish line
        #@NLexpression(mdl.mdl, costZ, 1 + (costZ_h-1) * (0.5+0.5*tanh(50*(mdl.z_Ol[1,N+1]+s_start-s_target))))
        @NLexpression(mdl.mdl, costZ, 1)
    end
    
    ## Cost to avoid obstacle. increases when car is near obstacle currently implemented as : a *1/(0.1+cost)
    @NLexpression(mdl.mdl, costObstacle, sum{0.02*1/(0.1+ ( (mdl.z_Ol[i,1]-s_obst)/rs )^2 + ( (mdl.z_Ol[i,2]-sy_obst)/ry )^2 - 1 )^2,i=1:N+1})


    #objective formulation, minimize the sum of all parts of the objective
    @NLobjective(mdl.mdl, Min, costZ + costZTerm + constZTerm + derivCost + controlCost + laneCost + costObstacle)

    #println("Model formulation:")
    #println(mdl.mdl)
    # Solve Problem and return solution
    #t_make_mdl= toq()


    tic()
    sol_status  = solve(mdl.mdl)
    ttt= toq()

    sol_u       = getvalue(mdl.u_Ol)
    sol_z       = getvalue(mdl.z_Ol)
    mpcSol.lambda = getvalue(mdl.lambda)

    # c_print = getvalue(mdl.c)
    # println("curvature: $c_print")
    
    #println("Predicting until s = $(sol_z[end,1])") 
    
    


    # COST PRINTS: ********************************************************
    # print("************************COSTS*****************")
    # println("coeff: $(getvalue(mdl.coeff))")
    # println("z0: $(getvalue(mdl.z0))")

    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    mpcSol.cost = zeros(7)

    #tic()
    #mpcSol.cost = [getvalue(costZ);getvalue(costZTerm);getvalue(constZTerm);getvalue(derivCost);getvalue(controlCost);getvalue(laneCost); getvalue(costObstacle)]
    objvl = getobjectivevalue(mdl.mdl)
    #t_get2 =toq()


    if iter%50 == 0 
        println("$iter th iteration****************************")
        println("pure solver time: $ttt")
        #println("get cost time: $t_get2")
        #println("creat model time: $t_make_mdl")
     

    end


    # if lapStatus.currentLap > 100
    #     ss = collect(zCurr[1]:.01:zCurr[1]+0.3)
    #     #p  = [ss.^4 ss.^3 ss.^2 ss.^1 ss.^0] *coeffTermCost[:,:,1]
    #     p  = [ss.^4 ss.^3 ss.^2 ss.^1 ss.^0] *coeffTermConst[:,:,1,1]
    #     plot(ss,p,getvalue(mdl.z_Ol[1,N+1]),getvalue(costZTerm),"o")
    #     grid()
    #     readline()
    # end

    # println(getvalue(mdl.z_Ol))
    # println("==============")
    # println(getvalue(mdl.u_Ol))
    nothing #nothing to return apprently syntax is like this
end
