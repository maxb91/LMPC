# Variable definitions
# mdl.z_Ol[i,j] = z_OpenLoop, open loop prediction of the state, i = state, j = step

# States in path following mode:
# i = 1 -> s
# i = 2 -> ey
# i = 3 -> epsi
# i = 4 -> v

# States in LMPC and system ID mode:
# i = 1 -> xDot
# i = 2 -> yDot
# i = 3 -> psiDot
# i = 4 -> ePsi
# i = 5 -> eY
# i = 6 -> s

function solveMpcProblem(mdl::MpcModel,mpcSol::MpcSol,mpcCoeff::MpcCoeff,mpcParams::MpcParams,trackCoeff::TrackCoeff,lapStatus::LapStatus,posInfo::PosInfo,modelParams::ModelParams,zCurr::Array{Float64},uCurr::Array{Float64})

    println("--------------- MPC START -----------------------------------------------")
    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    N               = mpcParams.N
    Q               = mpcParams.Q
    Q_term          = mpcParams.Q_term
    R               = mpcParams.R
    coeffTermCost   = mpcCoeff.coeffCost::Array{Float64,2}
    coeffTermConst  = mpcCoeff.coeffConst::Array{Float64,3}
    order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
    s_target        = posInfo.s_target
    ey_max          = trackCoeff.width/2

    QderivZ         = mpcParams.QderivZ::Array{Float64,1}
    QderivU         = mpcParams.QderivU::Array{Float64,1}

    v_ref           = mpcParams.vPathFollowing

    dt = modelParams.dt

    sol_status::Symbol
    sol_u::Array{Float64,2}
    sol_z::Array{Float64,2}

    # println("************************************** MPC SOLVER **************************************")
    # println("zCurr    = $(zCurr')")
    # println("s_start  = $s_start")
    # println("s_target = $s_target")
    # println("s_total  = $((zCurr[1]+s_start)%s_target)")

    # Create function-specific parameters
    z_Ref::Array{Float64,2}
    z_Ref           = cat(2,s_target*ones(N+1,1),zeros(N+1,2),v_ref*ones(N+1,1))       # Reference trajectory: path following -> stay on line and keep constant velocity
    u_Ref           = zeros(N,2)

    # Update current initial condition, curvature and System ID coefficients
    #z_set = [zCurr[6],zCurr[5],zCurr[4],zCurr[1],zCurr[6],zCurr[5]]
    setvalue(mdl.z0,zCurr)
    setvalue(mdl.coeff,coeffCurvature)
    setvalue(mdl.c_Vx,mpcCoeff.c_Vx)
    setvalue(mdl.c_Vy,mpcCoeff.c_Vy)
    setvalue(mdl.c_Psi,mpcCoeff.c_Psi)

    println("mdl.z0    = $(getvalue(mdl.z0))")
    println("coeffCurvature = $(getvalue(mdl.coeff))")
    println("mdl.c_Vx  = $(getvalue(mdl.c_Vx))")
    println("mdl.c_Vy  = $(getvalue(mdl.c_Vy))")
    println("mdl.c_Psi = $(getvalue(mdl.c_Psi))")

    @NLexpression(mdl.mdl, costZTerm,   0)
    @NLexpression(mdl.mdl, constZTerm,  0)
    @NLexpression(mdl.mdl, derivCost,   0)
    @NLexpression(mdl.mdl, controlCost, 0)
    @NLexpression(mdl.mdl, laneCost,    0)

    for i=1:N
        #setlowerbound(mdl.u_Ol[i,1],dt/mpcCoeff.c_Vx[3]*0.1)
        #setupperbound(mdl.u_Ol[i,1],dt/mpcCoeff.c_Vx[3]*1.2)
        #setlowerbound(mdl.u_Ol[i,2],-Inf)
        #setupperbound(mdl.u_Ol[i,2],Inf)
    end

    # # Derivative cost
    # # ---------------------------------
    # @NLexpression(mdl.mdl, derivCost, sum{QderivZ[j]*((zCurr[j]-mdl.z_Ol[1,j])^2+sum{(mdl.z_Ol[i,j]-mdl.z_Ol[i+1,j])^2,i=1:N}),j=1:6} +
    #                                  sum{QderivU[j]*((uCurr[j]-mdl.u_Ol[1,j])^2+sum{(mdl.u_Ol[i,j]-mdl.u_Ol[i+1,j])^2,i=1:N-1}),j=1:2})

    # # Lane cost
    # # ---------------------------------
    # @NLexpression(mdl.mdl, laneCost, 10*sum{mdl.z_Ol[i,5]^2*((0.5+0.5*tanh(10*(mdl.z_Ol[i,5]-ey_max))) + (0.5-0.5*tanh(10*(mdl.z_Ol[i,5]+ey_max)))),i=1:N+1})

    # # Control Input cost
    # # ---------------------------------
    # @NLexpression(mdl.mdl, controlCost, 0.5*sum{R[j]*sum{(mdl.u_Ol[i,j]-u_Ref[i,j])^2,i=1:N},j=1:2})

    # # Terminal constraints (soft), starting from 2nd lap
    # # ---------------------------------
    # if lapStatus.currentLap > 2    # if at least in the 3rd lap
    #     @NLexpression(mdl.mdl, constZTerm, (sum{Q_term[j]*(mdl.ParInt[1]*sum{coeffTermConst[i,1,j]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1}+
    #                                     (1-mdl.ParInt[1])*sum{coeffTermConst[i,2,j]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1}-mdl.z_Ol[N+1,j])^2,j=1:5}))
    # elseif lapStatus.currentLap == 2        # if in the 2nd lap
    #     @NLexpression(mdl.mdl, constZTerm, sum{Q_term[j]*(sum{coeffTermConst[i,1,j]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1}-mdl.z_Ol[N+1,j])^2,j=1:5})
    # end

    # # Terminal cost
    # # ---------------------------------
    # # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
    # if lapStatus.currentLap > 2     # if at least in the 3rd lap
    #     @NLexpression(mdl.mdl, costZTerm, mdl.ParInt[1]*sum{coeffTermCost[i,1]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1}+
    #                               (1-mdl.ParInt[1])*sum{coeffTermCost[i,2]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1})
    # elseif lapStatus.currentLap == 2         # if we're in the second second lap
    #     @NLexpression(mdl.mdl, costZTerm, sum{coeffTermCost[i,1]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1})
    # end
    
    # @NLobjective(mdl.mdl, Min, costZTerm + constZTerm + derivCost + controlCost + laneCost)

    @NLexpression(mdl.mdl, costZ, sum{(mdl.z_Ol[j,5])^2+(mdl.z_Ol[j,1]-0.7)^2,j=2:N+1}+0.1*sum{(mdl.u_Ol[j,1])^2,j=1:N}+0.1*sum{(mdl.u_Ol[j,2])^2,j=1:N})    # Follow trajectory (only minimize eY and v = 0.2m/s)
    @NLobjective(mdl.mdl, Min, costZ)

    # println("Model formulation:")
    # println(mdl.mdl)
    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)
    println("Solved")

    println("z_Ol      = $(getvalue(mdl.z_Ol))")
    println("u_Ol      = $(getvalue(mdl.u_Ol))")
    sol_u       = getvalue(mdl.u_Ol)
    sol_z       = getvalue(mdl.z_Ol)
    println("Predicting until z = $(sol_z[end,1])")
    #println("curvature = $(getvalue(mdl.c))")

    # COST PRINTS: ********************************************************
    # println("coeff: $(getvalue(mdl.coeff))")
    # println("z0: $(getvalue(mdl.z0))")
    # println("Solution status: $sol_status")
    # println("Objective value: $(getobjectivevalue(mdl.mdl))")
    # println("Control Cost: $(getvalue(controlCost))")
    # println("CostZ:        $(getvalue(costZ))")
    # println("DerivCost:    $(getvalue(derivCost))")
    # println("LaneCost:     $(getvalue(laneCost))")
    # println("costZTerm:    $(getvalue(costZTerm))")
    # println("constZTerm:   $(getvalue(constZTerm))")

    # println("cost_ey:      $(0.5*sum(sol_z[2,:].^2)*Q[2])")
    # println("cost_ePsi:    $(0.5*sum(sol_z[3,:].^2)*Q[3])")
    # println("cost_V:       $(0.5*sum((sol_z[4,:]-z_Ref[:,4]').^2)*Q[4])")

    #println("z:")
    #println(getvalue(mdl.z_Ol))
    #println("u:")
    #println(getvalue(mdl.u_Ol))
    #mpcSol      = MpcSol(sol_u[1,1],sol_u[2,1],sol_status,getvalue(mdl.u_Ol),getvalue(mdl.z_Ol),[getvalue(costZ),getvalue(costZTerm),getvalue(constZTerm),getvalue(derivCost),getvalue(controlCost),getvalue(laneCost)])
    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    mpcSol.cost = zeros(6)
    mpcSol.cost = [0,getvalue(costZTerm),getvalue(constZTerm),getvalue(derivCost),getvalue(controlCost),getvalue(laneCost)]
    #mpcSol = MpcSol(sol_u[1,1],sol_u[2,1]) # Fast version without logging
    
    println("--------------- MPC END ------------------------------------------------")
    nothing
end

function solveMpcProblem_pathFollow(mdl::MpcModel_pF,mpcSol::MpcSol,mpcParams::MpcParams,trackCoeff::TrackCoeff,posInfo::PosInfo,modelParams::ModelParams,zCurr::Array{Float64},uCurr::Array{Float64})

    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    N               = mpcParams.N
    Q               = mpcParams.Q
    R               = mpcParams.R
    s_target        = posInfo.s_target

    QderivZ         = mpcParams.QderivZ::Array{Float64,1}
    QderivU         = mpcParams.QderivU::Array{Float64,1}

    v_ref           = mpcParams.vPathFollowing

    sol_status::Symbol
    sol_u::Array{Float64,2}
    sol_z::Array{Float64,2}

    # Create function-specific parameters
    z_Ref::Array{Float64,2}
    z_Ref           = cat(2,s_target*ones(N+1,1),zeros(N+1,2),v_ref*ones(N+1,1))       # Reference trajectory: path following -> stay on line and keep constant velocity
    u_Ref           = zeros(N,2)

    # Update current initial condition, curvature and System ID coefficients
    setvalue(mdl.z0,zCurr)
    setvalue(mdl.coeff,coeffCurvature)

    # Derivative cost
    # ---------------------------------
    @NLexpression(mdl.mdl, derivCost, sum{QderivZ[j]*((zCurr[j]-mdl.z_Ol[1,j])^2+sum{(mdl.z_Ol[i,j]-mdl.z_Ol[i+1,j])^2,i=1:N}),j=1:4} +
                                      sum{QderivU[j]*((uCurr[j]-mdl.u_Ol[1,j])^2+sum{(mdl.u_Ol[i,j]-mdl.u_Ol[i+1,j])^2,i=1:N-1}),j=1:2})

    # Control Input cost
    # ---------------------------------
    @NLexpression(mdl.mdl, controlCost, 0.5*sum{R[j]*sum{(mdl.u_Ol[i,j]-u_Ref[i,j])^2,i=1:N},j=1:2})

    # State cost
    # ---------------------------------
    @NLexpression(mdl.mdl, costZ, 0.5*sum{Q[i]*sum{(mdl.z_Ol[j,i]-z_Ref[j,i])^2,j=2:N+1},i=1:4})    # Follow trajectory


    @NLobjective(mdl.mdl, Min, costZ + derivCost + controlCost)

    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)
    sol_u       = getvalue(mdl.u_Ol)
    sol_z       = getvalue(mdl.z_Ol)

    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    #mpcSol.cost = zeros(6)
    mpcSol.cost = [getvalue(costZ),0,0,getvalue(derivCost),getvalue(controlCost),0]

    nothing
end
