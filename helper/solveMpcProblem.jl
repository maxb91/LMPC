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
    sol_status::Symbol
    sol_u::Array{Float64,2}
    sol_z::Array{Float64,2}

    N = mpcParams.N

    # println("************************************** MPC SOLVER **************************************")
    # println("zCurr    = $(zCurr')")
    # println("s_start  = $s_start")
    # println("s_target = $s_target")
    # println("s_total  = $((zCurr[1]+s_start)%s_target)")

    # Update current initial condition, curvature and System ID coefficients

    # println("previous values: ")
    # println("z0    = $(getvalue(mdl.z0))")
    # println("uCurr = $(getvalue(mdl.uCurr))")
    # println("c_Vx  = $(getvalue(mdl.c_Vx))")
    # println("c_Vy  = $(getvalue(mdl.c_Vy))")
    # println("c_Psi = $(getvalue(mdl.c_Psi))")
    # println("coeffTermCost  = $(getvalue(mdl.coeffTermCost))")
    # println("coeffTermConst = $(getvalue(mdl.coeffTermConst))")
    # println("trackCoeff     = $(getvalue(mdl.coeff))")
    
    Q_term          = mpcParams.Q_term
    R               = mpcParams.R
    order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
    ey_max          = trackCoeff.width/2

    QderivZ         = mpcParams.QderivZ::Array{Float64,1}
    QderivU         = mpcParams.QderivU::Array{Float64,1}
    Q_term_cost     = mpcParams.Q_term_cost::Float64

    n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation

    setvalue(mdl.z0,zCurr)
    setvalue(mdl.uCurr,uCurr[:])

    setvalue(mdl.c_Vx,mpcCoeff.c_Vx)            # System ID coefficients
    setvalue(mdl.c_Vy,mpcCoeff.c_Vy)
    setvalue(mdl.c_Psi,mpcCoeff.c_Psi)

    setvalue(mdl.coeff,trackCoeff.coeffCurvature)       # Track curvature
    setvalue(mdl.coeffTermCost,mpcCoeff.coeffCost)      # Terminal cost
    setvalue(mdl.coeffTermConst,mpcCoeff.coeffConst)    # Terminal constraints

    # Cost functions

    # Derivative cost
    # ---------------------------------
    #@NLexpression(mdl, derivCost, sum{QderivZ[j]*(sum{(z_Ol[i,j]-z_Ol[i+1,j])^2,i=1:N}),j=1:6} +
    #                                  sum{QderivU[j]*((uCurr[j]-u_Ol[1,j])^2+sum{(u_Ol[i,j]-u_Ol[i+1,j])^2,i=1:N-1}),j=1:2})

    # Lane cost
    # ---------------------------------
    #@NLexpression(mdl, laneCost, 10*sum{z_Ol[i,5]^2*((0.5+0.5*tanh(10*(z_Ol[i,5]-ey_max))) + (0.5-0.5*tanh(10*(z_Ol[i,5]+ey_max)))),i=1:N+1})

    # Control Input cost
    # ---------------------------------
    #@NLexpression(mdl, controlCost, 0.5*sum{R[j]*sum{(u_Ol[i,j])^2,i=1:N},j=1:2})

    # Terminal constraints (soft), starting from 2nd lap
    # ---------------------------------

    @NLexpression(mdl.mdl, constZTerm, sum{Q_term[j]*(mdl.ParInt*sum{mdl.coeffTermConst[i,1,j]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1}+
                                        (1-mdl.ParInt)*sum{mdl.coeffTermConst[i,2,j]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1}-mdl.z_Ol[N+1,j])^2,j=1:5})
    
    # Terminal cost
    # ---------------------------------
    # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
    
    @NLexpression(mdl.mdl, costZTerm, Q_term_cost*(mdl.ParInt*sum{mdl.coeffTermCost[i,1]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1}+
                                  (1-mdl.ParInt)*sum{mdl.coeffTermCost[i,2]*mdl.z_Ol[N+1,6]^(order+1-i),i=1:order+1}))
    
    
    #@NLobjective(mdl, Min, costZTerm + constZTerm + derivCost + controlCost + laneCost)
    #@NLobjective(mdl,Min,sum{(z_Ol[i,1]-0.5)^2,i=1:N+1})

    @NLobjective(mdl.mdl, Min, constZTerm + costZTerm + mdl.derivCost + mdl.controlCost + mdl.laneCost)



    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)
    println("Solved")

    println("derivCost: $(getvalue(mdl.derivCost))")
    println("controlCost: $(getvalue(mdl.controlCost))")
    println("termCost: $(getvalue(mdl.costZTerm))")
    println("mdl.termConst: $(getvalue(mdl.constZTerm))")
    println("termConst: $(getvalue(constZTerm))")
    println("laneCost: $(getvalue(mdl.laneCost))")

    sol_u       = getvalue(mdl.u_Ol)
    sol_z       = getvalue(mdl.z_Ol)
    println("Solution u = $sol_u")
    println("Predicting until z = $(sol_z[end,6])")

    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    mpcSol.cost = zeros(6)
    mpcSol.cost = [0,getvalue(mdl.costZTerm),getvalue(mdl.constZTerm),getvalue(mdl.derivCost),getvalue(mdl.controlCost),getvalue(mdl.laneCost)]
    
    println("--------------- MPC END ------------------------------------------------")
    nothing
end

function solveMpcProblem_pathFollow(mdl::MpcModel_pF,mpcSol::MpcSol,mpcParams::MpcParams,trackCoeff::TrackCoeff,posInfo::PosInfo,modelParams::ModelParams,zCurr::Array{Float64},uCurr::Array{Float64})

    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}

    sol_status::Symbol
    sol_u::Array{Float64,2}
    sol_z::Array{Float64,2}

    # Update current initial condition, curvature and previous input
    setvalue(mdl.z0,zCurr)
    setvalue(mdl.uCurr,uCurr[:])
    setvalue(mdl.coeff,coeffCurvature)

    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)
    sol_u       = getvalue(mdl.u_Ol)
    sol_z       = getvalue(mdl.z_Ol)

    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    mpcSol.cost = zeros(6)
    #mpcSol.cost = [getvalue(mdl.costZ),0,0,getvalue(mdl.derivCost),getvalue(mdl.controlCost),0]

    nothing
end
