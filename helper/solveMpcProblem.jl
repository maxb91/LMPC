# Variable definitions
# mdl.z_Ol[i,j] = z_OpenLoop, open loop prediction of the state, i = state, j = step

# States:
# i = 1 -> s
# i = 2 -> ey
# i = 3 -> epsi
# i = 4 -> v

function solveMpcProblem!(mdl::classes.MpcModel,mpcSol::classes.MpcSol,mpcCoeff::classes.MpcCoeff,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,lapStatus::classes.LapStatus,posInfo::classes.PosInfo,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64}, obstacle::classes.Obstacle, iter::Int64)
    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    s_obst     = obstacle.s_obstacle[iter]
    sy_obst    = obstacle.sy_obstacle[iter]

   local sol_u::Array{Float64,2} 
   local sol_z::Array{Float64,2} 

    # Update current initial condition
    setvalue(mdl.z0,zCurr')

    # Update model values
    setvalue(mdl.coeff,coeffCurvature)
    setvalue(mdl.coeffTermCost,mpcCoeff.coeffCost)
    setvalue(mdl.coeffTermConst,mpcCoeff.coeffConst)
    setvalue(mdl.uCurr,uCurr)
    setvalue(mdl.sCoord_obst[1],s_obst)
    setvalue(mdl.sCoord_obst[2],sy_obst)
    #println("Model formulation:")
    #println(mdl.mdl)
    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)
    sol_u       = getvalue(mdl.u_Ol)
    sol_z       = getvalue(mdl.z_Ol)
    mpcSol.lambda = getvalue(mdl.lambda)

    # c_print = getvalue(mdl.c)

    if iter%50 ==0
    # println("curvature: $c_print")
    #println("Predicting until s = $(sol_z[end,1])") 
    end
    
    
    # safe data to solution class
    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    mpcSol.cost = zeros(7)
    mpcSol.cost = [getvalue(costZ);getvalue(costZTerm);getvalue(constZTerm);getvalue(derivCost);getvalue(controlCost);getvalue(laneCost); getvalue(costObstacle)]
    #objvl = getobjectivevalue(mdl.mdl)

    nothing #nothing to return
end
