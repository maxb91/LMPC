# Variable definitions
# m.z_Ol[i,j] = z_OpenLoop, open loop prediction of the state, i = state, j = step

# States:
# i = 1 -> s
# i = 2 -> ey
# i = 3 -> epsi
# i = 4 -> v

function solveMpcProblem!(m::classes.MpcModel,mpcSol::classes.MpcSol,mpcCoeff::classes.MpcCoeff,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,lapStatus::classes.LapStatus,posInfo::classes.PosInfo,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64}, obstacle::classes.Obstacle, iter::Int64)
    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    s_obst     = obstacle.s_obstacle[iter]
    sy_obst    = obstacle.sy_obstacle[iter]

    local sol_u::Array{Float64,2} 
    local sol_z::Array{Float64,2} 
  
    

    # Update current initial condition
    setvalue(m.z0,zCurr')

    # Update model values
    setvalue(m.coeff,coeffCurvature)
    setvalue(m.coeffTermCost,mpcCoeff.coeffCost)
    setvalue(m.coeffTermConst,mpcCoeff.coeffConst)
    setvalue(m.uCurr,uCurr)
    setvalue(m.sCoord_obst[1],s_obst)
    setvalue(m.sCoord_obst[2],sy_obst)


    #println("Model formulation:")
    #println(m.m)
    # Solve Problem and return solution
    sol_status  = solve(m.mdl)
    sol_u       = getvalue(m.u_Ol)
    sol_z       = getvalue(m.z_Ol)
    mpcSol.lambda = getvalue(m.lambda)
    mpcSol.ssOn = getvalue(m.ssOn)
    # c_print = getvalue(m.c)

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
    mpcSol.cost = [getvalue(m.costZ);getvalue(m.costZTerm);getvalue(m.constZTerm);getvalue(m.derivCost);getvalue(m.controlCost);getvalue(m.laneCost); getvalue(m.costObstacle)]
    #objvl = getobjectivevalue(m.m)

    nothing #nothing to return
end
