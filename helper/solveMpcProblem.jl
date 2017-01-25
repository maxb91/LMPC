# Variable definitions
# m.z_Ol[i,j] = z_OpenLoop, open loop prediction of the state, i = state, j = step

# States:
# i = 1 -> s
# i = 2 -> ey
# i = 3 -> epsi
# i = 4 -> v

function solveLearningMpcProblem!(m::initLearningModel,mpcSol::classes.MpcSol,mpcCoeff::classes.MpcCoeff,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,
    lapStatus::classes.LapStatus,posInfo::classes.PosInfo,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64}, pred_obst, iter::Int64)
    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    

    # @show()


    local sol_u::Array{Float64,2} 
    local sol_z::Array{Float64,2} 
    # Update current initial condition
    setvalue(m.z0,zCurr')


    # # Update model values
    setvalue(m.coeff,coeffCurvature)
    setvalue(m.coeffTermCost,mpcCoeff.coeffCost[iter,:,:])
    setvalue(m.coeffTermConst,mpcCoeff.coeffConst[iter,:,:,:])
    setvalue(m.uCurr,uCurr)
    setvalue(m.sCoord_obst[:,1],pred_obst[:,1])
    setvalue(m.sCoord_obst[:,2],pred_obst[:,2])

    #println("Model formulation:")
    #println(m.m)
    # Solve Problem and return solution
    sol_status  = solve(m.mdl)
    counter = 0
    while sol_status != :Optimal && counter <= 10
        sol_status  = solve(mdl.mdl)
        counter += 1
        println("Not solved optimally, trying again...")
        println("state = ",zCurr)
    end


    sol_u       = getvalue(m.u_Ol)
    sol_z       = getvalue(m.z_Ol)
    mpcSol.lambda = getvalue(m.lambda)
    mpcSol.ssInfOn = getvalue(m.ssInfOn)
    mpcSol.eps[:,iter] = getvalue(m.eps)
    # c_print = getvalue(m.c)
    if iter%20 ==0
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
    mpcSol.cost = [getvalue(m.costZ);getvalue(m.costZTerm);getvalue(m.constZTerm);getvalue(m.derivCost);getvalue(m.controlCost);getvalue(m.laneCost);getvalue(m.velocityCost); getvalue(m.costObstacle)]
    #objvl = getobjectivevalue(m.m)
    nothing #nothing to return
end





function solvePathFollowMpc!(m::initPathFollowingModel,mpcSol::classes.MpcSol,mpcCoeff::classes.MpcCoeff,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,
    lapStatus::classes.LapStatus,posInfo::classes.PosInfo,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64}, pred_obst, iter::Int64)
    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    

    local sol_u::Array{Float64,2} 
    local sol_z::Array{Float64,2} 
    # Update current initial condition
    setvalue(m.z0,zCurr')

    # Update model values
    setvalue(m.coeff,coeffCurvature)
    setvalue(m.uCurr,uCurr)
    setvalue(m.sCoord_obst[:,1],pred_obst[:,1])
    setvalue(m.sCoord_obst[:,2],pred_obst[:,2])

    #println("Model formulation:")
    #println(m.m)
    # Solve Problem and return solution
    sol_status  = solve(m.mdl)
    sol_u       = getvalue(m.u_Ol)
    sol_z       = getvalue(m.z_Ol)
    mpcSol.eps[3,iter] = getvalue(m.eps)
    # c_print = getvalue(m.c)

    if iter%20 ==0
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
    mpcSol.cost = [getvalue(m.costPath);0;0;getvalue(m.derivCost);getvalue(m.controlCost);0;getvalue(m.velocityCost); getvalue(m.costObstacle)]
    #objvl = getobjectivevalue(m.m)

    nothing #nothing to return
end
