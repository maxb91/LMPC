function computeObstaclePos!(obstacle, dt::Float64, i::Int64, x_track::Array{Float64,2}, trackCoeff::classes.TrackCoeff, sCurr::Float64)
    
    # s_obstacle
    # sy_obstacle
    # rs
    # ry
    tol_mid = 0.001 # tolrance to make sure the center of the obstacle should not be on the center of the track because this leads to two solutions with the same optimality
    #keep the form and values of current obstacle and update s and sy values
    s_length_track =size(x_track)[2]*trackCoeff.ds

    obstacle.s_obstacle[i+1,1] = obstacle.s_obstacle[i,1] + obstacle.v[i,1] * dt
    #!! obstacle.s_obstacle[i+1,1] = obstacle.s_obstacle[i+1,1]%s_length_track# if the obstale makes morfe than one round it start at the begining again
    
    obstacle.sy_obstacle[i+1,1] = obstacle.sy_obstacle[i,1]
    obstacle.v[i+1,1] = obstacle.v[i,1]


    # if obstacle.s_obstacle[i,1] >6.5
    #      obstacle.v[i+1,1] = obstacle.v[i,1]-0.2
    #     if obstacle.v[i+1,1] <0
    #         obstacle.v[i+1,1] =0
    #     end
    # end


    # if sCurr >8.55
    #     obstacle.v[i+1,1] = obstacle.v[i,1]+1.3
    #     if obstacle.v[i+1,1] >1.3
    #         obstacle.v[i+1,1] =1.3
    #     end
    # end


    # obstacleNext.sy_obstacle = obstacleNext.sy_obstacle + rand(1,1)[1]/2*dt #!! this give only pos rand values
    # if - tol_mid <= obstacleNext.sy_obstacle <= tol_mid
    #     obstacleNext.sy_obstacle  = obstacleNext.sy_obstacle-0.01
    # end

 nothing
end

function predictObstaclePos(obstacle, modelParams::classes.ModelParams,mpcParams::classes.MpcParams, iter::Int64)
    s_sy_obst = zeros(mpcParams.N+1,3)
    s_sy_obst[1,1] = obstacle.s_obstacle[iter,1]
    s_sy_obst[1,2] = obstacle.sy_obstacle[iter,1]
    s_sy_obst[1,3] = obstacle.v[iter,1]
    for i = 1:mpcParams.N
        s_sy_obst[i+1,1] = s_sy_obst[i,1]+obstacle.v[iter,1]*modelParams.dt
        s_sy_obst[i+1,2] = s_sy_obst[i,2]
        s_sy_obst[i+1,3] = s_sy_obst[i,3]
    end
    return s_sy_obst
end

