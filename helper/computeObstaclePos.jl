function computeObstaclePos!(obstacle::classes.Obstacle, dt::Float64, i::Int64, j::Int64, x_track::Array{Float64,2}, trackCoeff::classes.TrackCoeff)
    
    # s_obstacle
    # sy_obstacle
    # rs
    # ry
    v = 0.0# [m/s]
    tol_mid = 0.001 # tolrance to make sure the center of the obstacle should not be on the center of the track because this leads to two solutions with the same optimality
    #keep the form and values of current obstacle and update s and sy values
    s_length_track =size(x_track)[2]*trackCoeff.ds

    obstacle.s_obstacle[i+1,j] = obstacle.s_obstacle[i,j] + v * dt
    obstacle.s_obstacle[i+1,j] = obstacle.s_obstacle[i+1,j]%s_length_track# if the obstale makes morfe than one round it start at the begining again
    
    obstacle.sy_obstacle[i+1,j] = obstacle.sy_obstacle[i,j]

    # obstacleNext.sy_obstacle = obstacleNext.sy_obstacle + rand(1,1)[1]/2*dt #!! this give only pos rand values
    # if - tol_mid <= obstacleNext.sy_obstacle <= tol_mid
    #     obstacleNext.sy_obstacle  = obstacleNext.sy_obstacle-0.01
    # end

 nothing
end
