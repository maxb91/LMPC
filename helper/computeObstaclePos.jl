function computeObstaclePos!(obstacle::Obstacle, dt::Float64, i::Int64)
    
    # s_obstacle
    # sy_obstacle
    # rs
    # ry
    v = 0.4# [m/s]
    tol_mid = 0.001 # tolrance to make sure the center of the obstacle should not be on the center of the track because this leads to two solutions with the same optimality
    #keep the form and values of current obstacle and update s and sy values
    

    obstacle.s_obstacle[i+1] = obstacle.s_obstacle[i] + v * dt
    obstacle.sy_obstacle[i+1] = obstacle.sy_obstacle[i]

    # obstacleNext.sy_obstacle = obstacleNext.sy_obstacle + rand(1,1)[1]/2*dt #!! this give only pos rand values
    # if - tol_mid <= obstacleNext.sy_obstacle <= tol_mid
    #     obstacleNext.sy_obstacle  = obstacleNext.sy_obstacle-0.01
    # end

 nothing
end
