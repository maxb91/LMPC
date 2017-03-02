function deleteInfeasibleTrajectories!(m::initLearningModel,oldTraj,distance2obst::Float64,obstacle, close_pred_obst::Array{Float64,2}, i::Int64, zCurr_s::Array{Float64,2},dt::Float64)

    v_ego = zCurr_s[i,4]
    v_obst = close_pred_obst[1,3] # close_pred_obst[1,3] current velocity ofclosest car
    courage_factor = 1.0 # courage factor between 0 and 1 smaller factor excludes infeas traj late


    index_first = Array{Int64}(oldTraj.n_oldTraj)
    index_last = Array{Int64}(oldTraj.n_oldTraj)
    if (distance2obst < 2.0 && distance2obst > - 2*obstacle.rs)#if our car is closer than two meter to obstacle and not fully after it
        for k =1:oldTraj.n_oldTraj #search in all trajectories in the safe set 
            index_first[k]  = findfirst(x -> x > posInfo.s, oldTraj.oldTraj[:,1,k])# find the next greater point in the safe set to the current position for the current regarded old trajectory 
            index_last[k] = findfirst(y -> y > close_pred_obst[end,1]+obstacle.rs, oldTraj.oldTraj[:,1,k])# find the next point in the old trajetory after the predicted obstacle position
            if index_first[k] == 0 || index_last[k] == 0# if this is true probably an error has occured
                warn("no s value in old Traj greater than current position.")
                break
            end
            for ii = index_first[k]:index_last[k], kk = 1:mpcParams.N+1
                v_diff = v_ego - v_obst
                            

                t2collision = (distance2obst/v_diff)/courage_factor
                if v_diff < 0 && distance2obst > 0.0 # if the obstacle is faster than the ego vehicle and the ego vehicle is still behind the obstacle
                    t2collision = 100.0#the cars wont collied soon. in the distance2obst is still a safety marging the time to collision is set to a high  value so this trajectory wont be excluded
                    # println("ego vehicle is slower. s =  $(zCurr_s[i,1])")
                end
                if distance2obst <= 0 # trajectories should always be excluded if we are currently next to the obstacle if both the distance and v_diff are negativ the time2collision might become a big number so we set collision to 0 instead
                    # println("ego vehicle is next to obst. s =  $(zCurr_s[i,1])")
                    t2collision = 0.0
                end
                if ((oldTraj.oldTraj[ii,1,k]-close_pred_obst[kk,1])/obstacle.rs )^2 + ( (oldTraj.oldTraj[ii,2,k]-close_pred_obst[kk,2])/obstacle.ry )^2 <= 1 && t2collision<=1.0
                
                    setvalue(m.ssInfOn[k],1500)#set ssInOn to a high value to exclude the corresponding trajectory from the safe set
                    break
                end
            end
        end
    end
end


function addOldtoNewPos(oldTraj, distance2obst::Float64, obstacle, iter::Int64, close_pred_obst::Array{Float64}, mpcParams::classes.MpcParams, zCurr_s::Array{Float64},dt::Float64, mpcCoeff::classes.MpcCoeff)

    pLength = mpcCoeff.pLength
    v_feas_traj = zeros(oldTraj.n_oldTraj)
    feasible_starting_indeces= zeros(oldTraj.n_oldTraj)
    exist_feas_traj = 0
    s_horizon = zeros(mpcParams.N+1)
    s_horizon[1] = zCurr_s[iter,1] # s position at the first prediction step is the current position
    v_ego = zCurr_s[iter,4]
    N_points= size(oldTraj.oldTraj,1) 

    eps0 = 0.05 #tol distance
    eps1 = 0.3 #tol v_obst
    eps2 = 0.01 #tol e_y
    eps3 = 0.1#tol curvature

    for i = 1:mpcParams.N
        s_horizon[i+1] = s_horizon[i]+v_ego*dt #we assume the ego vehicle will continue to travel with its current speed
    end
    s_diff = s_horizon[end]-s_horizon[1]# calculate how many meters the ego vehicle will probably travel in the prediction horizon



                    #####plot test to see current approx curvature and the curvature in the safe set
                    # curv_points = 30
                    # curvature_curr = zeros(curv_points)
                    # s_curv = linspace(zCurr_s[iter,1],zCurr_s[iter,1]+s_diff,curv_points )
                    # for inp = 1:curv_points
                    #     curvature_curr[inp] = trackCoeff.coeffCurvature[1]*s_curv[inp]^4+trackCoeff.coeffCurvature[2]*s_curv[inp]^3+trackCoeff.coeffCurvature[3]*s_curv[inp]^2+trackCoeff.coeffCurvature[4]*s_curv[inp] + trackCoeff.coeffCurvature[5]
                    # end
                    # plot(s_curv,curvature_curr, color = "blue") 

                    # distance_s =(oldTraj.oldTraj[:,1,1]-zCurr_s[iter,1]).^2
                    # index_s = findmin(distance_s,1)[2]
                    # ind_last_old  = findfirst(x -> x > (oldTraj.oldTraj[index_s,1,1]+s_diff)[1], oldTraj.oldTraj[:,1,1])
                    # plot_old = oldTraj.curvature[index_s[1]:ind_last_old,1]
                    # plot_old_s= linspace(oldTraj.oldTraj[index_s,1,1][1],(oldTraj.oldTraj[index_s,1,1]+s_diff)[1],size(plot_old)[1])
                    # plot(plot_old_s,plot_old, color = "red")
                     #######end test

    if distance2obst < 2.0 && distance2obst > - 2*obstacle.rs
        for k =1:oldTraj.n_oldTraj-2 #do not search last trajectory as it to be replaced. do not search second last as it contains path following round
            for i = 1:oldTraj.oldNIter[k], l = 1:obstacle.n_obstacle #search for all time steps of the rounds and for all obstacles in the track
                if distance2obst <= oldTraj.distance2obst[i,k,l] + eps0 && distance2obst >= oldTraj.distance2obst[i,k,l] - eps0 # set if the distance two closet obstacle is similar
                        compare_v_ey = 0 #this counter will be increased if the obsatcle has sufficient similarities. Only if teh similarities are present for all prediction steps the trajectory of this round can be copied
                        for i_pred =0:mpcParams.N  
                            if close_pred_obst[i_pred+1,3] <= obstacle.v[i+i_pred,k,l] + eps1 && close_pred_obst[i_pred+1,3] >= obstacle.v[i+i_pred,k,l] - eps1 &&
                                close_pred_obst[i_pred+1,2] <= obstacle.sy_obstacle[i+i_pred,k,l] +eps2 && close_pred_obst[i_pred+1,2] >= obstacle.sy_obstacle[i+i_pred,k,l] -eps2
                                compare_v_ey+=1
                            end
                        end
                    last_considered  = findfirst(x -> x > oldTraj.oldTraj[i,1,k]+s_diff, oldTraj.oldTraj[:,1,k]) #find the index forthe old trajectory whose s values is 2 meters greater than the current s value
                    
                    curv_points = last_considered-i+1 #the number of curvature points which have to be similar in the old traj to the current position. 
                    if curv_points<1 # will just be active if something goes wrong. is applied to not lose all data through an error in thath case
                        warn("something went wrong couldnt find an s position $s_diff m infornt of car")
                        break
                    end
                    curvature_curr = zeros(curv_points)
                    s_curv = linspace(zCurr_s[iter,1],zCurr_s[iter,1]+s_diff,curv_points )
                    compare_curv = 0
                    for inp = 1:curv_points
                        #see if the curvature is similar for all points in the old trajectory which could be reached in the prediction horzion
                        curvature_curr[inp] = trackCoeff.coeffCurvature[1]*s_curv[inp]^4+trackCoeff.coeffCurvature[2]*s_curv[inp]^3+trackCoeff.coeffCurvature[3]*s_curv[inp]^2+trackCoeff.coeffCurvature[4]*s_curv[inp] + trackCoeff.coeffCurvature[5]
                        if curvature_curr[inp] <= oldTraj.curvature[i+inp-1,k]+eps3 && curvature_curr[inp] >=  oldTraj.curvature[i+inp-1,k]-eps3 
                            compare_curv+=1
                        end       
                    end

                    if curv_points == compare_curv && compare_v_ey == mpcParams.N+1 # if the obstacle and curvature properties are fullfilled this round could be used for copying
                        feasible_starting_indeces[k]=i # we therefore save the index along the track which oculd be used for copying
                        v_feas_traj[k] = oldTraj.oldTraj[i+mpcParams.N,4,k]# if several traj can be copied we need to have a criteria to determine which traj will be copied, here we use high velocity at the end of the prediction horzion in the safe set and not low cost because trajectories in diffrent point of track have different magnitufe of cost despite the speed of completion
                        exist_feas_traj = 1#if we found at least one trajectroy eligible for copying the next if statement is initiated
                        #############
                        #just for debugging to see how similar curvature profiles are
                        # if k ==1
                        #     clf()
                        #     plot(s_curv,curvature_curr, color = "blue")#, label ="curr curvature")
                        #     plot(s_curv,oldTraj.curvature[i:i+curv_points-1,k], color= "red")#, label ="old curvature")
                        #     readline()
                        # end
                        ################
                        break #break out of for loop to check next trajectory once we found a possible tarj to copy in one round all other states in this round wont be searched 
                    end
                end
            end
        end
    end
    copy_info = zeros(4)
    if exist_feas_traj == 1  #&& 0 == 1 #activate the comment to deactive the copying functionality
        index_of_traj_2_copy = findmax(v_feas_traj)[2] #the fastest traj out of all that are eligible for copying is chosen. If velocities are equal takes first value in array. least old trajectory
        ind_start = convert(Int64,feasible_starting_indeces[index_of_traj_2_copy])
        s_start = oldTraj.oldTraj[ind_start,1,index_of_traj_2_copy] #the s value from which the trajectory onwards will be copied
        oldTraj.oldTraj[:,1,end] = oldTraj.oldTraj[:,1,index_of_traj_2_copy]-s_start+zCurr_s[iter,1] #copy the s coordinates of the fastest trajectory but shift it towards the current s postion of the ego vehicle
        oldTraj.oldTraj[:,2:4,end] = oldTraj.oldTraj[:,2:4,index_of_traj_2_copy] #copy the furter states of the trajectory in order for the solver to be able to rely on them

        # we calculate the points with the closest distance and the corrsponding cost for all old trajectories. This was to be able to adjust the cost not to copied old traje but to the cheapeast old trajectory. now this would not be necessary anymore  as we just need the cost of the copied traj
        distance_s = (oldTraj.oldTraj[:,1,1:end-1]-zCurr_s[iter,1]).^2
        index_s = findmin(distance_s,1)[2]
        costs = zeros(oldTraj.n_oldTraj-1)
        for k = 1:oldTraj.n_oldTraj-1
            costs[k]= oldTraj.cost2Target[index_s[k]-N_points*(k-1),k]
        end 

        #########################################
        #######adjust the cost to curent position
        #the cost on the old trajectory at the first index that gets copied
        costCopiedBeginCopying         = oldTraj.cost2Target[ind_start,index_of_traj_2_copy][1]
        #the cost on the old trajectory at the position closest to the current vehicle position.
        costCopiedCloseCurrentPosition = oldTraj.cost2Target[index_s[index_of_traj_2_copy]-N_points*(index_of_traj_2_copy-1),index_of_traj_2_copy][1] #the term N_points*(index_of_traj_2_copy-1) is used because of the array indexing in julia that comes form the find function
        #we adjust the cost along the trajectory to start at the absolute value that the copied trajectry has at the current ego vehicle position
        oldTraj.cost2Target[ind_start:ind_start+pLength,end] = oldTraj.cost2Target[ind_start:ind_start+pLength,index_of_traj_2_copy]-(costCopiedBeginCopying-costCopiedCloseCurrentPosition)
        #########################################

        copy_info =[index_of_traj_2_copy,s_start, zCurr_s[iter,1],0]# save the data of the copied traj for future plotting
        # println("copied s :$(oldTraj.oldTraj[ind_start,1,index_of_traj_2_copy]) from old round : $index_of_traj_2_copy, curr s: $(zCurr_s[iter,1]), iterartion : $iter")
       

    else
    # if we found no simialr trajectory the last trajectory is filled with the pathfollowing round and (most probably) wont be used in the solution
        oldTraj.oldTraj[:,:,end]  = oldTraj.oldTraj[:,:,end-1]
        oldTraj.cost2Target[:,end] = oldTraj.cost2Target[:,end-1]
    end
    return copy_info
end


