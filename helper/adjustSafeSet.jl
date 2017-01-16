function deleteInfeasibleTrajectories!(oldTraj,posInfo::classes.PosInfo,obstacle, pred_obst::Array{Float64,2}, i::Int64, zCurr_s::Array{Float64,2},dt::Float64)

    v_ego = zCurr_s[i,4]

    courage_factor = 1.0 # courage factor between 0 and 1 smaller factor excludes infeas traj late
    a_breakmax = 0.0

    index_first = Array{Int64}(oldTraj.n_oldTraj)
    index_last = Array{Int64}(oldTraj.n_oldTraj)
    if (obstacle.s_obstacle[i,1]-obstacle.rs) - posInfo.s <=2.2 && (obstacle.s_obstacle[i,1]+obstacle.rs)>= posInfo.s #if our car is closer than one meter to obstacle and not fully after it
        for k =1:oldTraj.n_oldTraj
            index_first[k]  = findfirst(x -> x > posInfo.s, oldTraj.oldTraj[:,1,k])
            index_last[k] = findfirst(y -> y > pred_obst[end,1]+obstacle.rs, oldTraj.oldTraj[:,1,k])
            for ii = index_first[k]:index_last[k], kk = 1:mpcParams.N+1
                obstacle_v = obstacle.v[i,k]-a_breakmax*dt*kk #robust control assume maximum break acceleration
                if obstacle_v < 0
                    obstacle_v = 0
                end
                v_diff = v_ego - obstacle_v
                premeter2collision = (oldTraj.oldTraj[ii,1,k] - oldTraj.oldTraj[index_first[k],1,k])
                
                # meter2collision = (premeter2collision - (kk-1)*v_ego*dt)
                # t2collision = (meter2collision/v_diff +(kk-1)*dt)/safety_factor
                meter2collision = (premeter2collision - (kk-1)*obstacle.v[i,k]*dt)
                t2collision = (meter2collision/v_diff)/courage_factor
                if ((oldTraj.oldTraj[ii,1,k]-pred_obst[kk,1])/obstacle.rs )^2 + ( (oldTraj.oldTraj[ii,2,k]-pred_obst[kk,2])/obstacle.ry )^2 <= 1 && t2collision<=1.0
                  
                    # @show kk
                    # @show k
                    # @show t2collision
                    # println("++++++++++++++++++++++")
                    # println("                         ")

                    # problem with obstacle and safe set
                    # first approach for distance calculation dicuss
                    # do i need stochastic mpc ? 
                    # kinematic model ok or dynamic?

                    setvalue(m.ssInfOn[k],1500)#1500
                    break
                end
            end
        end
    end
end


function addOldtoNewPos(oldTraj, distance2obst::Float64, obstacle, iter::Int64, pred_obst::Array{Float64}, mpcParams::classes.MpcParams, zCurr_s::Array{Float64},dt::Float64, mpcCoeff::classes.MpcCoeff)
#take old trajectory with obstacle near car and add traj for better learning
pLength = mpcCoeff.pLength

v_feas_traj = zeros(oldTraj.n_oldTraj)
feasible_starting_indeces= zeros(oldTraj.n_oldTraj)
exist_feas_traj = 0
s_horizon = zeros(mpcParams.N+1)
s_horizon[1] = zCurr_s[iter,1]
v_ego = zCurr_s[iter,4]
N_points        = size(oldTraj.oldTraj,1) 

eps0 = 0.1 #tol distance
eps1 = 0.01 #tol v_obst
eps2 = 0.01 #tol e_y
eps3 = 0.01#tol curvature

for i = 1:mpcParams.N
    s_horizon[i+1] = s_horizon[i]+v_ego*dt
end
s_diff = s_horizon[end]-s_horizon[1]
#s_diff = 0.3


            #####plot test
                # curv_points = 30
                # curvature_curr = zeros(curv_points)
                # s_curv = linspace(zCurr_s[iter,1],zCurr_s[iter,1]+s_diff,curv_points )
                # for inp = 1:curv_points
                #     curvature_curr[inp] = trackCoeff.coeffCurvature[1]*s_curv[inp]^4+trackCoeff.coeffCurvature[2]*s_curv[inp]^3+trackCoeff.coeffCurvature[3]*s_curv[inp]^2+trackCoeff.coeffCurvature[4]*s_curv[inp] + trackCoeff.coeffCurvature[5]
                # end
                # plot(s_curv,curvature_curr, color = "blue") 
                # distance_s =(oldTraj.oldTraj[:,1,1]-zCurr_s[iter,1]).^2
                # index_s = findmin(distance_s,1)[2]
                # plot_old = oldTraj.curvature[index_s[1]:index_s[1]+curv_points-1,1]
                # plot(s_curv,plot_old, color = "red")
                 #######end test

if distance2obst < 2.0 && distance2obst > - 2*obstacle.rs
    for k =1:oldTraj.n_oldTraj-1 #do not search last trajectory as it to be replaced.
        for i = 1:oldTraj.oldNIter[k]
            if distance2obst <= oldTraj.distance2obst[i,k] + eps0 && distance2obst >= oldTraj.distance2obst[i,k] - eps0
                compare_v_ey = 0
                for i_pred =0:mpcParams.N  
                    #if obstacle.sy_obstacle[iter,1]  <= obstacle.sy_obstacle[l,k] +eps2   && obstacle.sy_obstacle[iter,1]  >= obstacle.sy_obstacle[l,k] -eps2 #the current 
                    if pred_obst[i_pred+1,3] <= obstacle.v[i+i_pred,k] + eps1 && pred_obst[i_pred+1,3] >= obstacle.v[i+i_pred,k] - eps1 &&
                        pred_obst[i_pred+1,2] <= obstacle.sy_obstacle[i+i_pred,k] +eps2 && pred_obst[i_pred+1,2] >= obstacle.sy_obstacle[i+i_pred,k] -eps2
                        compare_v_ey+=1
                    end
                end
                last_considered  = findfirst(x -> x > oldTraj.oldTraj[i,1,k]+s_diff, oldTraj.oldTraj[:,1,k]) #find the index forthe old trajectory whose s values is 2 meters greater than the current s value,where the distance the obstacle is the same
                curv_points = last_considered-i+1
                curvature_curr = zeros(curv_points)
                s_curv = linspace(zCurr_s[iter,1],zCurr_s[iter,1]+s_diff,curv_points )
                compare_curv = 0
                for inp = 1:curv_points
                    curvature_curr[inp] = trackCoeff.coeffCurvature[1]*s_curv[inp]^4+trackCoeff.coeffCurvature[2]*s_curv[inp]^3+trackCoeff.coeffCurvature[3]*s_curv[inp]^2+trackCoeff.coeffCurvature[4]*s_curv[inp] + trackCoeff.coeffCurvature[5]
                    if curvature_curr[inp] <= oldTraj.curvature[i+inp-1,k]+eps3 && curvature_curr[inp] >=  oldTraj.curvature[i+inp-1,k]-eps3 
                        compare_curv+=1
                    end       
                end

                if curv_points == compare_curv && compare_v_ey == mpcParams.N+1
                    feasible_starting_indeces[k]=i
                    v_feas_traj[k] = oldTraj.oldTraj[i+mpcParams.N,4,k] # use high velocity and not low cost because trajectories in diffrent point of track have different magnitufe of cost despite the speed of completion
                    exist_feas_traj = 1
                    # if k ==1
                    #     clf()
                    #     plot(s_curv,curvature_curr, color = "blue")#, label ="curr curvature")
                    #     plot(s_curv,oldTraj.curvature[i:i+curv_points-1,k], color= "red")#, label ="old curvature")
                    #     readline()
                    # end
                    break #break out of for loop to check next trajectory
                end
            end
        end
    end
end

if exist_feas_traj == 1 && 0 == 1 #!! change 
    index_of_traj_2_copy = findmax(v_feas_traj)[2] #finds the maximum velocity at a later point of the trajectory . If velocities are equal takes first value in array. least old trajectory
    ind_start = convert(Int64,feasible_starting_indeces[index_of_traj_2_copy])
    s_start = oldTraj.oldTraj[ind_start,1,index_of_traj_2_copy]
    oldTraj.oldTraj[:,1,end] = oldTraj.oldTraj[:,1,index_of_traj_2_copy]-s_start+zCurr_s[iter,1] #copy a s postion 
    oldTraj.oldTraj[:,2:4,end] = oldTraj.oldTraj[:,2:4,index_of_traj_2_copy]

    distance_s = (oldTraj.oldTraj[:,1,1:end-1]-zCurr_s[iter,1]).^2
    index_s = findmin(distance_s,1)[2]
    costs = zeros(oldTraj.n_oldTraj-1)
    for k = 1:oldTraj.n_oldTraj-1
        costs[k]= oldTraj.cost2Target[index_s[k]-N_points*(k-1),k]
    end
    curr_min_cost, traj_min = findmin(costs,1)
    #implement dnamic model
    # check curvature for copying probelems in curve
    # see if assigned costz make sense to use trajectories
    # plot oldTRaj aprox
    # @show curr_min_cost
    # @show oldTraj.cost2Target[ind_start,traj_min][1]
     #oldTraj.cost2Target[ind_start:ind_start+pLength,traj_min]
    #@show oldTraj.cost2Target[index_s[1]:index_s[1]+pLength,1]
    # @show ind_start
    # @show oldTraj.cost2Target[ind_start:ind_start+pLength,index_of_traj_2_copy]
    oldTraj.cost2Target[ind_start:ind_start+pLength,end] = oldTraj.cost2Target[ind_start:ind_start+pLength,index_of_traj_2_copy]-(oldTraj.cost2Target[ind_start,index_of_traj_2_copy][1]-oldTraj.cost2Target[index_s[index_of_traj_2_copy]-N_points*(index_of_traj_2_copy-1),index_of_traj_2_copy][1])
    println("copied s :$(oldTraj.oldTraj[ind_start,1,index_of_traj_2_copy]) from old round : $index_of_traj_2_copy, curr s: $(zCurr_s[iter,1]), iterartion : $iter")
    # println("cost copied $(oldTraj.cost2Target[ind_start,end])")
    # l=1
    # println("cost last traj $(oldTraj.cost2Target[index_s[l]-N_points*(l-1),l])")
    # println("$(oldTraj.oldTraj[ind_start:ind_start+pLength+3,1,end])")
    # println("$(zCurr_s[iter,1])")
else

    oldTraj.oldTraj[:,:,end]  = oldTraj.oldTraj[:,:,end-1]
    oldTraj.cost2Target[:,end] = oldTraj.cost2Target[:,end-1]
end
end


