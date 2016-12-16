function deleteInfeasibleTrajectories!(oldTraj,posInfo::classes.PosInfo,obstacle, pred_obst::Array{Float64,2}, i::Int64, zCurr_x::Array{Float64,2})

    v_ego = zCurr_x[i,4]
    ds = 0.1# get as input
    dt = 0.1
    safety_factor = 0.95

    index_first = Array{Int64}(oldTraj.n_oldTraj)
    index_last = Array{Int64}(oldTraj.n_oldTraj)
    if (obstacle.s_obstacle[i,1]-obstacle.rs) - posInfo.s <=2.2 && (obstacle.s_obstacle[i,1]+obstacle.rs)>= posInfo.s #if our car is closer than one meter to obstacle and not fully after it
        for k =1:oldTraj.n_oldTraj
            index_first[k]  = findfirst(x -> x > posInfo.s, oldTraj.oldTraj[:,1,k])
            index_last[k] = findfirst(y -> y > pred_obst[end,1]+obstacle.rs, oldTraj.oldTraj[:,1,k])
            for ii = index_first[k]:index_last[k], kk = 1:mpcParams.N+1
                v_diff = v_ego - obstacle.v
                premeter2collision = (oldTraj.oldTraj[ii,1,k] - oldTraj.oldTraj[index_first[k],1,k])
                
                # meter2collision = (premeter2collision - (kk-1)*v_ego*dt)
                # t2collision = (meter2collision/v_diff +(kk-1)*dt)/safety_factor
                meter2collision = (premeter2collision - (kk-1)*obstacle.v*dt)
                t2collision = (meter2collision/v_diff)/safety_factor
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


function addOldtoNewPos(oldTraj, distance2obst::Float64, obstacle, iter::Int64)
#take old trajectory with obstacle near car and add traj for better learning


# if distance2obst < 2.0 && distance2obst > -2*obstacle.rs
#     for k =1:oldTraj.n_oldTraj-1 #do not search last trajectory as it to be replaced.
#         for j= 1:oldTraj.oldNIter[k]
#             if distance2obst  <= oldTraj.distance2obst[j,k] + eps && distance2obst  >= oldTraj.distance2obst[j,k] - eps       
#                 if obstacle.sy_obstacle[iter,1]  <= obstacle.sy_obstacle[l,k] +eps2   && obstacle.sy_obstacle[iter,1]  >= obstacle.sy_obstacle[l,k] -eps2                # for whole prediction horizon
# #             for current_pos:curent_pos+max_pred_horizon  check if curvature = curvature_old +- tol_c
# #                 copy all safe set states to dummy safe set with shift in  state s
# #                 activate dummy traj in mpc solver
#                 end
#             end
#         end
#     end
end
# solve mpc problem
# deactivate dummy trajectory


# also safe the curvature at all iterations
