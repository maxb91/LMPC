function deleteInfeasibleTrajectories!(oldTraj,posInfo::classes.PosInfo,obstacle, pred_obst::Array{Float64,2}, i::Int64)

    index_first = Array{Int64}(oldTraj.n_oldTraj)
    index_last = Array{Int64}(oldTraj.n_oldTraj)
    if (obstacle.s_obstacle[i,1]-obstacle.rs) - posInfo.s <=2.0 && (obstacle.s_obstacle[i,1]+obstacle.rs)>= posInfo.s #if our car is closer than one meter to obstacle and not fully after it
        for k =1:oldTraj.n_oldTraj
            index_first[k]  = findfirst(x -> x>posInfo.s, oldTraj.oldTraj[:,1,k])
            index_last[k] = findfirst(y -> y>pred_obst[end,1]+obstacle.rs, oldTraj.oldTraj[:,1,k])
            for ii = index_first[k]:index_last[k], kk = 1:mpcParams.N+1
                if ((oldTraj.oldTraj[ii,1,k]-pred_obst[kk,1])/obstacle.rs )^2 + ( (oldTraj.oldTraj[ii,2,k]-pred_obst[kk,2])/obstacle.ry )^2 <= 1
                    setvalue(m.ssInfOn[k],1000)#1500
                    break
                end
            end
        end
    end
end


function addOldtoNewPos(oldTraj)
#take old trajectory with obstacel near car and add traj for better learning
end

# function addOldTraj()
# if dist2obstacleS < 2*v_diff && dist2obstacleS > -length_obstacle/2
#     search in oldTrajectories: states with dist2obstacleS_old  = dist2obstacleS+- eps
#         check if y_obstacle  == y_obstacle_old +- tol for whole prediction horizon
#             for current_pos:curent_pos+max_pred_horizon  check if curvature = curvature_old +- tol_c
#                 copy all safe set states to dummy safe set with shift in  state s
#                 activate dummy traj in mpc solver
# end
# solve mpc problem
# deactivate dummy trajectory


# also safe the curvature at all iterations