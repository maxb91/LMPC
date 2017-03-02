using JuMP
using Ipopt
using JLD

include("helper/classes.jl")
include("helper/functions.jl")
include("helper/initializeModel.jl")
include("helper/coeffConstraintCost.jl")
include("helper/solveMpcProblem.jl")
include("helper/simModel.jl")
include("helper/localizeVehicleCurvAbs.jl")
include("helper/computeObstaclePos.jl")
include("helper/adjustSafeSet.jl")

#just loads one specified track with distances betwwen points =1 m and returnx x and y coordinatates
println("loadMap.......")
include("helper/loadTestMap.jl")
x_track, y_track = loadTestMap() #load a track to test controller
println("loaded")


# Define and initialize variables
println("define types ........")
oldTraj                     = classes.OldTrajectory()
posInfo                     = classes.PosInfo()
mpcCoeff                    = classes.MpcCoeff()
lapStatus                   = classes.LapStatus(1,1)
mpcSol                      = classes.MpcSol()
obstacle                    = classes.Obstacle()
trackCoeff                  = classes.TrackCoeff()      # info about track (at current position, approximated)
modelParams                 = classes.ModelParams()
mpcParams                   = classes.MpcParams()

buffersize                  = 1101 # buffersize has to be big enough to fit all iterations for the path following trajectory
close("all")


#######################################
# specify values to setup vehicle in first round in the xy space
z_Init = zeros(6)
z_Init[1] = 0 #x
z_Init[2] = 0 #y
z_Init[3] = 0.4*cos(z_Init[5]) #v_x
z_Init[4]  = 0.4*sin(z_Init[5]) #v_y
z_Init[5] = 0 #psi
z_Init[6] = 0.0 #psi_dot






load_safeset = true # true if we want to simulate on the basis of a safe set that was created before.
safeset = "data/2017-02-10-16-00-Data.jld" #specify filename which will be loaded if load_safeset == true
n_rounds = 1 # how many rounds we want to simulate in the current executing of the script
active_obstacle = true #if true the obstacle behavior can be specified below, if false the obsatcles will be placed behind the finish line of the track so that they dont affect the solving
continue_obstacle = true #if true the obstacles get initilaized with the values from the loaded safe set, if false that obsatcles get initilaed as stated below
 

obstacle.n_obstacle = 8
s_obst_init=[3, 9, 15, 21, 27, 32, 39, 47]#[3, 9, 15, 21, 25, 32, 39, 45]#[2.5, 10, 17.5, 25, 32.5, 40, 47.5, 55]#[3, 10, 17, 24, 31, 38, 48]#[4, 11, 18, 25, 32, 39, 49] # [4, 13, 22, 31, 40, 49, 58] #
sy_obst_init = -0.22*ones(obstacle.n_obstacle)
v_obst_init = 1.7*ones(obstacle.n_obstacle)



######### 
InitializeParameters(mpcParams,trackCoeff,modelParams,oldTraj,mpcCoeff,mpcSol,lapStatus,obstacle,buffersize) #initillize value for the classes specified above
#########

@show posInfo.s_target            = (size(x_track)[2]-1)*trackCoeff.ds#calculates the length of a given track or can be set to a fixd value if the finish line should not be at the end of the track

#####################################
println("Initialize Model........")
#the controller relies on two independent models for the pathfollowing round and all following learning mpc rounds
mdl_Path = initPathFollowingModel(mpcParams,modelParams,trackCoeff,mpcCoeff, obstacle,oldTraj.n_oldTraj)
mdl_LMPC = initLearningModel(mpcParams,modelParams,trackCoeff,mpcCoeff, obstacle,oldTraj.n_oldTraj)

println("Initial solve........")
solve(mdl_LMPC.mdl)# intial solve is necessary for LMPC model to prevent an invalid number error. seems to be caused by problems with the coefficients of the terminal cost and terminal set in the controller.
println("Initial solve done!")
println("*******************************************************")
# Simulation parameters
dt                          = modelParams.dt
t                           = collect(0:dt:(buffersize-1)*dt)

Pcurvature = zeros(length(t),2) # just for debugging process. Collects the approximated curvature at all positions from the localizing function. 


trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]        # polynomial coefficients for curvature approximation (zeros for straight line)
trackCoeff.nPolyCurvature = 4 # has to be 4 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbs
trackCoeff.nPolyXY = 6  # has to be 6 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbs

z_pred_log = zeros(mpcParams.N+1,4,length(t),n_rounds) #array to log the predicted states for all time steps for later plotting
u_pred_log = zeros(mpcParams.N,2,length(t),n_rounds) #array to log the predicted inputs for all time steps for later plotting
lambda_log = zeros(oldTraj.n_oldTraj,length(t),n_rounds) #array to log the predicted lambdas
ssInfOn_log = zeros(oldTraj.n_oldTraj, length(t), n_rounds) #array to log the exclusion of trajectories. "InfinityCost" for those trajectories
cost        = zeros(8,length(t),n_rounds) #logs the MPC costs with values for all expresions in the cost function (terminal cost, lane cost, input cost etc.)

pred_obst = zeros(mpcParams.N+1,3,obstacle.n_obstacle) # array contains the states of all obstacles for the whole prediction horizon

if load_safeset == true
    SafeSetData = load(safeset)
    oldTraj = SafeSetData["oldTraj"]
    obstacle = SafeSetData["obstacle"]
end

j = 1 #set j to one as we start at the first round j
z_final_x = zeros(1,4)::Array{Float64,2} # the final states 
u_final = zeros(1,2)# and final contorl input at the end every round that are used for initialiezing the car in the new round





for j=1:n_rounds #loop over all rounds
    
    no_solution_found = 0 #counts number of unsuccesful attempts if too many unsuccesful attempts the round will be terminated

    lapStatus.currentLap = j
    tt          = zeros(length(t),1)     #logs the time used by the solveLearningMpcProblem
    tt1         = zeros(length(t),1)    #logs the time used by safe set operations ( addTrajectory, deleteTrajectory, interpolateTrajectory)
    tt_total    = zeros(length(t),1)    #logs the total compuational time for one iteration of the inner loop
    zCurr_s     = zeros(length(t),4)          # s, ey, epsi, v
    zCurr_x     = zeros(length(t),6)          # x, y, v_x,v_y, psi, psi_dot
    uCurr       = zeros(length(t),2)
    distance2obst = 1000*ones(length(t), obstacle.n_obstacle) #distance to the closest obstacle to ego vehicle at all time steps
    curvature_curr = zeros(length(t))
    copyInfo    = zeros(length(t),4) #contains the index of the trajectory that was copied, the s position which had similarities and which the copying was initialized from, the current s position along the track, and the lambda value of the copied trajectory generated by the solver
    

    #setup point for vehicle on track in first round. gets overwritten in other rounds
    zCurr_x[1,1] = z_Init[1] 
    zCurr_x[1,2] = z_Init[2] 
    zCurr_x[1,3] = z_Init[3]
    zCurr_x[1,4] = z_Init[4]
    zCurr_x[1,5] = z_Init[5]
    zCurr_x[1,6] = z_Init[6]

    #increase the index of all obsatcle tarjectories in th safe set as we are about to create new values for the oldTraj with index 1
    for k = oldTraj.n_oldTraj-1:-1:1
        obstacle.s_obstacle[:,k+1,:]  = obstacle.s_obstacle[:,k,:]
        obstacle.sy_obstacle[:,k+1,:] = obstacle.sy_obstacle[:,k,:]
        obstacle.v[:,k+1,:] = obstacle.v[:,k,:]
    end
    
    
    
    if j == 1 && load_safeset == true # if we load a safe set setup the ego vehicle at the position where it last was after it corssed the finish line
        zCurr_x[1,1] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],1,1]
        zCurr_x[1,2] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],2,1]
        zCurr_x[1,3] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],3,1]
        zCurr_x[1,4] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],4,1]
        zCurr_x[1,5] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],5,1]
        zCurr_x[1,6] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],6,1]
    end
    if j>1             #setup point for vehicle after first round                   # if we are in the second or higher lap
        zCurr_x[1,:] = z_final_x
        uCurr[1,:] = u_final
    end
    if (continue_obstacle == true && j == 1) || j>1 # use values in safe set or 
        obstacle.s_obstacle[1,1,:]  = obstacle.s_obstacle[oldTraj.oldNIter[1],2,:]
        obstacle.sy_obstacle[1,1,:] = obstacle.sy_obstacle[oldTraj.oldNIter[1],2,:]
        obstacle.v[1,1,:]           = obstacle.v[oldTraj.oldNIter[1],2,:]
    elseif continue_obstacle == false && j == 1 #use values which were specified above
        obstacle.s_obstacle[1,1,:] = s_obst_init
        obstacle.sy_obstacle[1,1,:] = sy_obst_init
        obstacle.v[1,1,:] = v_obst_init
    end

    ###########iterations learning
    finished    = false #is set to true if finish line is reaches and round is finished
    i = 1 #index for iterations in one round
    while i<=length(t)-1 && !finished # as long as we have not reached the maximal iteration time for one round or ended the round

   
        # the input i in localizeVehicleCurvAbs  is solely used for debugging purposes plots not needed for control
        # localize takes some time see ProfileView.view()
        tic()
        #approximates the states in the s -e?y frame for the controller and approximates a curvature around the current postion and returns corresponding polynomial coefficients
        zCurr_s[i,:], trackCoeff.coeffCurvature = localizeVehicleCurvAbs(zCurr_x[i,:],x_track,y_track,trackCoeff, i, mpcParams.N, modelParams.dt, Pcurvature)
        if i == 1 && zCurr_s[1,1] > 2 # if the xy coordinates from the loaded states are closest to a point at the end of the track not at the beginning of the track in the first step the s coordinate is reset to zero
            # warn("closest point was before finish line, forced s =0")
            # println("x:$(zCurr_x[i,1]), y:$(zCurr_x[i,2])")
            zCurr_s[1,1]= 0
        end
        posInfo.s   = zCurr_s[i,1]
        curvature_curr[i] = trackCoeff.coeffCurvature[1]*posInfo.s^4+trackCoeff.coeffCurvature[2]*posInfo.s^3+trackCoeff.coeffCurvature[3]*posInfo.s^2+trackCoeff.coeffCurvature[4]*posInfo.s +trackCoeff.coeffCurvature[5]# calculate teh curvature just for debugging
        

        # if the ego vehicle is close to the finish line so some of the prediction steps might lie behind it, the s coordinates of obstacles just behind the finish line with s just over zero will be augmented wwith the value of the tarck length in order for the solver to see the proximity of those obsatcles
        if posInfo.s > posInfo.s_target-(modelParams.v_max*(mpcParams.N+1)*dt)
            for l=1:obstacle.n_obstacle
                if obstacle.s_obstacle[i,1,l]< (modelParams.v_max*(mpcParams.N+1)*dt)
                    obstacle.s_obstacle[i,1,l]+=posInfo.s_target
                end
            end
        end
        
        

        #t_absci =toq()

        if zCurr_s[i,1] >= posInfo.s_target #if the car has crossed the finish line
            println("Reached finish line at step $i")
            finished = true
            #we count up here as the first round ends just as we would do if the loop get terminiated because i is >= length(t). we count it down later on again to get right index
            i = i + 1
            lapStatus.currentIt = i
            break
        elseif i >1 && zCurr_s[i-1,1] >= posInfo.s_target-1 && zCurr_s[i,1]<1 # on a closed track the finish line lies close to the starting point. it can be that the discrete calculated s postion is bigger than the end of the track and already in the new round. to still detect the finish line we added this statement
            println("alternative detection of finish line")
            println("Reached finish line at step $i")
            finished = true
            #we count up here as the first round ends just as we would do if the loop get terminiated because i is >= length(t). we count it down later on again to get right index
            i = i + 1
            lapStatus.currentIt = i
            break
        end            

        ##################
        pred_obst = predictObstaclePos(obstacle, modelParams, mpcParams, i) #predict the state of the obsatcel for the whole prediction horizon of the solver
        distance2obst[i,:] = (pred_obst[1,1,:]-obstacle.rs) - posInfo.s #calculate the distacne to all obstacles 
        ind_closest_obst = findmin(abs(distance2obst[i,:]))[2] #find the closest obstacle which is the only one regarded in the solver
        ##################
        tic()
        if j > 1 || load_safeset == true 3 # if the learning MPC controller is active we add traj, deleet, trajectoreis, and interpolate them for the solver
            copyInfo[i,:] = addOldtoNewPos(oldTraj, distance2obst[i,ind_closest_obst],obstacle,i,pred_obst[:,:,ind_closest_obst], mpcParams,zCurr_s,modelParams.dt,mpcCoeff)
            deleteInfeasibleTrajectories!(mdl_LMPC, oldTraj,distance2obst[i,ind_closest_obst],obstacle, pred_obst[:,:,ind_closest_obst], i, zCurr_s,modelParams.dt) # sets the ssInfOn parameter
	        
            coeffConstraintCost!(oldTraj,mpcCoeff,posInfo,mpcParams,i)
        end
        tt1[i] = toq()
        tic()


        #####################
        if j == 1 && load_safeset == false #use pathfollowing controller in the first round with a pre-loaded safe set
            solvePathFollowMpc!(mdl_Path,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:]',[mpcSol.a_x;mpcSol.d_f], pred_obst[:,:,ind_closest_obst],i)
        elseif j >= 2 || load_safeset == true #use lmpc for all other rounds
            solstat = solveLearningMpcProblem!(mdl_LMPC,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:]',[mpcSol.a_x;mpcSol.d_f], pred_obst[:,:,ind_closest_obst],i)
            if solstat == false
                no_solution_found += 1
            else
                no_solution_found = 0
            end
            if no_solution_found >= 3
                warn("Over 3 unsuccesful iterations. Abort solving!")
                break
            end    
            setvalue(mdl_LMPC.ssInfOn,ones(oldTraj.n_oldTraj))# reset the all trajectories to on position with ssInfOn parameter
        end
        tt[i]= toq()
        

        uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f] 
        zCurr_x[i+1,:]  = simModel_exact_dyn_x(zCurr_x[i,:],uCurr[i,:],modelParams.dt,modelParams)#simulate the model with the inputs generated by the controller

        #update Position of the Obstacle car for the next iteration in this round     
        computeObstaclePos!(obstacle, dt, i, x_track, trackCoeff, zCurr_s[i,1], active_obstacle) #this funciton computes values for row i+1
        
        #save the solution of the controller for future plotting
        cost[:,i,j]         = mpcSol.cost
        lambda_log[:,i,j]   = mpcSol.lambda
        z_pred_log[:,:,i,j] = mpcSol.z
        u_pred_log[:,:,i,j] = mpcSol.u
        ssInfOn_log[:,i,j]  = mpcSol.ssInfOn
        copyInfo[i,4] = mpcSol.lambda[end]#the last traj for the solver is the copied trajectory.
        tt_total[i]= toq()


        #just for debugging
        # if mpcSol.lambda[end] > 0.6
        #     # println("copied s :$(copyInfo[i,2]) from old round : $(copyInfo[i,1]), curr s: $(copyInfo[i,3]), lambda : $(mpcSol.lambda[end])")
        # end

        #print every 50th iteration to see that script is still running
        if i%50 == 0 
            println(" Time: $(tt[i]) s, Solving step $i of $(length(t)), s = $(zCurr_s[i,1]) - Status: $(mpcSol.solverStatus)")
            #println(" Time for whole iteration: $(tt2) s")
            # if j > 1 || load_safeset == true
            #     println("ssInfOn time: $tt1 s")
            # end
        end
        # if tt[i] >0.08 #if solving takes long
        #     println(" Time: $(tt[i]) s, Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus)")
        #     println(" Time: $(tt_total[i]) s for the whole step of the loop ")
        #     println("                           ")
        # end
        # @show mpcSol.solverStatus
        if string(mpcSol.solverStatus) != "Optimal" #print every time the solver fails to find an optimal solution
            println(" Time: $(tt[i]) s, Solving step $i of $(length(t)), s = $(zCurr_s[i,1]) - Status: $(mpcSol.solverStatus)")
        end
        

        i = i + 1#count up i for next iteration in this round
        lapStatus.currentIt = i
    

    end
    #at this point either the finish line is crossed or the array to safe data is full. 
    if i >= length(t)
        println("used whole length t. either adapt array size to track or problems in solving process")
    end
    
    i = i-1 
    lapStatus.currentIt -= 1 # has finished lap already so we need to count both coutners one down as we counted them up after we have crossed the finish line
    z_final_x = zCurr_x[i,:]
    u_final = uCurr[i,:]
    println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i])) s")
    println("Finished Lap Nr. $j")

    # Save states in oldTraj:
    # --------------------------------
    tic()

    #save the obsatcle and ego vehicle data into for the safe set and for future plotting
    saveOldTraj(oldTraj,zCurr_s, zCurr_x,uCurr,lapStatus,buffersize,modelParams.dt, load_safeset, cost[:,:,j],
      lambda_log[:,:,j],z_pred_log[:,:,:,j],u_pred_log[:,:,:,j],ssInfOn_log[:,:,j], mpcSol, obstacle,distance2obst, curvature_curr, copyInfo)
    tt3= toq()



    println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print_with_color(:grey,"Max time to delete/add traj:")
    print_del_time = formatFloat!(maximum(tt1),3)
    print_with_color(:red," $print_del_time\n")

    print_with_color(:grey,"Max time to solve:")
    print_solve_time= formatFloat!(maximum(tt),3)
    print_with_color(:red," $print_solve_time\n")

    print_with_color(:grey,"Max time for calc of one loop")
    print_max_time = formatFloat!(maximum(tt_total),3)
    print_with_color(:red," $(print_max_time)\n")
    
    print_with_color(:grey,"Number of copied Trajectories: ") 
    print_with_color(:yellow,"$(length(find(f->f!=0,copyInfo[:,1])))\n")

    print_with_color(:grey,"Number of used copied Trajectories: ") 
    print_with_color(:yellow,"$(length(find(f->f!=0,copyInfo[:,4])))\n")

    println(" Time to save and overwrite trajectories: $(tt3) s")
    ###############
    println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    oldTraj.oldNIter[1] = i#safe how many iterations it took to reach the finish line (or until the array was full)
    if j>1 && oldTraj.cost2Target[1,2] <= oldTraj.cost2Target[1,1]
        warn("round was not faster. cost : $(oldTraj.cost2Target[1,1])")
        println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    else
        println("cost of round: $(oldTraj.cost2Target[1,1])")
    end

    println("*************************************************************************")


    # figure()
    # plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],oldTraj.curvature[1:oldTraj.oldNIter[j],j], color = "red")
    # plot(Pcurvature[1:oldTraj.oldNIter[j],1], Pcurvature[1:oldTraj.oldNIter[j],2], color = "green")
end






#set filename to safe data
#### numbering to generate results to keep 
filename = string("data/"string(Dates.today()),"-",Dates.format(now(), "HH-MM"),"-Data.jld")
if isfile(filename)
    filename = string("data/"string(Dates.today()),"-",Dates.format(now(), "HH-MM"),"-Data-2.jld")
    warn("File already exists. Added extension \"-2\" ")
end
println("Save data to $filename .......")


#safe data to file
jldopen(filename, "w") do file
    #addrequire(file, classes) #ensures that custom data types are working when loaded

    JLD.write(file, "x_track", x_track)
    JLD.write(file, "y_track", y_track)
    JLD.write(file, "trackCoeff", trackCoeff)
    JLD.write(file, "obstacle", obstacle)
    JLD.write(file, "modelParams", modelParams)
    JLD.write(file, "mpcParams", mpcParams)
    JLD.write(file, "buffersize", buffersize)
    JLD.write(file, "oldTraj", oldTraj)
    JLD.write(file, "mpcCoeff",mpcCoeff)

end
println("finished")
