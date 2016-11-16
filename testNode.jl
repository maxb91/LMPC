    using JuMP
    using Ipopt
    using PyPlot
    using JLD

    include("helper/classes.jl")
    include("helper/functions.jl")
    include("helper/coeffConstraintCost.jl")
    include("helper/solveMpcProblem.jl")
    include("helper/simModel.jl")
    include("helper/localizeVehicleCurvAbs.jl")
    include("helper/computeObstaclePos.jl")

    #just loads one specified track with distances betwwen points =1 m and returnx x and y coordinatates
    println("loadMap.......")
    include("helper/loadTestMap.jl")
    x_track, y_track = loadTestMap() #load a track to test controller
    println("loaded")

    # Load Variables and create Model:
    #println("Loading and defining variables...")
    #include("helper/createModel.jl")

    # DEFINE PARAMETERS
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
    mdl                         = classes.MpcModel()

    buffersize                  = 601
    close("all")


    #######################################
    z_Init = zeros(4)
    z_Init[1] = 1.81 # x = 1.81 for s = 32     14 in curve
    z_Init[2] = 2.505 # y = 2.505 for s = 32  12.6
    z_Init[3] = 0.94
    z_Init[4]  = 0.4
    n_rounds = 10 
    InitializeParameters(mpcParams,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,obstacle,buffersize,n_rounds)

    posInfo.s_start             = 0.0 #does not get changed with the current version
    posInfo.s_target            = 5.2 #has to be fitted to track , current test track form ugo has 113.2 meters
     
    ##define obstacle x and xy vlaues not used at the moment 
    #for a clean definition of the x,y points the value of s_obstacle has to be the same as one of the points of the source map. 
    # the end semi axes are approximated over the secant of the points of the track. drawing might not be 100% accurate
    s_obst_init = 80.0 
    sy_obst_init = 0.1
    obstacle.s_obstacle[1,:] = s_obst_init#gets overwritten in loop (moving)
    obstacle.sy_obstacle[1,:]  = sy_obst_init#gets overwritten in loop (moving)
    obstacle.rs = 0.5
    obstacle.ry = 0.2

    
    #####################################

    println("Initialize Model........")
    InitializeModel(mdl,mpcParams,modelParams,trackCoeff,z_Init, obstacle)
    println("Initial solve........")
    solve(mdl.mdl)
    println("Initial solve done!")
    println("*******************************************************")
    println("*******************************************************")
    # Simulation parameters
    dt                          = modelParams.dt
    t                           = collect(0:dt:(buffersize-1)*dt) #40


    
    trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]        # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.nPolyCurvature = 4 # has to be 4 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs
    trackCoeff.nPolyXY = 6  # has to be 6 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs

 

   

     #T
    
    i_final = Array{Int64}(n_rounds)
    for j = 1:n_rounds
        i_final[j]= buffersize
    end
    z_pred_log = zeros(mpcParams.N+1,4,length(t),n_rounds)
    u_pred_log = zeros(mpcParams.N,2,length(t),n_rounds)
    lambda_log = zeros(5,length(t),n_rounds)
    cost        = zeros(7,length(t),n_rounds)

    xStates_log = zeros(length(t),4, n_rounds)
    sStates_log = zeros(length(t),4, n_rounds)
    uAppl_log = zeros(length(t),2, n_rounds) 
    j= 1
    for j=1:n_rounds #10
        

        lapStatus.currentLap = j
        tt          = zeros(length(t),1)
        zCurr_s     = zeros(length(t),4)          # s, ey, epsi, v
        zCurr_x     = zeros(length(t),4)          # x, y, psi, v
        uCurr       = zeros(length(t),2)
        z_final_x = zeros(1,4)::Array{Float64,2}
        u_final = zeros(1,2)::Array{Float64,2}
        #T
      
        obstacle.s_obstacle[1,:] = s_obst_init
        obstacle.sy_obstacle[1,:] = sy_obst_init
        #setup point for vehicle on track in first round. gets overwritten in other rounds
        zCurr_x[1,1] = z_Init[1] # x = 1.81 for s = 32     14 in curve
        zCurr_x[1,2] = z_Init[2] # y = 2.505 for s = 32  12.6
        zCurr_x[1,3] = z_Init[3]
        zCurr_x[1,4] = z_Init[4]  # compare value to v_pathfollowing
        
        if j>1               #setup point for vehicle after first round                   # if we are in the second or higher lap
            zCurr_x[1,:] = z_final_x
            # because the track is not closed we always set up the car at the same place each round
            zCurr_x[1,1] = z_Init[1] # x = 1.81 for s = 32     14 in curve
            zCurr_x[1,2] = z_Init[2] # y = 2.505 for s = 32  12.6
            zCurr_x[1,3] = z_Init[3]
            uCurr[1,:] = u_final    
        end
        

        finished    = false

        i = 1
        #############plot
        # figure(8)
        # clf()
        # ax10= subplot(1,1,1)
        # ax10[:plot](x_track',y_track', linestyle="--", color = "yellow", linewidth = 0.5)#plot the racetrack
        # car_plot = ax10[:plot](zCurr_x[i,1], zCurr_x[i,2], color = "blue")
        # obstacle_plot = ax10[:plot](obstacle.xy_vector[i,1], obstacle.xy_vector[i,2], color = "white")
        # ax10[:plot](boundary_up[1,:], boundary_up[2,:],color="green",linestyle="--")
        # ax10[:plot](boundary_down[1,:], boundary_down[2,:],color="green",linestyle="--")
        # grid()
        ##################
        while i<=length(t) && !finished # as long as we have not reached the maximal iteration time for one round or ended the round
        
            # to make it work s start has to grow over time actual it is just always at 0
       
            # the argument i in localizeVehicleCurvAbs  is solely used for debugging purposes plots not needed for control
            # localize takes some time see ProfileView.view()
            zCurr_s[i,:], trackCoeff.coeffCurvature = localizeVehicleCurvAbs(zCurr_x[i,:],x_track,y_track,trackCoeff, i)
            #if the car has crossed the finish line
            if zCurr_s[i,1] >= posInfo.s_target
                println("Reached finish line at step $i")
                finished = true
                #we count up here as the first round ends just as we would do if the loop get terminiated because i is >= length(t). we count it down later on again to get right index
                i = i + 1
                lapStatus.currentIt = i
                break
            end
            #!! println("s = $(zCurr_s[i,1])")
     
            posInfo.s   = zCurr_s[i,1]
            if j > 1
                       tic()
                coeffConstraintCost!(oldTraj,mpcCoeff,posInfo,mpcParams)
                tt1 = toq()
                #println("coeffConstr: $tt1 s")
            end


            
            tic()
	      
            #solve with zCurr_s containing s ey values 
            solveMpcProblem!(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:]',[mpcSol.a_x;mpcSol.d_f], obstacle,i)
        
                
            tt[i]       = toq()
            uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f]
            cost[:,i,j]   = mpcSol.cost

            lambda_log[:,i,j] = mpcSol.lambda
            z_pred_log[:,:,i,j] = mpcSol.z
            u_pred_log[:,:,i,j] = mpcSol.u
          
            #have Zcurr as states XY and simulate from there return XY values of states 
            zCurr_x[i+1,:]  = simModel_x(zCurr_x[i,:],uCurr[i,:],modelParams.dt,modelParams) #!! @show

            #update Position of the Obstacle car        
            computeObstaclePos!(obstacle, dt, i, j)#this funciton computes values for row i+1
            
            
            if i%50 == 0 
                println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")
            end

           ###T
            # N_plot = 20
            # if i%N_plot == 0 
               
            #     car_plot[1][:remove]()
            #     obstacle_plot[1][:remove]()
            #     ax10[:plot](zCurr_x[i-N_plot+1:i,1], zCurr_x[i-N_plot+1:i,2], color = "blue")
            #     car_plot = ax10[:plot](zCurr_x[i,1], zCurr_x[i,2], color = "blue", marker="o")
            #     ax10[:plot](obstacle.xy_vector[i-N_plot+1:i,1], obstacle.xy_vector[i-N_plot+1:i,2], color = "red")
            #     obstacle_plot = ax10[:plot](obstacle.xy_vector[i,1], obstacle.xy_vector[i,2], color = "red", marker="o")          
            # end 
            ###endT

            i = i + 1
            lapStatus.currentIt = i
        

        end
        if i >= length(t)
            println("took too long used whole length t")
        end
        
        i = i-1 
        lapStatus.currentIt -= 1 # has finished lap already so we need to count both coutners one down as we counted them up after we have crossed the finish line
        z_final_x = zCurr_x[i,:]
        u_final = uCurr[i,:]
        println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i])) s")
        println("Finished Lap Nr. $j")

        # Save states in oldTraj:
        # --------------------------------

        saveOldTraj(oldTraj,zCurr_s, zCurr_x,uCurr,lapStatus,buffersize,modelParams.dt)
        xStates_log[:,:,j] = zCurr_x
        uAppl_log[:,:,j] = uCurr
        sStates_log[:,:,j]= zCurr_s


        
        i_final[j] = i
        if j >1 && i_final[j-1] <= i_final[j]
            println("round was not faster. no learning")
            println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        end

        println("*************************************************************************")
        println("Press c to cancel MPC")
        println("Press Enter for next round solving")
        
        a = ' '
        a = readline()
        if a == "c\r\n" 
                break
        end
    end
    n_rounds = j #update n_rounds to represent actual number of simualated rounds
########save date
          println("Save data to file.......")
    filename = string("../LMPCdata/"string(Dates.today()),"-Data.jld")
    jldopen(filename, "w") do file
        addrequire(file, classes) #ensures that custom data types are working when loaded
        
        JLD.write(file, "sStates_log", sStates_log)
        JLD.write(file, "xStates_log", xStates_log)
        JLD.write(file, "uAppl_log", uAppl_log)
        JLD.write(file, "z_pred_log", z_pred_log)
        JLD.write(file, "u_pred_log", u_pred_log)
        JLD.write(file, "lambda_log", lambda_log)
        JLD.write(file, "cost", cost)
        JLD.write(file, "i_final", i_final)
        JLD.write(file, "x_track", x_track)
        JLD.write(file, "y_track", y_track)
        JLD.write(file, "n_rounds", n_rounds)
        JLD.write(file, "trackCoeff", trackCoeff)
        JLD.write(file, "obstacle", obstacle)
        JLD.write(file, "modelParams.dt", modelParams.dt)
        JLD.write(file, "mpcParams", mpcParams)
        JLD.write(file, "buffersize", buffersize)
    end
    println("finished")