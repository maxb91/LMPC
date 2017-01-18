    using JuMP
    using Ipopt
    using PyPlot
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
    m                           = classes.MpcModel()

    buffersize                  = 501 #1501
    close("all")


    #######################################
    z_Init = zeros(6)
    z_Init[1] = 0 # x = 1.81 for s = 32     14 in curve
    z_Init[2] = 0 # y = 2.505 for s = 32  12.6
    z_Init[5] = 0 #psi
    z_Init[3] = 0.4*cos(z_Init[5]) #v_x
    z_Init[4]  = 0.4*sin(z_Init[5]) #v_y
    z_Init[6] = 0.0
   
    load_safeset = true#currently the safe set has to contain the same number of trajectories as the oldTraj class we initialize
    safeset = "data/2017-01-17-23-29-Data.jld"

    #########
    InitializeParameters(mpcParams,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,obstacle,buffersize)
    mpcSol.u  = zeros(mpcParams.N,2)
    mpcSol.z  = zeros(mpcParams.N+1,4) 
    mpcSol.lambda = zeros(oldTraj.n_oldTraj)
    mpcSol.lambda[1] = 1
    mpcSol.eps = zeros(2,buffersize)
    #########



    posInfo.s_start             = 0.0 #does not get changed with the current version
    posInfo.s_target            = (size(x_track)[2]-1)*trackCoeff.ds#59.5 #has to be fitted to track , current test track form ugo has 113.2 meters
     
    ##define obstacle x and xy vlaues not used at the moment 
    #for a clean definition of the x,y points the value of s_obstacle has to be the same as one of the points of the source map. 
    # the end semi axes are approximated over the secant of the points of the track. drawing might not be 100% accurate
    s_obst_init = 9.56 
    sy_obst_init = -0.2
    v_obst_init = 0#1.8#1.5#1.5##1.8
    obstacle.rs = 0.5 # if we load old trajecory these values get overwritten
    obstacle.ry = 0.19 # if we load old trajecory these values get overwritten
    
    


    
    #####################################
    println("Initialize Model........")
    InitializeModel(m,mpcParams,modelParams,trackCoeff,z_Init, obstacle,oldTraj)
    println("Initial solve........")
    solve(m.mdl)
    println("Initial solve done!")
    println("*******************************************************")
    println("*******************************************************")
    # Simulation parameters
    dt                          = modelParams.dt
    t                           = collect(0:dt:(buffersize-1)*dt) #40

    Pcurvature = zeros(length(t),2)

    
    trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]        # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.nPolyCurvature = 4 # has to be 4 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs
    trackCoeff.nPolyXY = 6  # has to be 6 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs
    n_rounds = 1
    z_pred_log = zeros(mpcParams.N+1,4,length(t),n_rounds)
    u_pred_log = zeros(mpcParams.N,2,length(t),n_rounds)
    lambda_log = zeros(oldTraj.n_oldTraj,length(t),n_rounds)
    cost        = zeros(7,length(t),n_rounds)

    ssInfOn_log = zeros(oldTraj.n_oldTraj, length(t), n_rounds)
    curv_approx = zeros(mpcParams.N,length(t), n_rounds)
    pred_obst = zeros(mpcParams.N+1,3)

    if load_safeset == true
        SafeSetData = load(safeset)
        oldTraj = SafeSetData["oldTraj"]
        obstacle = SafeSetData["obstacle"]
    end

    j = 1
    z_final_x = zeros(1,4)::Array{Float64,2}
    u_final = zeros(1,2)
    for j=1:n_rounds #10
        

        lapStatus.currentLap = j
        tt          = zeros(length(t),1)
        tt1         = zeros(length(t),1)
        zCurr_s     = zeros(length(t),4)          # s, ey, epsi, v
        zCurr_x     = zeros(length(t),6)          # x, y, v_x,v_y, psi, psi_dot
        uCurr       = zeros(length(t),2)
        distance2obst = 1000*ones(length(t))
        curvature_curr = zeros(length(t))
        
        #T
        for k = oldTraj.n_oldTraj-1:-1:1
            obstacle.s_obstacle[:,k+1]  = obstacle.s_obstacle[:,k]
            obstacle.sy_obstacle[:,k+1] = obstacle.sy_obstacle[:,k]
            obstacle.v[:,k+1] = obstacle.v[:,k]
        end
        obstacle.s_obstacle[1,1] = s_obst_init
        obstacle.sy_obstacle[1,1] = sy_obst_init
        obstacle.v[1,1] = v_obst_init
        #setup point for vehicle on track in first round. gets overwritten in other rounds
        zCurr_x[1,1] = z_Init[1] # x = 1.81 for s = 32     14 in curve
        zCurr_x[1,2] = z_Init[2] # y = 2.505 for s = 32  12.6
        zCurr_x[1,3] = z_Init[3]
        zCurr_x[1,4] = z_Init[4]  # compare value to v_pathfollowing
        zCurr_x[1,5] = z_Init[5]
        zCurr_x[1,6] = z_Init[6]
        
        if j == 1 && load_safeset == true
            zCurr_x[1,1] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],1,1]#+convert(Float64,trackCoeff.ds)*0.5
            zCurr_x[1,2] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],2,1]
            zCurr_x[1,3] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],3,1]
            zCurr_x[1,4] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],4,1]
            zCurr_x[1,5] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],5,1]
            zCurr_x[1,6] =oldTraj.oldTrajXY[oldTraj.oldNIter[1],6,1]
        end
        if j>1             #setup point for vehicle after first round                   # if we are in the second or higher lap
            zCurr_x[1,:] = z_final_x
    
            # because the track is not closed we always set up the car at the same place each round
            #zCurr_x[1,1] = z_Init[1]#+convert(Float64,trackCoeff.ds)*0.5 # x = 1.81 for s = 32     14 in curve
            # zCurr_x[1,1] = z_Init[1]
            # zCurr_x[1,2] = z_Init[2] # y = 2.505 for s = 32  12.6
            # zCurr_x[1,5] = z_Init[5]
            # zCurr_x[1,6] = z_Init[6]
            # zCurr_x[1,3] = z_final_x[3]
            # zCurr_x[1,4] = z_final_x[4]
            # @show z_final_x[3]
            # @show z_final_x[4]
            # v_abs = sqrt(z_final_x[3]^.2+z_final_x[4]^.2)
            # zCurr_x[1,3] = v_abs * cos(z_Init[5])
            # zCurr_x[1,4] = v_abs * sin(z_Init[5])

            uCurr[1,:] = u_final
        end

            
        if j == 1 && load_safeset == false
            # path following cost in first round
            @NLobjective(m.mdl, Min, m.costPath + m.derivCost + m.controlCost + m.costObstacle)
        elseif j == 2 || (load_safeset == true && j == 1)
            #learning objective formulation, minimize the sum of all parts of the objective
            @NLobjective(m.mdl, Min, m.costZ + m.costZTerm + m.constZTerm + m.derivCost + m.controlCost + m.laneCost + m.costObstacle)
        end
        
        ###########iterations learning
        finished    = false
        i = 1
        while i<=length(t)-1 && !finished # as long as we have not reached the maximal iteration time for one round or ended the round
        
            # to make it work s start has to grow over time actual it is just always at 0
       
            # the argument i in localizeVehicleCurvAbs  is solely used for debugging purposes plots not needed for control
            # localize takes some time see ProfileView.view()
            tic()
            zCurr_s[i,:], trackCoeff.coeffCurvature = localizeVehicleCurvAbs(zCurr_x[i,:],x_track,y_track,trackCoeff, i, mpcParams.N, modelParams.dt, Pcurvature)
            if i == 1 && zCurr_s[1,1] > 2
                warn("closest point was before finish line, forced s =0")
                println("x:$(zCurr_x[i,1]), y:$(zCurr_x[i,2])")
                zCurr_s[1,1]= 0
            end
            posInfo.s   = zCurr_s[i,1]
            curvature_curr[i] = trackCoeff.coeffCurvature[1]*posInfo.s^4+trackCoeff.coeffCurvature[2]*posInfo.s^3+trackCoeff.coeffCurvature[3]*posInfo.s^2+trackCoeff.coeffCurvature[4]*posInfo.s +trackCoeff.coeffCurvature[5]
            distance2obst[i] = (obstacle.s_obstacle[i,1]-obstacle.rs) - posInfo.s
            #t_absci =toq()
            #if the car has crossed the finish line
            if zCurr_s[i,1] >= posInfo.s_target
                println("Reached finish line at step $i")
                finished = true
                #we count up here as the first round ends just as we would do if the loop get terminiated because i is >= length(t). we count it down later on again to get right index
                i = i + 1
                lapStatus.currentIt = i
                break
            elseif i >1 && zCurr_s[i-1,1] >= posInfo.s_target-1 && zCurr_s[i,1]<1 # on a closed trakc the finish line lies close to the startinf point. it can be that the discrete calculated s postion is bigger than the end of the track and already in the new round. to still detect the finish line we added this statement
                println("alternative detection of finish line")
                println("Reached finish line at step $i")
                finished = true
                #we count up here as the first round ends just as we would do if the loop get terminiated because i is >= length(t). we count it down later on again to get right index
                i = i + 1
                lapStatus.currentIt = i
                break
            end
            #!! println("s = $(zCurr_s[i,1])")
            

            ##################
            pred_obst = predictObstaclePos(obstacle, modelParams, mpcParams, i)
            ##################
            
            if j > 1 || load_safeset == true
                tic()
                addOldtoNewPos(oldTraj, distance2obst[i],obstacle,i,pred_obst, mpcParams,zCurr_s,modelParams.dt,mpcCoeff)
                deleteInfeasibleTrajectories!(oldTraj,posInfo,obstacle, pred_obst, i, zCurr_s,modelParams.dt)
		        
                tt1[i] = toq()
                #println(" time to add/remove traj $(tt1[i])")
                coeffConstraintCost!(oldTraj,mpcCoeff,posInfo,mpcParams,i)
            end
            tic()
            #solve with zCurr_s containing s ey values 

            #####warm start
            # if i ==1 && j == 1 && load_safeset == false
            #     for k = 1:mpcParams.N+1
            #         mpcSol.z[k,1] = z_Init[4]*k*dt+0.32
            #         mpcSol.z[k,4] = z_Init[4]
            #     end
            # elseif i==1 && j>1
            #     mpcSol.z[:,1] = mpcSol.z[:,1]-mpcSol.z[1,1]+0.32
            #     mpcSol.lambda = zeros(oldTraj.n_oldTraj)
            #     mpcSol.lambda[1] = 1
            # end
            # setvalue(m.u_Ol, mpcSol.u)
            # setvalue(m.z_Ol, mpcSol.z)
            # setvalue(m.lambda, mpcSol.lambda)
            # setvalue(m.eps, mpcSol.eps)
            #####################
          
            solveMpcProblem!(m,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:]',[mpcSol.a_x;mpcSol.d_f], pred_obst,i)
            tt[i]       = toq()
            setvalue(m.ssInfOn,ones(oldTraj.n_oldTraj))# reset the all trajectories to on position

            uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f]
            #have Zcurr as states XY and simulate from there return XY values of states 
            zCurr_x[i+1,:]  = simModel_exact_dyn_x(zCurr_x[i,:],uCurr[i,:],modelParams.dt,modelParams) #!! @show

            #update Position of the Obstacle car        
            computeObstaclePos!(obstacle, dt, i, x_track, trackCoeff, zCurr_s[i,1]) #this funciton computes values for row i+1
            

            cost[:,i,j]         = mpcSol.cost
            lambda_log[:,i,j]   = mpcSol.lambda
            z_pred_log[:,:,i,j] = mpcSol.z
            u_pred_log[:,:,i,j] = mpcSol.u
            ssInfOn_log[:,i,j]  = mpcSol.ssInfOn
            curv_approx[:,i,j]  = getvalue(m.c)
            

            tt2= toq()
            if i%50 == 0 
                println(" Time: $(tt[i]) s, Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus)")
                #println(" Time for whole iteration: $(tt2) s")
                # if j > 1 || load_safeset == true
                #     println("ssInfOn time: $tt1 s")
                # end
                # println("calculate abs: $t_absci s")
                # println("get curve-approx = $t_curv")
            end
            # if tt[i] >0.08 #if solving takes long
            #     println(" Time: $(tt[i]) s, Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus)")
            #     println(" Time: $(tt2) s for the whole step of the loop ")
            #     println("                           ")
            # end

            

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
        tic()
        saveOldTraj(oldTraj,zCurr_s, zCurr_x,uCurr,lapStatus,buffersize,modelParams.dt, load_safeset, cost[:,:,j],
          lambda_log[:,:,j],z_pred_log[:,:,:,j],u_pred_log[:,:,:,j],ssInfOn_log[:,:,j], mpcSol, obstacle,distance2obst, curvature_curr)
        tt3= toq()
        println("Max time to calculate infeasible traj: $(maximum(tt1))")
        println(" Time to save and overwrite trajectories: $(tt3) s")
        ###############

        oldTraj.oldNIter[1] = i
        if j>1 && oldTraj.cost2Target[1,2] <= oldTraj.cost2Target[1,1]
            warn("round was not faster. cost : $(oldTraj.cost2Target[1,1])")

            println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        else
            println("cost of round: $(oldTraj.cost2Target[1,1])")
        end

        println("*************************************************************************")


        # println("Press c to cancel MPC")
        # println("Press Enter for next round solving")
        # a = ' '
        # a = readline()
        # if a == "c\r\n" 
        #         break
        # end
        # figure()
        # plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],oldTraj.curvature[1:oldTraj.oldNIter[j],j], color = "red")
        # plot(Pcurvature[1:oldTraj.oldNIter[j],1], Pcurvature[1:oldTraj.oldNIter[j],2], color = "green")
    end
    n_rounds = j #update n_rounds to represent actual number of simualated rounds

    #filename = string("../LMPCdata/"string(Dates.today()),"-Data.jld")
    #### alternative numbering to generate results to keep 
    filename = string("data/"string(Dates.today()),"-",Dates.format(now(), "HH-MM"),"-Data.jld")
    if isfile(filename)
        filename = string("data/"string(Dates.today()),"-",Dates.format(now(), "HH-MM"),"-Data-2.jld")
        warn("File already exists. Added extension \"-2\" ")
    end
    println("Save data to $filename .......")
    jldopen(filename, "w") do file
        #addrequire(file, classes) #ensures that custom data types are working when loaded

        JLD.write(file, "x_track", x_track)
        JLD.write(file, "y_track", y_track)
        JLD.write(file, "trackCoeff", trackCoeff)
        JLD.write(file, "obstacle", obstacle)
        JLD.write(file, "modelParams", modelParams)
        JLD.write(file, "mpcParams", mpcParams)
        JLD.write(file, "buffersize", buffersize)
        JLD.write(file, "curv_approx", curv_approx) ###not same numberign as in oldTraj
        JLD.write(file, "oldTraj", oldTraj)
        JLD.write(file, "mpcCoeff",mpcCoeff)

    end
    println("finished")
