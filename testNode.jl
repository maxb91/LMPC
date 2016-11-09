module testNode

using JuMP
using Ipopt
using PyPlot

include("helper/classes.jl")
include("helper/functions.jl")
include("helper/coeffConstraintCost.jl")
include("helper/solveMpcProblem.jl")
include("helper/simModel.jl")
include("helper/localizeVehicleCurvAbs.jl")
include("helper/computeObstaclePos.jl")
include("helper/calculateObstacleXY.jl")

#just loads one specified track with distances betwwen points =1 m and returnx x and y coordinatates
println("loadMap.......")
include("helper/loadTestMap.jl")
x_track, y_track = loadTestMap() #load a track to test controller
println("loaded")

# Load Variables and create Model:
#println("Loading and defining variables...")
#include("helper/createModel.jl")


function run_sim()
    # DEFINE PARAMETERS
    # Define and initialize variables
    println("define types ........")
    oldTraj                     = OldTrajectory()
    posInfo                     = PosInfo()
    mpcCoeff                    = MpcCoeff()
    lapStatus                   = LapStatus(1,1)
    mpcSol                      = MpcSol()
    obstacle                    = Obstacle()
    trackCoeff                  = TrackCoeff()      # info about track (at current position, approximated)
    modelParams                 = ModelParams()
    mpcParams                   = MpcParams()
    mdl                         = MpcModel()

    buffersize                  = 800
    close("all")

    z_Init = zeros(4)# used in InitializeModel function not to setup problem
    z_Init[1] = 0 # x = 1.81 for s = 32     14 in curve
    z_Init[2] = 0# y = 2.505 for s = 32  12.6
    z_Init[3] = 0.94
    z_Init[4]  = 0.4

    
  
    InitializeParameters(mpcParams,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,obstacle,buffersize)
    ##define obstacle x and xy vlaues not used at the moment 
    #for a clean definition of the x,y points the value of s_obstacle has to be the same as one of the points of the source map. 
    # the end semi axes are approximated over the secant of the points of the track. drawing might not be 100% accurate
    s_obst_init = 5.0 
    sy_obst_init = 0.1
    obstacle.s_obstacle[1] = s_obst_init#gets overwritten in loop (moving)
    obstacle.sy_obstacle[1]  = sy_obst_init#gets overwritten in loop (moving)
    obstacle.rs = 0.5
    obstacle.ry = 0.2

    ##this part is to calculate the tracks boundaries and plot them later
    xy_track  = [x_track; y_track]
    trackL = size(xy_track,2)
    boundary_up = zeros(2,trackL)
    boundary_down = zeros(2,trackL)
    for kkk = 1:1:trackL
        if 1< kkk < trackL 
            xt_secant = x_track[kkk+1] - x_track[kkk-1]
            yt_secant = y_track[kkk+1] - y_track[kkk-1]
            normVec = [-yt_secant; xt_secant]/norm([-yt_secant; xt_secant])
            boundary_up[:,kkk] = xy_track[:,kkk] + normVec * trackCoeff.width /2
            boundary_down[:,kkk] = xy_track[:,kkk] - normVec * trackCoeff.width /2
        elseif kkk == 1
            xt_secant = x_track[kkk+1] - x_track[kkk]
            yt_secant = y_track[kkk+1] - y_track[kkk]
            normVec = [-yt_secant; xt_secant]/norm([-yt_secant; xt_secant])
            boundary_up[:,kkk] = xy_track[:,kkk] + normVec * trackCoeff.width /2
            boundary_down[:,kkk] = xy_track[:,kkk] - normVec * trackCoeff.width /2
        else
            xt_secant = x_track[kkk] - x_track[kkk-1]
            yt_secant = y_track[kkk] - y_track[kkk-1]
            normVec = [-yt_secant; xt_secant]/norm([-yt_secant; xt_secant])
            boundary_up[:,kkk] = xy_track[:,kkk] + normVec * trackCoeff.width /2
            boundary_down[:,kkk] = xy_track[:,kkk] - normVec * trackCoeff.width /2
        end
    end
    ######

    println("Initialize Model........")
    InitializeModel(mdl,mpcParams,modelParams,trackCoeff,z_Init, obstacle)
    println("Initial solve........")
    solve(mdl.mdl)
    println("Initial solve done!")
    println("*******************************************************")
    println("*******************************************************")
    # Simulation parameters
    dt                          = modelParams.dt
    t                           = collect(0:dt:60) #40
    cost                        = zeros(length(t),6)


    posInfo.s_start             = 0.0 #does not get changed with the current version
    posInfo.s_target            = 25.2 #has to be fitted to track , current test track form ugo has 113.2 meters
    trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]        # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.nPolyCurvature = 4 # has to be 4 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs
    trackCoeff.nPolyXY = 6  # has to be 6 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs

    z_final_x = zeros(1,4)::Array{Float64,2}
    u_final = zeros(1,2)::Array{Float64,2}

   

     #T
    
    z_log = zeros(11,4,length(t))
    u_log = zeros(10,2,length(t))
    i_final= 100000000 # high value so real values will be smaller
   
    

    for j=1:10 #10
        

        lapStatus.currentLap = j
        tt          = zeros(length(t),1)
        zCurr_s     = zeros(length(t)+1,4)          # s, ey, epsi, v
        zCurr_x     = zeros(length(t)+1,4)          # x, y, psi, v
        uCurr       = zeros(length(t),2)
        ParIntLog = zeros(length(t))
        #T
      
        obstacle.s_obstacle[1] = s_obst_init
        obstacle.sy_obstacle[1] = sy_obst_init
        #setup point for vehicle on track in first round. gets overwritten in other rounds
        zCurr_x[1,1] = 0 # x = 1.81 for s = 32     14 in curve
        zCurr_x[1,2] = 0 # y = 2.505 for s = 32  12.6
        zCurr_x[1,3] = 0.94
        zCurr_x[1,4] = 0.4  # compare value to v_pathfollowing
        
        if j>1               #setup point for vehicle after first round                   # if we are in the second or higher lap
            zCurr_x[1,:] = z_final_x
            # because the track is not closed we always set up the car at the same place each round
            zCurr_x[1,1] = 0 # x = 1.81 for s = 32     14 in curve
            zCurr_x[1,2] = 0 # y = 2.505 for s = 32  12.6
            zCurr_x[1,3] = 0.94
            uCurr[1,:] = u_final
            
        end
        cost        = zeros(length(t),7)
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
                calculateObstacleXY!(obstacle, trackCoeff, xy_track, i)# caluclate the postion of the obstacle at the end
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
            solveMpcProblem!(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:]',uCurr[i,:]', obstacle,i)
        
                
            tt[i]       = toq()
            cost[i,:]   = mpcSol.cost
            ParIntLog[i] = mpcSol.ParInt[1]
            uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f]
            z_log[:,:,i] = mpcSol.z
            u_log[:,:,i] = mpcSol.u
          
            #have Zcurr as states XY and simulate from there return XY values of states 
            zCurr_x[i+1,:]  = simModel_x(zCurr_x[i,:],uCurr[i,:],modelParams.dt,modelParams) #!! @show

            #update Position of the Obstacle car
            #!!
            #the computation in the XY plane is in a seperate function as these values are just used for the visualization and could be calculated in post procssing if necessary for real time
            calculateObstacleXY!(obstacle, trackCoeff, xy_track, i) #this funciton computes values for row i
            computeObstaclePos!(obstacle, dt, i)#this funciton computes values for row i+1
            
            
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


        i_final_old = i_final
        i_final = i
        if i_final_old <= i_final
            println("round was not faster. no learning")
            println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        end

          # Print results
        # --------------------------------
        #########################################################
        #########################################################
        #### plot states and cost
        figure(1)
        ax1=subplot(311)
        plot(t,zCurr_s[1:end-1,1],"y",t,zCurr_s[1:end-1,2],"r",t,zCurr_s[1:end-1,3],"g",t,zCurr_s[1:end-1,4],"b")
        grid(1)
        legend(["s","eY","ePsi","v"])
        title("States")
        ax2=subplot(312,sharex=ax1)
        plot(t,uCurr[:,1],"r",t,uCurr[:,2],"g")
        grid(1)
        title("Control input")
        legend(["a","d_f"])
        ax3=subplot(313,sharex=ax1)
        plot(t,cost[:,1],"r",t,cost[:,3],"b",t,cost[:,4],"y",t,cost[:,5],"m",t,cost[:,6],"c", t, cost[:,7])
        grid(1)
        title("Cost distribution")
        legend(["z","z_Term_const","deriv","control","lane", "Obstacle"])

        
        figure(2)
        clf()
        ax9= subplot(3,1,1)
        ax9[:plot](oldTraj.oldTraj[1:i_final,1,1],cost[1:i_final,2])  
        grid()
        xlabel("s in [m]")
        ylabel("Terminal cost ") 

        ax7= subplot(3,1,2,sharex=ax9)
        ax7[:plot](oldTraj.oldTraj[1:i_final,1,1],cost[1:i_final,3])  
        grid()
        xlabel("s in [m]")
        ylabel("cost constraint")    
        ax8= subplot(3,1,3,sharex=ax9)
        ax8[:plot](oldTraj.oldTraj[1:i_final,1,1],cost[1:i_final,7])  
        grid()
        xlabel("s in [m]")
        ylabel("cost Obstacle")     

         # ####T
        if j == 1
           
            # figure(4)
            # for i = 1:i_final
            #     clf()
            #     ax1=subplot(311)
            #     plot(zCurr_s[1:i_final,1],zCurr_s[1:i_final,4])
            #     ax1=subplot(312)
            #     plot(zCurr_s[1:i_final,1],zCurr_s[1:i_final,2])
            #     ax1=subplot(313)
            #     plot(zCurr_s[1:i_final,1],zCurr_s[1:i_final,3])
            #     subplot(311)
            #     plot(z_log[:,1,i], z_log[:,4,i] ,marker="o")
            #     subplot(312)
            #     plot(z_log[:,1,i], z_log[:,2,i] ,marker="o")
            #     subplot(313)
            #     plot(z_log[:,1,i], z_log[:,3,i] ,marker="o")
            #     i = i+1
            #     println("Press Enter for next plot step")
            #     println("Press c to cancel plot")
            #     a = ' '
            #     a = readline()
            #     if a == "c\r\n" 
            #         break
            #     end
            # end
        ####endT
        end

        #x-y Plot of the racetrack noundaries and initilaize obstacle and car
        figure(8)
        clf()
        ax10= subplot(1,1,1)
        ax10[:plot](x_track',y_track', linestyle="--", color = "yellow", linewidth = 0.5)#plot the racetrack
        car_plot = ax10[:plot](zCurr_x[1,1], zCurr_x[1,2], color = "blue",marker="o")
        obstacle_plot = ax10[:plot](obstacle.xy_vector[1,1], obstacle.xy_vector[1,2], color = "red",marker="o")
        y_obst_plot = ax10[:plot]([obstacle.axis_y_up[1,1],obstacle.axis_y_down[1,1]],[obstacle.axis_y_up[1,2],obstacle.axis_y_down[1,2]],color = "red")#plot the y semi axis
        s_obst_plot = ax10[:plot]([obstacle.axis_s_up[1,1],obstacle.axis_s_down[1,1]],[obstacle.axis_s_up[1,2],obstacle.axis_s_down[1,2]],color = "red")# plot the s semi axis
        ax10[:plot](boundary_up[1,:], boundary_up[2,:],color="green",linestyle="--")
        ax10[:plot](boundary_down[1,:], boundary_down[2,:],color="green",linestyle="--")
        grid()     
        ##########
        
       #not everything necessary if we do the live plot
       #plot of trajectory of car and obstacle
        # ax10[:plot](zCurr_x[1:i_final,1], zCurr_x[1:i_final,2], color = "blue", linewidth = 1)# plot the trajectory of the car
        # ax10[:plot](obstacle.xy_vector[1:i_final,1], obstacle.xy_vector[1:i_final,2], color = "red")
        # car_plot[1][:remove]()
        # obstacle_plot[1][:remove]()
        # car_plot = ax10[:plot](zCurr_x[i_final,1], zCurr_x[i_final,2], color = "blue", marker="o")
        # obstacle_plot = ax10[:plot](obstacle.xy_vector[i_final,1], obstacle.xy_vector[i_final,2], color = "red", marker="o")  
        
        # plot the two boundary lines
        # ax10[:plot](boundary_up[1,:], boundary_up[2,:],color="green",linestyle="--")
        # ax10[:plot](boundary_down[1,:], boundary_down[2,:],color="green",linestyle="--")
        #plot the obstacle  semi axes
        # ax4 = gca()
        # scatter(obstacle.xy_vector[1,1], obstacle.xy_vector[1,2],color = "red")#plot the center of the obstacle
        # plot([ey_up[1,1],ey_down[1,1]],[ey_up[1,2],ey_down[1,2]],color = "red")#plot the y semi axis
        # plot([es_up[1,1],es_down[1,1]],[es_up[1,2],es_down[1,2]],color = "red")# plot the s semi axis
        # #ax4[:set_xlim]([5.7,5.9])
        # #ax4[:set_ylim]([8.0,8.2])
        # grid(0.5)

        if j >= 5
            
            if j >= 3
                figure(7)
                clf()
                plot(t,ParIntLog)
                xlabel("t in [s]")
                ylabel("ParInt")
            end

            figure(4)   
            clf()
            ax4 = subplot(311)
            plot(oldTraj.oldTraj[1:i_final,1,1], oldTraj.oldTraj[1:i_final,4,1], color= "blue")
            plot(oldTraj.oldTraj[1:i_final_old,1,2], oldTraj.oldTraj[1:i_final_old,4,2], color= "yellow")
            grid()
            xlabel("s in [m]")
            ylabel("v in [m/s]")
            legend(["v current round","v old round"])
            p1 = plot(1,1)

            ax5 = subplot(312, sharex=ax4)
            plot(oldTraj.oldTraj[1:i_final,1,1], oldTraj.oldTraj[1:i_final,2,1], color= "blue")
            plot(oldTraj.oldTraj[1:i_final_old,1,2], oldTraj.oldTraj[1:i_final_old,2,2], color= "yellow")
            grid()
            xlabel("s in [m]")
            ylabel("e_y in [m]")
            legend(["e_y current round","e_y old round"])
            p2 = ax5[:plot](1,1)

            ax6 = subplot(313, sharex=ax4)
            hold(true)
            plot(oldTraj.oldTraj[1:i_final,1,1], oldTraj.oldTraj[1:i_final,3,1], color= "blue")
            plot(oldTraj.oldTraj[1:i_final_old,1,2], oldTraj.oldTraj[1:i_final_old,3,2], color= "yellow")
            grid()
            xlabel("s in [m]")
            ylabel("e_psi in [rad]")
            legend(["e_psi current round","e_psi old round"])
            p3 = ax6[:plot](1,1)

            
            
            
            for i = 1:i_final


                #plot predicted states over s
              
                p1[1][:remove]()
                p1 = ax4[:plot](z_log[:,1,i], z_log[:,4,i] ,marker="o") #plot predicted velocity


                p2[1][:remove]()
                p2 = ax5[:plot](z_log[:,1,i], z_log[:,2,i] ,marker="o") #plot predicted e_y
             
                p3[1][:remove]()
                p3 = ax6[:plot](z_log[:,1,i], z_log[:,3,i] ,marker="o")  #plot predicted e_psi
                
               
               #x-y plot
                car_plot[1][:remove]()
                obstacle_plot[1][:remove]()
                y_obst_plot[1][:remove]()
                s_obst_plot[1][:remove]()
                ax10[:plot](zCurr_x[1:i,1], zCurr_x[1:i,2], color = "blue")
                car_plot = ax10[:plot](zCurr_x[i,1], zCurr_x[i,2], color = "blue", marker="o")
                ax10[:plot](obstacle.xy_vector[1:i,1], obstacle.xy_vector[1:i,2], color = "red")
                obstacle_plot = ax10[:plot](obstacle.xy_vector[i,1], obstacle.xy_vector[i,2], color = "red", marker="o")  
                y_obst_plot = ax10[:plot]([obstacle.axis_y_up[i,1],obstacle.axis_y_down[i,1]],[obstacle.axis_y_up[i,2],obstacle.axis_y_down[i,2]],color = "red")#plot the y semi axis
                s_obst_plot = ax10[:plot]([obstacle.axis_s_up[i,1],obstacle.axis_s_down[i,1]],[obstacle.axis_s_up[i,2],obstacle.axis_s_down[i,2]],color = "red")# plot the s semi axis


                # #plot inputs
                # figure(5)
                # clf()
                # subplot(211) 
                # plot(oldTraj.oldTraj[1:i_final-1,1,1], oldTraj.oldInput[1:i_final-1,1,1], color= "green")
                # plot(z_log[1:end-1,1,i], u_log[:,1,i] ,marker="o")
                # xlabel("s in [m]")
                # ylabel("acceleration in [m/s^2]")
                # legend(["applied acceleration","predicted acceleration"])
                # subplot(212)
                #  plot(oldTraj.oldTraj[1:i_final-1,1,1], oldTraj.oldInput[1:i_final-1,2,1], color= "green")  
                # plot(z_log[1:end-1,1,i], u_log[:,2,i] ,marker="o")  
                #  xlabel("s in [m]")
                # ylabel("steering angle in [rad]")
                # legend(["applied steering angle","predicted steering angle"]) 


            
                # figure(6)
                # plot(zCurr_s[1:i_final,1], zCurr_x[1:i_final,3])
                # xlabel("s in [m]")
                # ylabel("absolute angle psi in [rad]")
                # legend(["absolute angle psi plotted over s"]) 

                i = i+10



                println("Press Enter for next plot step")
                println("Press c to cancel plot")
                a = ' '
                a = readline()
                if a == "c\r\n" 
                     break
                end
            end
            
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

end

end # end of module. we declare the main fail as a module to reduce warnings when including the file again