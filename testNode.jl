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

#just loads one specified track with distances betwwen points =1 m and returnx x and y coordinatates
println("loadMap.......")
include("helper/loadTestMap.jl")
x_track, y_track = loadTestMap() #load a track to test controller
println("loaded")

# Load Variables and create Model:
#println("Loading and defining variables...")
#include("helper/createModel.jl")

# Initialize model by solving it once #?? why don do the initial solve?
#println("Initial solve...")
#solve(mdl)

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

    z_Init = zeros(4)# used in InitializeModel function
    z_Init[1] = 0 # x = 1.81 for s = 32     14 in curve
    z_Init[2] = 0# y = 2.505 for s = 32  12.6
    z_Init[3] = 0.94
    z_Init[4]  = 0.4

    
  
    InitializeParameters(mpcParams,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,buffersize)
    ##define obstacle x and xy vlaues not used at the moment 
    #for a clean definition of the x,y points the value of s_obstacle has to be the same as one of the points of the source map. 
    # the end semi axes are approximated over the secant of the points of the track. drawing might not be 100% accurate
    obstacle.s_obstacle = 25
    obstacle.sy_obstacle = 0.1
    obstacle.rs = 0.5
    obstacle.ry = 0.2

    # all thes values are currently just used for plotting the semi axes of the obstacle
    obstacle.index = obstacle.s_obstacle/trackCoeff.ds+1 
    obstacle.xy_vector = [x_track[obstacle.index]; y_track[obstacle.index]]
    x_secant = x_track[obstacle.index+1] - x_track[obstacle.index-1]   
    y_secant = y_track[obstacle.index+1] - y_track[obstacle.index-1]
    secant_vec = [ x_secant ; y_secant]/norm([ x_secant ; y_secant])
    normal_vec = [-y_secant ; x_secant]/norm([-y_secant ; x_secant])
    vector_sy = normal_vec*obstacle.sy_obstacle
    obstacle.xy_vector = obstacle.xy_vector + vector_sy
    #points to plot the ry semi axis of the ellipsis
    vector_ry = normal_vec*obstacle.ry
    ey_up = obstacle.xy_vector + vector_ry
    ey_down = obstacle.xy_vector - vector_ry
    #points to plot the rs semi axis of the ellipsis
    vector_rs = secant_vec*obstacle.rs
    es_up = obstacle.xy_vector + vector_rs
    es_down = obstacle.xy_vector - vector_rs

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


    posInfo.s_start             = 0.0 #?? should be changed probably in the localizeabs function
    posInfo.s_target            = 30.2 #has to be fitted to track , current test track form ugo has 113.2 meters
    trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]        # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.nPolyCurvature = 4 # has to be 4 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs
    trackCoeff.nPolyXY = 6  # has to be 6 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs

    z_final_x = zeros(1,4)::Array{Float64,2}
    u_final = zeros(1,2)::Array{Float64,2}

   

     #T
    
    z_log = zeros(11,4,length(t))
    u_log = zeros(10,2,length(t))
    i_final= 100000000 # high value so real values will be smaller
   
    
    # figure(1)
    # plot(x_track',y_track')
    for j=1:2 #10
        lapStatus.currentLap = j
        tt          = zeros(length(t),1)
        zCurr_s                     = zeros(length(t)+1,4)          # s, ey, epsi, v
        zCurr_x                     = zeros(length(t)+1,4)          # x, y, psi, v
        uCurr       = zeros(length(t),2)
        ParIntLog = zeros(length(t))
        #T
      
        #setup point for vehicle on track in first round. gets overwritten in other rounds
        zCurr_x[1,1] = 0 # x = 1.81 for s = 32     14 in curve
        zCurr_x[1,2] = 0 # y = 2.505 for s = 32  12.6
        zCurr_x[1,3] = 0.94
        zCurr_x[1,4] = 0.4  #?? initialize with v_ref v_pathfollowing ?
        
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
        while i<=length(t) && !finished # as long as we have not reached the maximal iteration time for one round or ended the round? #?? error if i > length(t)?
        
           
            
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
                coeffConstraintCost(oldTraj,mpcCoeff,posInfo,mpcParams)
                tt1 = toq()
                #println("coeffConstr: $tt1 s")
            end
            
            tic()
	      
            #solve with zCurr_s containing s ey values 
            
            solveMpcProblem(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:]',uCurr[i,:]', obstacle)
        
                
            tt[i]       = toq()
            cost[i,:]   = mpcSol.cost
            ParIntLog[i] = mpcSol.ParInt[1]
            uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f]
            z_log[:,:,i] = mpcSol.z
            u_log[:,:,i] = mpcSol.u
          
            #have Zcurr as states XY and simulate from there return XY values of states 
            zCurr_x[i+1,:]  = simModel_x(zCurr_x[i,:],uCurr[i,:],modelParams.dt,modelParams) #!! @show
            
            if i%50 == 0 
                println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")
            end

           ###T
            # if i%20 == 0 
              
            # #   plot(xs,ys) # tangent to current point  
            #     if j == 1
            #         figure(1)
            #         scatter(zCurr_x[i+1,1], zCurr_x[i+1,2], color = "red")
            #     end
            #     if j == 2
            #         figure(1)
            #         scatter(zCurr_x[i+1,1], zCurr_x[i+1,2], color = "green")
            #     end
            #     if j >= 3
            #         figure(1)
            #         scatter(zCurr_x[i+1,1], zCurr_x[i+1,2], color = "blue")
            #     end
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

          # Print results
        # --------------------------------


        #### plot states and cost
        figure(2)
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
        plot(t,cost[:,1],"r",t,cost[:,2],"g",t,cost[:,3],"b",t,cost[:,4],"y",t,cost[:,5],"m",t,cost[:,6],"c", t, cost[:,7])
        grid(1)
        title("Cost distribution")
        legend(["z","z_Term","z_Term_const","deriv","control","lane"])

        i_final_old = i_final
        i_final = i
        if i_final_old <= i_final
            println("round was not faster. no learning")
            println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        end
       

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

        #x-y Plot of the racetrack and obstacle
        figure(8)
        clf()
        plot(x_track',y_track', linestyle="--", color = "yellow", linewidth = 0.5)#plot the racetrack
        plot(zCurr_x[1:i_final,1], zCurr_x[1:i_final,2], color = "blue", linewidth = 1)# plot the trajectory of the car
        scatter(obstacle.xy_vector[1], obstacle.xy_vector[2],color = "red")#plot the center of the obstacle
        # plot the two boundary lines
        plot(boundary_up[1,:], boundary_up[2,:],color="black",linestyle="--")
        plot(boundary_down[1,:], boundary_down[2,:],color="black",linestyle="--")
        #plot the obstacle  semi axes
        ax4= gca()
        plot([ey_up[1],ey_down[1]],[ey_up[2],ey_down[2]],color = "red")
        plot([es_up[1],es_down[1]],[es_up[2],es_down[2]],color = "red")
        #ax4[:set_xlim]([5.7,5.9])
        #ax4[:set_ylim]([8.0,8.2])
        grid(0.5)

        if j >= 2
            
            if j >= 3
                figure(7)
                clf()
                plot(t,ParIntLog)
                xlabel("t in [s]")
                ylabel("ParInt")
            end
            for i = 1:i_final


                #plot predicted states over s
                figure(4)   
                clf()
                ax4=subplot(311)
                plot(oldTraj.oldTraj[1:i_final,1,1], oldTraj.oldTraj[1:i_final,4,1], color= "blue")
                plot(oldTraj.oldTraj[1:i_final_old,1,2], oldTraj.oldTraj[1:i_final_old,4,2], color= "yellow")
                grid()
                xlabel("s in [m]")
                ylabel("v in [m/s]")
                legend(["v current round","v old round"])
                plot(z_log[:,1,i], z_log[:,4,i] ,marker="o")
                ax5=subplot(312, sharex=ax4)
                plot(oldTraj.oldTraj[1:i_final,1,1], oldTraj.oldTraj[1:i_final,2,1], color= "blue")
                plot(oldTraj.oldTraj[1:i_final_old,1,2], oldTraj.oldTraj[1:i_final_old,2,2], color= "yellow")
                grid()
                xlabel("s in [m]")
                ylabel("e_y in [m]")
                legend(["e_y current round","e_y old round"])
                plot(z_log[:,1,i], z_log[:,2,i] ,marker="o")
                ax6=subplot(313, sharex=ax4)
                plot(oldTraj.oldTraj[1:i_final,1,1], oldTraj.oldTraj[1:i_final,3,1], color= "blue")
                plot(oldTraj.oldTraj[1:i_final_old,1,2], oldTraj.oldTraj[1:i_final_old,3,2], color= "yellow")
                grid()
                xlabel("s in [m]")
                ylabel("e_psi in [rad]")
                legend(["e_psi current round","e_psi old round"])
                plot(z_log[:,1,i], z_log[:,3,i] ,marker="o")
                
               

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

                i = i+1



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