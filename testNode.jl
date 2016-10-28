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
include("helper/loadTestMap.jl")
x_track, y_track = loadTestMap() #load a track to test controller


# Load Variables and create Model:
#println("Loading and defining variables...")
#include("helper/createModel.jl")

# Initialize model by solving it once #?? why don do the initial solve?
#println("Initial solve...")
#solve(mdl)

function run_sim()
    # DEFINE PARAMETERS
    # Define and initialize variables
    oldTraj                     = OldTrajectory()
    posInfo                     = PosInfo()
    mpcCoeff                    = MpcCoeff()
    lapStatus                   = LapStatus(1,1)
    mpcSol                      = MpcSol()
    trackCoeff                  = TrackCoeff()      # info about track (at current position, approximated)
    modelParams                 = ModelParams()
    mpcParams                   = MpcParams()
    mdl                         = MpcModel()

    buffersize                  = 700

    z_Init = zeros(4)
    z_Init[1] = 140 # x = 1.81 for s = 32     14 in curve
    z_Init[2] = 126# y = 2.505 for s = 32  12.6
    z_Init[3] = -0.9
    z_Init[4]  = 2 

    InitializeParameters(mpcParams,trackCoeff,modelParams,posInfo,oldTraj,mpcCoeff,lapStatus,buffersize)
    InitializeModel(mdl,mpcParams,modelParams,trackCoeff,z_Init)

    # Simulation parameters
    dt                          = modelParams.dt
    t                           = collect(0:dt:40) #40
    zCurr_s                     = zeros(length(t)+1,4)          # s, ey, epsi, v
    zCurr_x                     = zeros(length(t)+1,4)          # x, y, psi, v
    uCurr                       = zeros(length(t),2)
    cost                        = zeros(length(t),6)

    #changeMetersforBARC
    posInfo.s_start             = 0.0 #?? should be changed probably in the localizeabs function
    posInfo.s_target            = 252 #has to be fitted to track , current test track form ugo has 113.2 meters
    trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]        # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.nPolyCurvature = 4 # has to be 4 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs
    trackCoeff.nPolyXY = 6  # has to be 6 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs

    z_final_x = zeros(1,4)::Array{Float64,2}
    u_final = zeros(1,2)::Array{Float64,2}

     #T
    z_log = zeros(11,4,length(t))
    u_log = zeros(10,2,length(t))
    i_final_1st = 0
   
    
     # figure(1)
     # plot(x_track',y_track')
    for j=1:4 #10
        lapStatus.currentLap = j
        tt          = zeros(length(t),1)
        zCurr_x       = zeros(length(t)+1,4)
        #T
      
        #changeMetersforBARC
        zCurr_x[1,1] = 140 # x = 1.81 for s = 32     14 in curve
        zCurr_x[1,2] = 126 # y = 2.505 for s = 32  12.6
        zCurr_x[1,3] = -0.9
        
        zCurr_x[1,4]  = 2  #?? initialize with v_ref v_pathfollowing ?
        uCurr       = zeros(length(t),2)
        if j>1                                  # if we are in the second or higher lap
            zCurr_x[1,:] = z_final_x
            # because the track is not closed we always set up the car at the same place each round
            zCurr_x[1,1] = 140 # x = 1.81 for s = 32     14 in curve
            zCurr_x[1,2] = 126 # y = 2.505 for s = 32  12.6
            zCurr_x[1,3] = -0.9
            uCurr[1,:] = u_final
            
        end
        cost        = zeros(length(t),6)
        finished    = false

        i = 1
        while i<=length(t) && !finished # as long as we have not reached the maximal iteration time for one round or ended the round? #?? error if i > length(t)?
            #!! set coeff according to a track elsewhere
           
            #!!todo func z curr oaut of xy
            # to make it work s start has to grow over time actual it is just always at 0
           #!!see what is defined in mpc params arguments vall there?
             # the argument i in localizeVehicleCurvAbs  is solely used for debugging purposes plots not needed for control
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
            println("s = $(zCurr_s[i,1])")
     
            posInfo.s   = zCurr_s[i,1]
            if j > 1
                       tic()
                coeffConstraintCost(oldTraj,mpcCoeff,posInfo,mpcParams)
                tt1 = toq()
                println("coeffConstr: $tt1 s")
            end
            
            tic()
	      
          #!!solve with zCurr_s containing s ey values 
            if i > 1
                solveMpcProblem(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:]',uCurr[i,:]')
            else
                solveMpcProblem(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:]',u_final') #in first round u_final is just a zero vector
            end
            tt[i]       = toq()
            cost[i,:]   = mpcSol.cost
            uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f]
            z_log[:,:,i] = mpcSol.z
            u_log[:,:,i] = mpcSol.u
            #!!todo sim model with xy cur
            #have Zcurr as states XY and simulate from there return XY values of states 
            @show zCurr_x[i+1,:]  = simModel_x(zCurr_x[i,:],uCurr[i,:],modelParams.dt,modelParams)
            println("Solving step $i of $(length(t)) - Status: $(mpcSol.solverStatus), Time: $(tt[i]) s")
              

           #T
            # if i%20 == 0 
              
            # #   plot(xs,ys) # tangent to current point  
            #    if j == 1
            #     figure(1)
            #     scatter(zCurr_x[i+1,1], zCurr_x[i+1,2], color = "red")
            # end
            # if j > 1
            #     figure(1)
            #     scatter(zCurr_x[i+1,1], zCurr_x[i+1,2], color = "green")
            # end
            # end 
            

            i = i + 1
            lapStatus.currentIt = i
        

        end
        if i >= length(t)
            println("took too long used whole buffersize")
        end
        
        i = i-1 #not necessary anymore as we jump out of the loop with the break point if all goes well so i is not augmented at end
        lapStatus.currentIt -= 1 # has finished lap already so we need to count both coutners one down as we counted them up after we have crossed the finish line
        z_final_x = zCurr_x[i,:]
        u_final = uCurr[i,:]
        println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i])) s")
        println("Finished Lap Nr. $j")

        # Save states in oldTraj:
        # --------------------------------
        #?? might be problematic with the s[i+1] issue above which we should aclculate at the end
        saveOldTraj(oldTraj,zCurr_s, zCurr_x,uCurr,lapStatus,buffersize,modelParams.dt)

          # Print results
        # --------------------------------

        # figure()
        # plot(zCurr[:,1],zCurr[:,2],"r",zCurr[:,1],zCurr[:,3],"g",zCurr[:,1],zCurr[:,4],"b")
        # grid(1)
        # legend(["eY","ePsi","v"])
        # title("States over s")

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
        plot(t,cost[:,1],"r",t,cost[:,2],"g",t,cost[:,3],"b",t,cost[:,4],"y",t,cost[:,5],"m",t,cost[:,6],"c")
        grid(1)
        title("Cost distribution")
        legend(["z","z_Term","z_Term_const","deriv","control","lane"])

  
        i_final = i
        #plot u
        if j == 1
            i_final_1st = i
            # ####T
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


        if j >= 2
            
            for i = 1:i_final

                figure(4)   
                clf()
                ax1=subplot(311)
                plot(oldTraj.oldTraj[1:i_final,1,1], oldTraj.oldTraj[1:i_final,4,1], color= "blue")
                plot(oldTraj.oldTraj[1:i_final_1st,1,2], oldTraj.oldTraj[1:i_final_1st,4,2], color= "yellow")
                legend(["v 2nd round","v 1st round"])

                ax1=subplot(312)
                plot(oldTraj.oldTraj[1:i_final,1,1], oldTraj.oldTraj[1:i_final,2,1], color= "blue")
                plot(oldTraj.oldTraj[1:i_final_1st,1,2], oldTraj.oldTraj[1:i_final_1st,2,2], color= "yellow")
                legend(["e_y 2nd round","e_y 1st round"])
                ax1=subplot(313)
                plot(oldTraj.oldTraj[1:i_final,1,1], oldTraj.oldTraj[1:i_final,3,1], color= "blue")
                plot(oldTraj.oldTraj[1:i_final_1st,1,2], oldTraj.oldTraj[1:i_final_1st,3,2], color= "yellow")
                legend(["e_psi 2nd round","e_psi 1st round"])
                subplot(311)
                plot(z_log[:,1,i], z_log[:,4,i] ,marker="o")
                subplot(312)
                plot(z_log[:,1,i], z_log[:,2,i] ,marker="o")
                 subplot(313)
                plot(z_log[:,1,i], z_log[:,3,i] ,marker="o")

                figure(5)
                clf()
                subplot(211) 
                plot(oldTraj.oldTraj[1:i_final-1,1,1], oldTraj.oldInput[1:i_final-1,1,1], color= "green")
                plot(z_log[1:end-1,1,i], u_log[:,1,i] ,marker="o")
                legend(["applied acceleration","predicted acceleration"])
                subplot(212)
                 plot(oldTraj.oldTraj[1:i_final-1,1,1], oldTraj.oldInput[1:i_final-1,2,1], color= "green")  
                plot(z_log[1:end-1,1,i], u_log[:,2,i] ,marker="o")  
                legend(["applied steering angle","predicted steering angle"]) 
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
        println("Press Enter for next round solving")
        println("Press c to cancel MPC")
        a = ' '
        a = readline()
        if a == "c\r\n" 
                break
        end
    end

end

end # end of module. we declare the main fail as a module to reduce warnings when including the file again