using JLD
using PyPlot
include("classes.jl")

function plots(j::Int64 = 1, interactive_plot::Int64 = 1)
    plot_costs = 1
    plot_states_over_t = 0
    plot_xy = 1
    plot_lambda = 1
    plot_states_over_s = 1
    plot_curvature_approx=0
    plot_inputs = 1
    interactive_plot_steps = 3
    file = "../LMPCdata/2016-11-21-00-04-Data.jld"
    close("all")


    ####load data from file
    Data = load(file)
    sStates_log = Data["sStates_log"]
    xStates_log = Data["xStates_log"]
    uAppl_log = Data["uAppl_log"]
    z_pred_log = Data["z_pred_log"]
    u_pred_log = Data["u_pred_log"]
    lambda_log = Data["lambda_log"]
    cost = Data["cost"]
    i_final = Data["i_final"]
    x_track = Data["x_track"]
    y_track = Data["y_track"]
    n_rounds = Data["n_rounds"]
    trackCoeff = Data["trackCoeff"]
    obstacle = Data["obstacle"]
    dt = Data["modelParams.dt"]
    mpcParams = Data["mpcParams"]
    buffersize = Data["buffersize"]
    curv_approx = Data["curv_approx"]
    #include("calculateObstacleXY.jl")
    include("helper/calculateObstacleXY.jl")
    #end of data loading

    #####create additional data for plotting
    xy_track  = [x_track; y_track]
    t   = collect(0:dt:(buffersize-1)*dt)
    xy_pred = zeros(mpcParams.N+1,2,length(t),n_rounds)

    println("Number of simulated rounds in data: $n_rounds")
    println("Load File located in: $file")

    for k = 1:n_rounds
        for i = 1:i_final[j]
            # caluclate the predicted XY postion of the car from the s-ey values
            xy_pred[:,:,i,k] = calculatePredictedXY(z_pred_log[:,:,:,k], mpcParams, trackCoeff, xy_track, convert(Int64,i),k)
            #calculate the obstacle postion from the s-ey values
            calculateObstacleXY!(obstacle, trackCoeff, xy_track,i,k) #this funciton computes values for row i
        end
    end
    ##this part is to calculate the tracks boundaries and plot them later
    convert(Array{Int64},i_final)
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

      # Print results
            # --------------------------------
    #########################################################
    #########################################################
    #### plot states and cost
    if plot_states_over_t == 1
        fig_1 = figure(1)
        fig_1[:canvas][:set_window_title]("States and Inputs over t")
        ax1=subplot(211)
        plot(t[1:i_final[j]],sStates_log[1:i_final[j],1,j],"y",t[1:i_final[j]],sStates_log[1:i_final[j],2,j],"r",t[1:i_final[j]],sStates_log[1:i_final[j],3,j],"g",t[1:i_final[j]],sStates_log[1:i_final[j],4,j],"b")
        grid(1)
        legend(["s","eY","ePsi","v"], bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        title("States")
        ax2=subplot(212,sharex=ax1)
        plot(t[1:i_final[j]-1],uAppl_log[1:i_final[j]-1,1,j],"r",t[1:i_final[j]-1],uAppl_log[1:i_final[j]-1,2,j],"g")
        grid(1)
        title("Control input")
        legend(["a","d_f"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        # ax3=subplot(313,sharex=ax1)
        # plot(t[1:i_final[j]-1],cost[1:i_final[j]-1,1,j],"r",t[1:i_final[j]-1],cost[1:i_final[j]-1,3,j],"b",t[1:i_final[j]-1],cost[1:i_final[j]-1,4,j],"y",t[1:i_final[j]-1],cost[1:i_final[j]-1,5,j],"m",t[1:i_final[j]-1],cost[1:i_final[j]-1,6,j],"c", t[1:i_final[j]-1], cost[1:i_final[j]-1,7,j])
        # grid(1)
        # title("Cost distribution")
        # legend(["z","z_Term_const","deriv","control","lane", "Obstacle"], bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    end

    ####plot of Cost components in seperated axis
    if plot_costs ==1
        plt_cost = figure(2)
        plt_cost[:canvas][:set_window_title]("Costs over s")
        ax9 = subplot(3,2,1)
        ax9[:plot](sStates_log[1:i_final[j]-1,1,j],cost[2,1:i_final[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("Terminal cost ") 

        act_s_plot = ax9[:plot]([sStates_log[i_final[j]-1,1,j],sStates_log[i_final[j]-1,1,j]],[0,30])  #vertical line to show actual s

        ax7= subplot(3,2,2,sharex=ax9)
        ax7[:plot](sStates_log[1:i_final[j]-1,1,j],cost[3,1:i_final[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("cost constraint")    
        ax8= subplot(3,2,3,sharex=ax9)
        ax8[:plot](sStates_log[1:i_final[j]-1,1,j],cost[7,1:i_final[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("cost Obstacle")     
        ax12= subplot(3,2,4,sharex=ax9)
        ax12[:plot](sStates_log[1:i_final[j]-1,1,j],cost[4,1:i_final[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("deriv cost")    
        ax13= subplot(3,2,5,sharex=ax9)
        ax13[:plot](sStates_log[1:i_final[j]-1,1,j],cost[6,1:i_final[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("lane Cost")     
        ax14= subplot(3,2,6,sharex=ax9)
        ax14[:plot](sStates_log[1:i_final[j]-1,1,j],cost[1,1:i_final[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("z Cost")  
    end

    #x-y Plot of the racetrack and obstacle and car
    if plot_xy == 1 
        f_xy_plot= figure(3)
        f_xy_plot[:canvas][:set_window_title]("Track and cars in XY plane")
        ax10 = subplot(1,1,1)
        ax10[:plot](x_track',y_track', linestyle="--", color = "yellow", linewidth = 0.5)#plot the racetrack
        car_plot = ax10[:plot](xStates_log[1:i_final[j],1,j], xStates_log[1:i_final[j],2,j], color = "blue")#always plot trajectory either to stay or delete in interactive mode
        if j >1#plot last tarjectory
            ax10[:plot](xStates_log[1:i_final[j-1],1,j-1], xStates_log[1:i_final[j-1],2,j-1], color = "green",linestyle=":")
        end
        
        if interactive_plot == 1
            ax10[:plot](xStates_log[1:i_final[j],1,j], xStates_log[1:i_final[j],2,j], color = "blue" ,linestyle=":")# plot trajectory less visible to overwrite it with interactiv simulation
            pred_plot = ax10[:plot](xy_pred[:,1,1,j],xy_pred[:,2,1,j],color = "yellow", marker="o") 
            obstacle_plot = ax10[:plot](obstacle.xy_vector[1,1,j], obstacle.xy_vector[1,2,j], color = "red",marker="o")
            y_obst_plot = ax10[:plot]([obstacle.axis_y_up[1,1,j],obstacle.axis_y_down[1,1,j]],[obstacle.axis_y_up[1,2,j],obstacle.axis_y_down[1,2,j]],color = "red")#plot the y semi axis
            s_obst_plot = ax10[:plot]([obstacle.axis_s_up[1,1,j],obstacle.axis_s_down[1,1,j]],[obstacle.axis_s_up[1,2,j],obstacle.axis_s_down[1,2,j]],color = "red")# plot the s semi axis
        else

            ax10[:plot](obstacle.xy_vector[1:i_final[j],1,j], obstacle.xy_vector[1:i_final[j],2,j], color = "red")
            y_obst_plot = ax10[:plot]([obstacle.axis_y_up[i_final[j],1,j],obstacle.axis_y_down[i_final[j],1,j]],[obstacle.axis_y_up[i_final[j],2,j],obstacle.axis_y_down[i_final[j],2,j]],color = "red")#plot the y semi axis
            s_obst_plot = ax10[:plot]([obstacle.axis_s_up[i_final[j],1,j],obstacle.axis_s_down[i_final[j],1,j]],[obstacle.axis_s_up[i_final[j],2,j],obstacle.axis_s_down[i_final[j],2,j]],color = "red")# plot the s semi axis
        end
        #plot the boundary lines
        ax10[:plot](boundary_up[1,:], boundary_up[2,:],color="green",linestyle="--")
        ax10[:plot](boundary_down[1,:], boundary_down[2,:],color="green",linestyle="--")
        gca()[:set_aspect]("equal", adjustable="box")
        grid() 
    end    
    ##########

#####old
    # #ax4[:set_xlim]([5.7,5.9])
    # #ax4[:set_ylim]([8.0,8.2])
    # figure(6)
        # plot(zCurr_s[1:i_final[j],1], zCurr_x[1:i_final[j],3])
        # xlabel("s in [m]")
        # ylabel("absolute angle psi in [rad]")
        # legend(["absolute angle psi plotted over s"]) 
######endold

    # plot the values of lambda over t
    if plot_lambda == 1 && j>1 
        f_lambda =figure(4)
        f_lambda[:canvas][:set_window_title]("Lambda values over t")
        ax11= subplot(1,1,1)
        plot(t[1:i_final[j]-1],lambda_log[1,1:i_final[j]-1,j])
        plot(t[1:i_final[j]-1],lambda_log[2,1:i_final[j]-1,j])
        plot(t[1:i_final[j]-1],lambda_log[3,1:i_final[j]-1,j])
        plot(t[1:i_final[j]-1],lambda_log[4,1:i_final[j]-1,j])
        plot(t[1:i_final[j]-1],lambda_log[5,1:i_final[j]-1,j])
        xlabel("t in [s]")
        ylabel("lambda")
        grid()
        ax11[:set_ylim]([-0.01,1.01])
        legend(["lambda1","lambda2","lambda3","lambda4","lambda5"])
    end


    ##plot the states over s and plot all old trajectories as well 
    if plot_states_over_s == 1
        f_states_over_s = figure(5)
        f_states_over_s[:canvas][:set_window_title]("States and safe set over s")   
        clf()
        ax4 = subplot(311)
        k = j
        while k>=1
            plot(sStates_log[1:i_final[k],1,k], sStates_log[1:i_final[k],4,k])
            k  = k-1
        end
        grid()
        xlabel("s in [m]")
        ylabel("v in [m/s]")
        legend(["v current","v 2nd", "3rd last","4th last ", "5th last"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        p1 = ax4[:plot](1,1)

        ax5 = subplot(312, sharex=ax4)
        k = j
        while k>=1
            plot(sStates_log[1:i_final[k],1,k], sStates_log[1:i_final[k],2,k])
            k  = k-1
        end
        grid()
        xlabel("s in [m]")
        ylabel("e_y in [m]")
        legend(["e_y current round","e_y 2nd", "3rd last","4th last ", "5th last"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        p2 = ax5[:plot](1,1)

        ax6 = subplot(313, sharex=ax4)
        hold(true)
        k = j
        while k>=1
            plot(sStates_log[1:i_final[k],1,k], sStates_log[1:i_final[k],3,k])
            k  = k-1
        end
        grid()
        xlabel("s in [m]")
        ylabel("e_psi in [rad]")
        legend(["e_psi current round","e_psi 2nd", "3rd last","4th last ", "5th last"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        p3 = ax6[:plot](1,1)
    end

    if plot_inputs == 1
        f_input = figure(7)
        f_input[:canvas][:set_window_title]("Inputs over s")  
        ax_Inp1 = subplot(211) 
        plot(sStates_log[1:i_final[j]-1,1,j], uAppl_log[1:i_final[j]-1,1,j], color= "green")
        xlabel("s in [m]")
        ylabel("acceleration in [m/s^2]")
        legend(["applied acceleration","predicted acceleration"])
        grid()
         ax_Inp2= subplot(212, sharex = ax_Inp1)
        plot(sStates_log[1:i_final[j]-1,1,j], uAppl_log[1:i_final[j]-1,2,j], color= "green")  
        xlabel("s in [m]")
        ylabel("steering angle in [rad]")
        legend(["applied steering angle","predicted steering angle"])
        grid() 
        inta_plot_a = ax_Inp1[:plot](1,1)
        inta_plot_df = ax_Inp2[:plot](1,1)
    end 

    #plot the estimated curvature for debugging
    if plot_curvature_approx==1
        f_curv_app = figure(6)
        f_curv_app[:canvas][:set_window_title]("Curvature approximation over s")  
        plot(sStates_log[1:i_final[j],1,j],curv_approx[1,1:i_final[j],j])
        grid()
    end
    
    if interactive_plot == 1
    for i = 1:interactive_plot_steps :i_final[j]# plot values at different times steps

        #plot predicted states over s
        if plot_states_over_s ==1
            p1[1][:remove]()
            p1 = ax4[:plot](z_pred_log[:,1,i,j], z_pred_log[:,4,i,j] ,marker="o") #plot predicted velocity

            p2[1][:remove]()
            p2 = ax5[:plot](z_pred_log[:,1,i,j], z_pred_log[:,2,i,j] ,marker="o") #plot predicted e_y
         
            p3[1][:remove]()
            p3 = ax6[:plot](z_pred_log[:,1,i,j], z_pred_log[:,3,i,j] ,marker="o")  #plot predicted e_psi
        end
       
       #x-y plot
       if plot_xy == 1 
            car_plot[1][:remove]()
            obstacle_plot[1][:remove]()
            y_obst_plot[1][:remove]()
            s_obst_plot[1][:remove]()
            pred_plot[1][:remove]()

            pred_plot = ax10[:plot](xy_pred[:,1,i,j],xy_pred[:,2,i,j], color = "yellow", marker="o") # plot predicted states
            ax10[:plot](xStates_log[1:i,1,j], xStates_log[1:i,2,j], color = "blue") # plot trajectory of this round for curent car
            car_plot = ax10[:plot](xStates_log[i,1,j], xStates_log[i,2,j], color = "blue", marker="o") #plot current position of car with a marker

            ax10[:plot](obstacle.xy_vector[1:i,1,j], obstacle.xy_vector[1:i,2,j], color = "red")
            obstacle_plot = ax10[:plot](obstacle.xy_vector[i,1,j], obstacle.xy_vector[i,2,j], color = "red", marker="o")  
            y_obst_plot = ax10[:plot]([obstacle.axis_y_up[i,1,j],obstacle.axis_y_down[i,1,j]],[obstacle.axis_y_up[i,2,j],obstacle.axis_y_down[i,2,j]],color = "red")#plot the y semi axis
            s_obst_plot = ax10[:plot]([obstacle.axis_s_up[i,1,j],obstacle.axis_s_down[i,1,j]],[obstacle.axis_s_up[i,2,j],obstacle.axis_s_down[i,2,j]],color = "red")# plot the s semi axis
        end
        ##plot a line in the cost function that always show the current s
        if plot_costs ==1
            act_s_plot[1][:remove]()
            act_s_plot= ax9[:plot]([sStates_log[i,1,j],sStates_log[i,1,j]],[-3,30])  
        end

        # #plot inputs and their predictions
        if plot_inputs ==1 
            inta_plot_a[1][:remove]()
            inta_plot_df[1][:remove]()
            inta_plot_a = ax_Inp1[:plot](z_pred_log[1:end-1,1,i,j], u_pred_log[:,1,i,j] ,marker="o")
            inta_plot_df = ax_Inp2[:plot](z_pred_log[1:end-1,1,i,j], u_pred_log[:,2,i,j] ,marker="o")   
        end

        println("Press Enter for next plot step")
        println("Press c to cancel plot")
        a = ' '
        a = readline()
        if a == "c\r\n" 
             break
        end
    end # end for loop interactive plots
    end #end if statement interactive plots
end