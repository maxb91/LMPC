using JLD
using PyPlot
include("classes.jl")

#function plots(j::Int64 = 1, interactive_plot::Int64 = 1)
    j = 1
    interactive_plot = 1

    plot_costs = 0
    plot_states_over_t = 0
    plot_xy = 1
    plot_lambda = 1
    plot_states_over_s = 1
    plot_curvature_approx=0
    plot_inputs = 0
    plot_eps = 0
    interactive_plot_steps = 4
    n_oldTrajPlots = 5
    file = "data/2017-01-16-12-14-Data.jld"
    close("all")

    ####load data from file
    Data = load(file)

   
    x_track = Data["x_track"]
    y_track = Data["y_track"]
    trackCoeff = Data["trackCoeff"]
    obstacle = Data["obstacle"]
    modelParams = Data["modelParams"]
    mpcParams = Data["mpcParams"]
    buffersize = Data["buffersize"]
    curv_approx = Data["curv_approx"]
    oldTraj     = Data["oldTraj"]
    mpcCoeff = Data["mpcCoeff"]
    include("calculateObstacleXY.jl")
    include("colorModule.jl")
    # include("helper/calculateObstacleXY.jl")
    # include("helper/colorModule.jl")
    #end of data loading

    #####create additional data for plotting
    dt = modelParams.dt
    xy_track  = [x_track; y_track]
    t   = collect(0:dt:(buffersize-1)*dt)
    xy_pred = zeros(mpcParams.N+1,2,length(t),oldTraj.n_oldTraj)

    println("Number of simulated rounds in data: $(oldTraj.n_oldTraj)")
    println("Load File located in: $file")

    for k = 1:oldTraj.n_oldTraj
        for i = 1:oldTraj.oldNIter[j]
            # caluclate the predicted XY postion of the car from the s-ey values
            xy_pred[:,:,i,k] = calculatePredictedXY(oldTraj.z_pred_sol[:,:,:,k], mpcParams, trackCoeff, xy_track, convert(Int64,i),k)
            #calculate the obstacle postion from the s-ey values
            calculateObstacleXY!(obstacle, trackCoeff, xy_track,i,k) #this funciton computes values for row i
        end
    end
    ################################
    ##this part is to calculate the tracks boundaries and plot them later
    ################################
    convert(Array{Int64},oldTraj.oldNIter)
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
    ################################
    ##calculate interpolated values to check for differences with trajectories
    # ################################
    # state_approx =zeros(buffersize,oldTraj.n_oldTraj,3)
    # s_vec = zeros(mpcCoeff.order+1)
    # for k = 2:oldTraj.n_oldTraj
    #     for i = 1:oldTraj.oldNIter[k]-1
    #         for l = 1:mpcCoeff.order+1
    #             s_vec[l] = oldTraj.z_pred_sol[11,1,i,k]^(mpcCoeff.order+1-l)
    #         end
    #         for m = 1:3
    #                state_approx[i,k,m] = dot(mpcCoeff.coeffConst[i,:,k-1,m],  s_vec) 
    #         end
    #     end
    # end
        

      # Print results
            # --------------------------------
    #########################################################
    #########################################################
    #### plot states and cost
    if plot_states_over_t == 1
        fig_1 = figure(1)
        fig_1[:canvas][:set_window_title]("States and Inputs over t")
        ax1=subplot(211)
        plot(t[1:oldTraj.oldNIter[j]],oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],"y",t[1:oldTraj.oldNIter[j]],oldTraj.oldTraj[1:oldTraj.oldNIter[j],2,j],"r",t[1:oldTraj.oldNIter[j]],oldTraj.oldTraj[1:oldTraj.oldNIter[j],3,j],"g",t[1:oldTraj.oldNIter[j]],oldTraj.oldTraj[1:oldTraj.oldNIter[j],4,j],"b")
        grid(1)
        legend(["s","eY","ePsi","v"], bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        title("States")
        ax2=subplot(212,sharex=ax1)
        plot(t[1:oldTraj.oldNIter[j]-1],oldTraj.oldInput[1:oldTraj.oldNIter[j]-1,1,j],"r",t[1:oldTraj.oldNIter[j]-1],oldTraj.oldInput[1:oldTraj.oldNIter[j]-1,2,j],"g")
        grid(1)
        title("Control input")
        legend(["a","d_f"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        # ax3=subplot(313,sharex=ax1)
        # plot(t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,1,j],"r",t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,3,j],"b",t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,4,j],"y",t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,5,j],"m",t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,6,j],"c", t[1:oldTraj.oldNIter[j]-1], oldTraj.costs[1:oldTraj.oldNIter[j]-1,7,j])
        # grid(1)
        # title("Cost distribution")
        # legend(["z","z_Term_const","deriv","control","lane", "Obstacle"], bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    end

    ####plot of Cost components in seperated axis
    if plot_costs ==1
        plt_cost = figure(2)
        plt_cost[:canvas][:set_window_title]("Costs over s")
        ax9 = subplot(3,2,1)
        ax9[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[2,1:oldTraj.oldNIter[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("Terminal cost ") 

        act_s_plot = ax9[:plot]([oldTraj.oldTraj[oldTraj.oldNIter[j]-1,1,j],oldTraj.oldTraj[oldTraj.oldNIter[j]-1,1,j]],[0,30])  #vertical line to show actual s

        ax7= subplot(3,2,2,sharex=ax9)
        ax7[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[3,1:oldTraj.oldNIter[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("cost constraint")    
        ax8= subplot(3,2,3,sharex=ax9)
        ax8[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[7,1:oldTraj.oldNIter[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("cost Obstacle")     
        ax12= subplot(3,2,4,sharex=ax9)
        ax12[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[4,1:oldTraj.oldNIter[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("deriv cost")    
        ax13= subplot(3,2,5,sharex=ax9)
        ax13[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[6,1:oldTraj.oldNIter[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("lane Cost")     
        ax14= subplot(3,2,6,sharex=ax9)
        ax14[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[1,1:oldTraj.oldNIter[j]-1,j])  
        grid()
        xlabel("s in [m]")
        ylabel("z Cost")  
    end

    #x-y Plot of the racetrack and obstacle and car
    if plot_xy == 1 
        f_xy_plot= figure(3)
        f_xy_plot[:canvas][:set_window_title]("Track and cars in XY plane")
        ax10 = subplot(1,1,1)
        ax10[:plot](x_track',y_track', linestyle="--", color = "#FF9900", linewidth = 0.5, label="_nolegend_")#plot the racetrack
        
        #plot older trajectory
        colorObjectXY= colorModule.ColorManager()
        for k= 1:n_oldTrajPlots
            colorXYold = colorModule.getColor(colorObjectXY)
            ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+k],1,j+k], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+k],2,j+k],linestyle="--", color = colorXYold, label= "$k old Traj")
        end
        
        if interactive_plot == 1
            car_plot = ax10[:plot](oldTraj.oldTrajXY[1,1,j], oldTraj.oldTrajXY[1,2,j], color = "black") # just dummy to use remove func later
            ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],2,j], color = "black" ,linestyle=":", label="current Traj")# plot trajectory less visible to overwrite it with interactiv simulation
            pred_plot = ax10[:plot](xy_pred[:,1,1,j],xy_pred[:,2,1,j],color = "yellow", marker="o") 
            obstacle_plot = ax10[:plot](obstacle.xy_vector[1,1,j], obstacle.xy_vector[1,2,j], color = "red",marker="o", label = "obstacle Traj")
            y_obst_plot = ax10[:plot]([obstacle.axis_y_up[1,1,j],obstacle.axis_y_down[1,1,j]],[obstacle.axis_y_up[1,2,j],obstacle.axis_y_down[1,2,j]],color = "red")#plot the y semi axis
            s_obst_plot = ax10[:plot]([obstacle.axis_s_up[1,1,j],obstacle.axis_s_down[1,1,j]],[obstacle.axis_s_up[1,2,j],obstacle.axis_s_down[1,2,j]],color = "red")# plot the s semi axis
        else
            car_plot = ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],2,j], color = "black", label="current Traj")
            ax10[:plot](obstacle.xy_vector[1:oldTraj.oldNIter[j],1,j], obstacle.xy_vector[1:oldTraj.oldNIter[j],2,j], color = "red")
            ax10[:plot](obstacle.xy_vector[oldTraj.oldNIter[j],1,j], obstacle.xy_vector[oldTraj.oldNIter[j],2,j], color = "red", marker = "o")
            y_obst_plot = ax10[:plot]([obstacle.axis_y_up[oldTraj.oldNIter[j],1,j],obstacle.axis_y_down[oldTraj.oldNIter[j],1,j]],[obstacle.axis_y_up[oldTraj.oldNIter[j],2,j],obstacle.axis_y_down[oldTraj.oldNIter[j],2,j]],color = "red")#plot the y semi axis
            s_obst_plot = ax10[:plot]([obstacle.axis_s_up[oldTraj.oldNIter[j],1,j],obstacle.axis_s_down[oldTraj.oldNIter[j],1,j]],[obstacle.axis_s_up[oldTraj.oldNIter[j],2,j],obstacle.axis_s_down[oldTraj.oldNIter[j],2,j]],color = "red")# plot the s semi axis
        end
        #plot the boundary lines
        ax10[:plot](boundary_up[1,:], boundary_up[2,:],color="green",linestyle=":")
        ax10[:plot](boundary_down[1,:], boundary_down[2,:],color="green",linestyle=":")
        for l=1:11
            ax10[:plot]([boundary_down[1,l*50+1],boundary_up[1,l*50+1]],[boundary_down[2,l*50+1],boundary_up[2,l*50+1]], color = "black", linestyle = ":", linewidth = 0.5)
            text(boundary_down[1,l*50+1]+1,boundary_down[2,l*50+1],"s = $(l*5)",fontsize=8)
        end
        gca()[:set_aspect]("equal", adjustable="box")
        # ax10[:set_ylim]([-5.1,5.1])
        grid() 
        legend(bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    end    

    # plot the values of lambda over t
    if plot_lambda == 1
        f_lambda =figure(4)
        f_lambda[:canvas][:set_window_title]("Lambda and ssOn values over t")
        ax11= subplot(2,2,1)
        colorObjectLambda = colorModule.ColorManager()
        for k=1: oldTraj.n_oldTraj//4
            k = convert(Int64,k)
            colorLambda = colorModule.getColor(colorObjectLambda)
            scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.lambda_sol[k,1:oldTraj.oldNIter[j]-1,j].*oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorLambda, marker = "x", label = string(L"\lambda","k"))
        end
        act_lambda_plot = ax11[:plot]([],[], label="_nolegend_")
        xlabel("s in [m]")
        ylabel("lambda")
        grid()
        #ax11[:set_ylim]([-0.1,1.1])
        legend()
        #legend([L"\lambda 1",L"\lambda 2",L"\lambda 3",L"\lambda 4",L"\lambda 5", L"\lambda 6",L"\lambda 7",L"\lambda 8",L"\lambda 9",L"\lambda 10"],bbox_to_anchor=(-0.201, 1), loc=2, borderaxespad=0.)
        
        axlambda2= subplot(2,2,2,sharex= ax11, sharey =ax11)
        for k=oldTraj.n_oldTraj//4+1: oldTraj.n_oldTraj//2
            k = convert(Int64,k)
            colorLambda = colorModule.getColor(colorObjectLambda)
            scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.lambda_sol[k,1:oldTraj.oldNIter[j]-1,j].*oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorLambda, marker = "x")
        end
        xlabel("s in [m]")
        ylabel("lambda")
        grid()
        legend([L"\lambda /4+1",L"\lambda /4+2",L"\lambda 3",L"\lambda 4",L"\lambda 5", L"\lambda 6"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)

        axlambda3= subplot(2,2,3,sharex= ax11, sharey =ax11)
        for k=oldTraj.n_oldTraj//2+1: oldTraj.n_oldTraj//4*3
            k = convert(Int64,k)
            colorLambda = colorModule.getColor(colorObjectLambda)
            scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.lambda_sol[k,1:oldTraj.oldNIter[j]-1,j].*oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorLambda, marker = "x")
        end
        xlabel("s in [m]")
        ylabel("lambda")
        grid()
        legend([L"\lambda /2+1",L"\lambda /2+2",L"\lambda 3",L"\lambda ",L"\lambda 5"],bbox_to_anchor=(-0.201, 1), loc=2, borderaxespad=0.)

        axlambda4= subplot(2,2,4,sharex= ax11, sharey =ax11)
        for k=oldTraj.n_oldTraj*3//4+1: oldTraj.n_oldTraj
            k = convert(Int64,k)
            colorLambda = colorModule.getColor(colorObjectLambda)
            scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.lambda_sol[k,1:oldTraj.oldNIter[j]-1,j].*oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorLambda, marker = "x")
        end
        xlabel("s in [m]")
        ylabel("lambda")
        grid()
        legend([L"\lambda *3/4+1",L"\lambda *3/4+2",L"\lambda 3",L"\lambda 4",L"\lambda 5"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        
    # plot values of ssOn over t    
        f_ssOn =figure(10)
        f_ssOn[:canvas][:set_window_title](" ssOn values over t")
        axssOn = subplot(2,2,1,sharex= ax11)
        colorObjectSafeSet= colorModule.ColorManager()
        for k=1: oldTraj.n_oldTraj//4
            k = convert(Int64,k)
            colorSafeSet= colorModule.getColor(colorObjectSafeSet)
            scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorSafeSet, marker = "x")
        end
        xlabel("s in [m]")
        ylabel("ssOn")
        grid()
        legend(["ssOn1","ssOn2","ssOn3","ssOn4","ssOn5", "ssOn6","ssOn7","ssOn8","ssOn9","ssOn10"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        
        axssOn2 = subplot(2,2,2,sharex= ax11)
        for k=oldTraj.n_oldTraj//4+1: oldTraj.n_oldTraj//2
            k = convert(Int64,k)
            colorSafeSet= colorModule.getColor(colorObjectSafeSet)
            scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j],color = colorSafeSet, marker = "x")
        end
        xlabel("s in [m]")
        ylabel("ssOn")
        grid()
        legend(["ssOn/4+1","ssOn/4+2","ssOn3","ssOn4","ssOn5", "ssOn6","ssOn7","ssOn8","ssOn9","ssOn10"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        axssOn3 = subplot(2,2,3,sharex= ax11)
        for k=oldTraj.n_oldTraj//2+1: oldTraj.n_oldTraj*3//4
            k = convert(Int64,k)
            colorSafeSet= colorModule.getColor(colorObjectSafeSet)
            scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j],color = colorSafeSet, marker = "x")
        end
        xlabel("s in [m]")
        ylabel("ssOn")
        grid()
        legend(["ssOn/2+1","ssOn/2+2","ssOn3","ssOn4","ssOn5", "ssOn6","ssOn7","ssOn8","ssOn9","ssOn10"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        axssOn4 = subplot(2,2,4,sharex= ax11)
        for k=oldTraj.n_oldTraj*3//4+1: oldTraj.n_oldTraj
            k = convert(Int64,k)
            colorSafeSet= colorModule.getColor(colorObjectSafeSet)
            scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j],color = colorSafeSet, marker = "x")
        end
        xlabel("s in [m]")
        ylabel("ssOn")
        grid()
        legend(["ssOn*3/4+1","ssOn/2+2","ssOn3","ssOn4","ssOn5", "ssOn6","ssOn7","ssOn8","ssOn9","ssOn10"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    end


    ##plot the states over s and plot all old trajectories as well 
    if plot_states_over_s == 1
        f_states_over_s = figure(5)
        f_states_over_s[:canvas][:set_window_title]("States and safe set over s")   
        clf()
        ax4 = subplot(311)
        plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTraj[1:oldTraj.oldNIter[j],4,j], color = "black")
        k = j+1
        colorObjectV= colorModule.ColorManager()
        while k<=j+n_oldTrajPlots
            colorV = colorModule.getColor(colorObjectV)
            plot(oldTraj.oldTraj[1:oldTraj.oldNIter[k],1,k], oldTraj.oldTraj[1:oldTraj.oldNIter[k],4,k], color = colorV)
            #plot(oldTraj.z_pred_sol[11,1,1:oldTraj.oldNIter[k]-1,k],state_approx[1:oldTraj.oldNIter[k]-1,k,3],color = colorV,label="_nolegend_", linestyle = "--")
            k  = k+1
        end
        grid()
        xlabel("s in [m]")
        ylabel("v in [m/s]")
        legend(["v current","v 2nd", "3rd last","4th last ", "5th last"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        p1 = ax4[:plot](1,1)
        ax4[:set_ylim](0.0,2.1)

        ax5 = subplot(312, sharex=ax4)
        plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTraj[1:oldTraj.oldNIter[j],2,j], color = "black")
        k = j+1
        colorObject_eY= colorModule.ColorManager()
        while k<=j+n_oldTrajPlots
            color_eY = colorModule.getColor(colorObject_eY)
            plot(oldTraj.oldTraj[1:oldTraj.oldNIter[k],1,k], oldTraj.oldTraj[1:oldTraj.oldNIter[k],2,k], color = color_eY)
            #plot(oldTraj.z_pred_sol[11,1,1:oldTraj.oldNIter[k]-1,k],state_approx[1:oldTraj.oldNIter[k]-1,k,1],color = color_eY,label="_nolegend_", linestyle = "--")
            k  = k+1
        end
        grid()
        xlabel("s in [m]")
        ylabel("e_y in [m]")
        legend(["e_y current round","e_y 2nd", "3rd last","4th last ", "5th last"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        p2 = ax5[:plot](1,1)

        ax6 = subplot(313, sharex=ax4)
        hold(true)
        plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTraj[1:oldTraj.oldNIter[j],3,j],color = "black")
        k = j+1
        colorObject_ePsi= colorModule.ColorManager()
        while k<=j+n_oldTrajPlots
            color_ePsi = colorModule.getColor(colorObject_ePsi)
            plot(oldTraj.oldTraj[1:oldTraj.oldNIter[k],1,k], oldTraj.oldTraj[1:oldTraj.oldNIter[k],3,k],color = color_ePsi)
            #plot(oldTraj.z_pred_sol[11,1,1:oldTraj.oldNIter[k]-1,k],state_approx[1:oldTraj.oldNIter[k]-1,k,2],color = color_ePsi,label="_nolegend_", linestyle = "--")
            k  = k+1
        end
        grid()
        xlabel("s in [m]")
        ylabel("e_psi in [rad]")
        legend(["e_psi current round","e_psi 2nd", "3rd last","4th last ", "5th last"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        p3 = ax6[:plot](1,1)
    end

    if plot_inputs == 1
        f_input = figure(6)
        f_input[:canvas][:set_window_title]("Inputs over s")  
        ax_Inp1 = subplot(211) 
        plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j], oldTraj.oldInput[1:oldTraj.oldNIter[j]-1,1,j], color= "green")
        xlabel("s in [m]")
        ylabel("acceleration in [m/s^2]")
        legend(["applied acceleration","predicted acceleration"])
        grid()
         ax_Inp2= subplot(212, sharex = ax_Inp1)
        plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j], oldTraj.oldInput[1:oldTraj.oldNIter[j]-1,2,j], color= "green")  
        xlabel("s in [m]")
        ylabel("steering angle in [rad]")
        legend(["applied steering angle","predicted steering angle"])
        grid() 
        inta_plot_a = ax_Inp1[:plot](1,1)
        inta_plot_df = ax_Inp2[:plot](1,1)
    end 

    #plot the estimated curvature for debugging
    if plot_curvature_approx==1
        f_curv_app = figure(7)
        f_curv_app[:canvas][:set_window_title]("Curvature approximation over s")  
        #plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],curv_approx[1,1:oldTraj.oldNIter[j],j])
        plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],oldTraj.curvature[1:oldTraj.oldNIter[j],j])
        grid()
    end

    if plot_eps ==1
        f_eps = figure(8)
        f_eps[:canvas][:set_window_title]("values of Eps over s") 
    
        plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],oldTraj.eps[1,1:oldTraj.oldNIter[j],j])
        plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],oldTraj.eps[2,1:oldTraj.oldNIter[j],j])
        xlabel("s in [m]")
        ylabel("value of epsilons")
        legend(["epsilon left boundary","epsilon right boundary", "epsilon obstacle"])
        grid()
    end
    
    if interactive_plot == 1
    for i = 1:interactive_plot_steps :oldTraj.oldNIter[j]# plot values at different times steps

        #plot predicted states over s
        if plot_states_over_s ==1
            p1[1][:remove]()
            p1 = ax4[:plot](oldTraj.z_pred_sol[:,1,i,j], oldTraj.z_pred_sol[:,4,i,j] ,marker="o", color = "yellow") #plot predicted velocity

            p2[1][:remove]()
            p2 = ax5[:plot](oldTraj.z_pred_sol[:,1,i,j], oldTraj.z_pred_sol[:,2,i,j] ,marker="o", color = "yellow") #plot predicted e_y
         
            p3[1][:remove]()
            p3 = ax6[:plot](oldTraj.z_pred_sol[:,1,i,j], oldTraj.z_pred_sol[:,3,i,j] ,marker="o", color = "yellow")  #plot predicted e_psi
        end
       
       #x-y plot
       if plot_xy == 1 
            car_plot[1][:remove]()
            obstacle_plot[1][:remove]()
            y_obst_plot[1][:remove]()
            s_obst_plot[1][:remove]()
            pred_plot[1][:remove]()

            pred_plot = ax10[:plot](xy_pred[:,1,i,j],xy_pred[:,2,i,j], color = "yellow", marker="o") # plot predicted states
            ax10[:plot](oldTraj.oldTrajXY[1:i,1,j], oldTraj.oldTrajXY[1:i,2,j], color = "black") # plot trajectory of this round for curent car
            car_plot = ax10[:plot](oldTraj.oldTrajXY[i,1,j], oldTraj.oldTrajXY[i,2,j], color = "black", marker="o") #plot current position of car with a marker

            ax10[:plot](obstacle.xy_vector[1:i,1,j], obstacle.xy_vector[1:i,2,j], color = "red", linestyle= ":")
            obstacle_plot = ax10[:plot](obstacle.xy_vector[i,1,j], obstacle.xy_vector[i,2,j], color = "red", marker="o")  
            y_obst_plot = ax10[:plot]([obstacle.axis_y_up[i,1,j],obstacle.axis_y_down[i,1,j]],[obstacle.axis_y_up[i,2,j],obstacle.axis_y_down[i,2,j]],color = "red")#plot the y semi axis
            s_obst_plot = ax10[:plot]([obstacle.axis_s_up[i,1,j],obstacle.axis_s_down[i,1,j]],[obstacle.axis_s_up[i,2,j],obstacle.axis_s_down[i,2,j]],color = "red")# plot the s semi axis
        end
        ##plot a line in the cost function that always show the current s
        if plot_costs ==1
            act_s_plot[1][:remove]()
            act_s_plot= ax9[:plot]([oldTraj.oldTraj[i,1,j],oldTraj.oldTraj[i,1,j]],[-3,oldTraj.oldNIter[j]], color ="black")  
        end
        if plot_lambda ==1
            act_lambda_plot[1][:remove]()
            act_lambda_plot= ax11[:plot]([t[i],t[i]],[0,1],label="_nolegend_", color ="black")  
        end

        # #plot inputs and their predictions
        if plot_inputs ==1 
            inta_plot_a[1][:remove]()
            inta_plot_df[1][:remove]()
            inta_plot_a = ax_Inp1[:plot](oldTraj.z_pred_sol[1:end-1,1,i,j], oldTraj.u_pred_sol[:,1,i,j] ,marker="o", color = "yellow")
            inta_plot_df = ax_Inp2[:plot](oldTraj.z_pred_sol[1:end-1,1,i,j], oldTraj.u_pred_sol[:,2,i,j] ,marker="o", color = "yellow")   
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
    #return mpcParams, oldTraj
#end