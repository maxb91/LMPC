using JLD
using PyPlot
using PyCall
@pyimport matplotlib.patches as patches


include("classes.jl")
include("plot_functions.jl")

#function plots(j::Int64 = 1, interactive_plot::Int64 = 1)
    newest2plot =1
    n_plot_rounds = 0

    interactive_plot = 1

    plot_costs = 0
    plot_states_over_t = 0
    plot_xy = 0
    plot_lambda = 0
    plot_states_over_s = 0
    plot_curvature_approx=1
    plot_inputs = 0
    plot_eps = 0
    plot_copied = 0
    interactive_plot_steps = 5
    n_oldTrajPlots = 2
    file = "data/2017-02-03-11-23-Data.jld"
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
    oldTraj     = Data["oldTraj"]
    mpcCoeff = Data["mpcCoeff"]
    include("calculateObstacleXY.jl")
    include("colorModule.jl")
    # include("helper/calculateObstacleXY.jl")
    # include("helper/colorModule.jl")
    #end of data loading
    println("Number of simulated rounds in data: $(oldTraj.n_oldTraj)")
    println("Load File located in: $file")

    #####create additional data for plotting
    ds = trackCoeff.ds
    dt = modelParams.dt
    xy_track  = [x_track; y_track]
    t   = collect(0:dt:(buffersize-1)*dt)




    colordefs=["#3b44ba",
            "#c8ce24",
            "#a756de",
            "#0ac753",
            "#ff53c4",
            "#01c3ba",
            "#f5218a",
            "#4cd5ff",
            "#f74b3a",
            "#016ed6",
            "#de8500",
            "#8a1e95",
            "#6c4d08",
            "#aca0ff",
            "#a2172f",
            "#abaae1",
            "#ffb38d",
            "#83317e",
            "#9c6f4b",
            "#ff91d1"]

    ################################
    ##this part is to calculate the tracks boundaries and plot them later
    ################################
    trackL = size(xy_track,2)
    boundary_up = zeros(2,trackL)
    boundary_down = zeros(2,trackL)
    trackCoeff.width = trackCoeff.width+0.2
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
    obstOrientation = zeros(size(obstacle.xy_vector)[1],size(obstacle.xy_vector)[3],size(obstacle.xy_vector)[4])
    if plot_xy == 1
        f_xy_plot= figure(3)
        f_xy_plot[:canvas][:set_window_title]("Track and cars in XY plane")
        ax10 = subplot(1,1,1)
        if interactive_plot == 1
            cartraj_plot = ax10[:plot](1,1)  #just for initialization
            #obsttraj_plot1 = ax10[:plot](1,1)  #just for initialization
            car_plot = ax10[:plot](oldTraj.oldTrajXY[1,1,1], oldTraj.oldTrajXY[1,2,1], color = "black") # just dummy to use remove func later
            pred_plot = ax10[:plot](oldTraj.oldTrajXY[1,1,1,1],oldTraj.oldTrajXY[1,2,1,1],color = "yellow", marker="o")

       
            if obstacle.n_obstacle >=1
                #obsttraj_plot1 = ax10[:plot](1,1)  #just for initialization
                obstacle_plot1 = ax10[:plot](obstacle.xy_vector[1,1,1,1], obstacle.xy_vector[1,2,1,1], color = "red",marker="o", label = "obstacle Traj")
                y_obst_plot1   = ax10[:plot]([obstacle.axis_y_up[1,1,1,1],obstacle.axis_y_down[1,1,1,1]],[obstacle.axis_y_up[1,2,1,1],obstacle.axis_y_down[1,2,1,1]],color = "red")#plot the y semi axis
                s_obst_plot1   = ax10[:plot]([obstacle.axis_s_up[1,1,1,1],obstacle.axis_s_down[1,1,1,1]],[obstacle.axis_s_up[1,2,1,1],obstacle.axis_s_down[1,2,1,1]],color = "red")# plot the s semi axis

                obst1 = patches.Ellipse([obstacle.xy_vector[1,1,1,1],obstacle.xy_vector[1,2,1,1]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[1,1,1]*180/pi)
                plt_obst1 = ax10[:add_patch](obst1)
            end
            if obstacle.n_obstacle >=2
                #obsttraj_plot2 = ax10[:plot](1,1)  #just for initialization
                obstacle_plot2 = ax10[:plot](obstacle.xy_vector[1,1,1,2], obstacle.xy_vector[1,2,1,2], color = "red",marker="o", label = "obstacle Traj")
                y_obst_plot2   = ax10[:plot]([obstacle.axis_y_up[1,1,1,2],obstacle.axis_y_down[1,1,1,2]],[obstacle.axis_y_up[1,2,1,2],obstacle.axis_y_down[1,2,1,2]],color = "red")#plot the y semi axis
                s_obst_plot2   = ax10[:plot]([obstacle.axis_s_up[1,1,1,2],obstacle.axis_s_down[1,1,1,2]],[obstacle.axis_s_up[1,2,1,2],obstacle.axis_s_down[1,2,1,2]],color = "red")# plot the s semi axis

                obst2 = patches.Ellipse([obstacle.xy_vector[1,1,1,2],obstacle.xy_vector[1,2,1,2]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[1,1,2]*180/pi)
                plt_obst2 = ax10[:add_patch](obst2)
            end
            if obstacle.n_obstacle >=3
                #obsttraj_plot3 = ax10[:plot](1,1)  #just for initialization
                obstacle_plot3 = ax10[:plot](obstacle.xy_vector[1,1,1,3], obstacle.xy_vector[1,2,1,3], color = "red",marker="o", label = "obstacle Traj")
                y_obst_plot3   = ax10[:plot]([obstacle.axis_y_up[1,1,1,3],obstacle.axis_y_down[1,1,1,3]],[obstacle.axis_y_up[1,2,1,3],obstacle.axis_y_down[1,2,1,3]],color = "red")#plot the y semi axis
                s_obst_plot3   = ax10[:plot]([obstacle.axis_s_up[1,1,1,3],obstacle.axis_s_down[1,1,1,3]],[obstacle.axis_s_up[1,2,1,3],obstacle.axis_s_down[1,2,1,3]],color = "red")# plot the s semi axis

                obst3 = patches.Ellipse([obstacle.xy_vector[1,1,1,3],obstacle.xy_vector[1,2,1,3]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[1,1,3]*180/pi)
                plt_obst3 = ax10[:add_patch](obst3)
            end
            if obstacle.n_obstacle >=4
                #obsttraj_plot4 = ax10[:plot](1,1)  #just for initialization
                obstacle_plot4 = ax10[:plot](obstacle.xy_vector[1,1,1,4], obstacle.xy_vector[1,2,1,4], color = "red",marker="o", label = "obstacle Traj")
                y_obst_plot4   = ax10[:plot]([obstacle.axis_y_up[1,1,1,4],obstacle.axis_y_down[1,1,1,4]],[obstacle.axis_y_up[1,2,1,4],obstacle.axis_y_down[1,2,1,4]],color = "red")#plot the y semi axis
                s_obst_plot4   = ax10[:plot]([obstacle.axis_s_up[1,1,1,4],obstacle.axis_s_down[1,1,1,4]],[obstacle.axis_s_up[1,2,1,4],obstacle.axis_s_down[1,2,1,4]],color = "red")# plot the s semi axis

                obst4 = patches.Ellipse([obstacle.xy_vector[1,1,1,4],obstacle.xy_vector[1,2,1,4]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[1,1,4]*180/pi)
                plt_obst4 = ax10[:add_patch](obst4)
            end
            if obstacle.n_obstacle >=5
                #obsttraj_plot5 = ax10[:plot](1,1)  #just for initialization
                obstacle_plot5 = ax10[:plot](obstacle.xy_vector[1,1,1,5], obstacle.xy_vector[1,2,1,5], color = "red",marker="o", label = "obstacle Traj")
                y_obst_plot5   = ax10[:plot]([obstacle.axis_y_up[1,1,1,5],obstacle.axis_y_down[1,1,1,5]],[obstacle.axis_y_up[1,2,1,5],obstacle.axis_y_down[1,2,1,5]],color = "red")#plot the y semi axis
                s_obst_plot5   = ax10[:plot]([obstacle.axis_s_up[1,1,1,5],obstacle.axis_s_down[1,1,1,5]],[obstacle.axis_s_up[1,2,1,5],obstacle.axis_s_down[1,2,1,5]],color = "red")# plot the s semi axis

                obst5 = patches.Ellipse([obstacle.xy_vector[1,1,1,5],obstacle.xy_vector[1,2,1,5]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[1,1,5]*180/pi)
                plt_obst5 = ax10[:add_patch](obst5)
            end
            if obstacle.n_obstacle >=6
                #obsttraj_plot6 = ax10[:plot](1,1)  #just for initialization
                obstacle_plot6 = ax10[:plot](obstacle.xy_vector[1,1,1,6], obstacle.xy_vector[1,2,1,6], color = "red",marker="o", label = "obstacle Traj")
                y_obst_plot6   = ax10[:plot]([obstacle.axis_y_up[1,1,1,6],obstacle.axis_y_down[1,1,1,6]],[obstacle.axis_y_up[1,2,1,6],obstacle.axis_y_down[1,2,1,6]],color = "red")#plot the y semi axis
                s_obst_plot6   = ax10[:plot]([obstacle.axis_s_up[1,1,1,6],obstacle.axis_s_down[1,1,1,6]],[obstacle.axis_s_up[1,2,1,6],obstacle.axis_s_down[1,2,1,6]],color = "red")# plot the s semi axis
                
                obst6 = patches.Ellipse([obstacle.xy_vector[1,1,1,6],obstacle.xy_vector[1,2,1,6]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[1,1,6]*180/pi)
                plt_obst6 = ax10[:add_patch](obst6)
            end
            if obstacle.n_obstacle >=7
                #obsttraj_plot6 = ax10[:plot](1,1)  #just for initialization
                obstacle_plot7 = ax10[:plot](obstacle.xy_vector[1,1,1,7], obstacle.xy_vector[1,2,1,7], color = "red",marker="o", label = "obstacle Traj")
                y_obst_plot7   = ax10[:plot]([obstacle.axis_y_up[1,1,1,7],obstacle.axis_y_down[1,1,1,7]],[obstacle.axis_y_up[1,2,1,7],obstacle.axis_y_down[1,2,1,7]],color = "red")#plot the y semi axis
                s_obst_plot7   = ax10[:plot]([obstacle.axis_s_up[1,1,1,7],obstacle.axis_s_down[1,1,1,7]],[obstacle.axis_s_up[1,2,1,7],obstacle.axis_s_down[1,2,1,7]],color = "red")# plot the s semi axis

                obst7 = patches.Ellipse([obstacle.xy_vector[1,1,1,7],obstacle.xy_vector[1,2,1,7]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[1,1,7]*180/pi)
                plt_obst7 = ax10[:add_patch](obst7)
            end
            ax10[:grid]() 
        end
    end
    stop_interactive = false
    for j=newest2plot+n_plot_rounds:-1:newest2plot

        ################################
        ##calculate predcted position of the car in XY plane
        # ################################
        xy_pred = zeros(mpcParams.N+1,2,size(t)[1],oldTraj.n_oldTraj)
        
        for k = 1:oldTraj.n_oldTraj
            for i = 1:oldTraj.oldNIter[j]
                # caluclate the predicted XY postion of the car from the s-ey values
                xy_pred[:,:,i,k] = calculatePredictedXY(oldTraj.z_pred_sol[:,:,:,k], mpcParams, trackCoeff, xy_track, convert(Int64,i),k)
                #calculate the obstacle postion from the s-ey values
                obstOrientation = calculateObstacleXY!(obstacle, trackCoeff, xy_track,i,k,obstOrientation ) #this funciton computes values for row i
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
           plotfct_states_over_t(oldTraj, t, j)
        end

        ####plot of Cost components in seperated axis
        if plot_costs ==1
            ax9 = plotfct_costs(oldTraj, j)
            act_s_plot = ax9[:plot]([oldTraj.oldTraj[oldTraj.oldNIter[j]-1,1,j],oldTraj.oldTraj[oldTraj.oldNIter[j]-1,1,j]],[0,30])  #vertical line to show actual s
        end

        #x-y Plot of the racetrack and obstacle and car
        if plot_xy == 1 
            # f_xy_plot= figure(3)
            # f_xy_plot[:canvas][:set_window_title]("Track and cars in XY plane")
            # ax10 = subplot(1,1,1)
            # ax10[:clear]()
            ax10[:plot](x_track',y_track', linestyle="--", color = "#FF9900", linewidth = 0.5, label="_nolegend_")#plot the racetrack
            
            #plot older trajectory
            colorObjectXY= colorModule.ColorManager()
            for k= 1:n_oldTrajPlots
                colorXYold = colorModule.getColor(colorObjectXY)
                ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+k],1,j+k], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+k],2,j+k],linestyle="--", color = colorXYold, label= "$k old Traj")
            end
            
            if interactive_plot == 1

                ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],2,j], color = "black" ,linestyle=":", label="current Traj")# plot trajectory less visible to overwrite it with interactiv simulation
 

            else
                car_plot = ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],2,j], color = "black", label="current Traj")
                ax10[:plot](obstacle.xy_vector[1:oldTraj.oldNIter[j],1,j,1], obstacle.xy_vector[1:oldTraj.oldNIter[j],2,j,1], color = "red")
                ax10[:plot](obstacle.xy_vector[oldTraj.oldNIter[j],1,j,1], obstacle.xy_vector[oldTraj.oldNIter[j],2,j,1], color = "red", marker = "o")
                y_obst_plot = ax10[:plot]([obstacle.axis_y_up[oldTraj.oldNIter[j],1,j,1],obstacle.axis_y_down[oldTraj.oldNIter[j],1,j,1]],[obstacle.axis_y_up[oldTraj.oldNIter[j],2,j,1],obstacle.axis_y_down[oldTraj.oldNIter[j],2,j,1]],color = "red")#plot the y semi axis
                s_obst_plot = ax10[:plot]([obstacle.axis_s_up[oldTraj.oldNIter[j],1,j,1],obstacle.axis_s_down[oldTraj.oldNIter[j],1,j,1]],[obstacle.axis_s_up[oldTraj.oldNIter[j],2,j,1],obstacle.axis_s_down[oldTraj.oldNIter[j],2,j,1]],color = "red")# plot the s semi axis
            end
            #plot the boundary lines
            ax10[:plot](boundary_up[1,:], boundary_up[2,:],color="green",linestyle=":")
            ax10[:plot](boundary_down[1,:], boundary_down[2,:],color="green",linestyle=":")
            
            for l=1:convert(Int64,trunc(trackL/51))
                ax10[:plot]([boundary_down[1,l*50+1],boundary_up[1,l*50+1]],[boundary_down[2,l*50+1],boundary_up[2,l*50+1]], color = "black", linestyle = ":", linewidth = 0.5)
                ax10[:text](boundary_down[1,l*50+1]+1,boundary_down[2,l*50+1],"s = $(convert(Int64,l*50*ds))",fontsize=8)
            end
            ax10[:set_aspect]("equal", adjustable="box")
            # ax10[:set_ylim]([-5.1,5.1])

            ax10[:legend](bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
        end    

        # plot the values of lambda over t
        if plot_lambda == 1

                ax11 = plotfct_lambda(oldTraj,j)
                plotfct_ssOn(oldTraj,j, ax11)
                act_lambda_plot = ax11[:plot]([],[], label="_nolegend_")
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
            plotfct_curvature(oldTraj,j)
        end

        if plot_eps ==1
            plotfct_epsilon(oldTraj, j)
        end

        if plot_copied == 1
            plotfct_copied(oldTraj,j)
        end

        ################################################################################################
        #                                                                                              #
        ################################################################################################
        if interactive_plot == 1
        for i = 1:interactive_plot_steps :oldTraj.oldNIter[j]-1# plot values at different times steps

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
                pred_plot[1][:remove]()
                cartraj_plot[1][:remove]()
                

                pred_plot = ax10[:plot](xy_pred[:,1,i,j],xy_pred[:,2,i,j], color = "yellow", marker="o") # plot predicted states
                cartraj_plot = ax10[:plot](oldTraj.oldTrajXY[1:i,1,j], oldTraj.oldTrajXY[1:i,2,j], color = "black") # plot trajectory of this round for curent car
                car_plot = ax10[:plot](oldTraj.oldTrajXY[i,1,j], oldTraj.oldTrajXY[i,2,j], color = "black", marker="o") #plot current position of car with a marker

                if obstacle.n_obstacle >=1
                    # obsttraj_plot1[1][:remove]()
                    obstacle_plot1[1][:remove]()
                    y_obst_plot1[1][:remove]()
                    s_obst_plot1[1][:remove]()
                    # obsttraj_plot1 = ax10[:plot](obstacle.xy_vector[1:i,1,j,1], obstacle.xy_vector[1:i,2,j,1], color = "red", linestyle= ":")
                    obstacle_plot1 = ax10[:plot](obstacle.xy_vector[i,1,j,1], obstacle.xy_vector[i,2,j,1], color = "red", marker="o")  
                    y_obst_plot1 = ax10[:plot]([obstacle.axis_y_up[i,1,j,1],obstacle.axis_y_down[i,1,j,1]],[obstacle.axis_y_up[i,2,j,1],obstacle.axis_y_down[i,2,j,1]],color = "red")#plot the y semi axis
                    s_obst_plot1 = ax10[:plot]([obstacle.axis_s_up[i,1,j,1],obstacle.axis_s_down[i,1,j,1]],[obstacle.axis_s_up[i,2,j,1],obstacle.axis_s_down[i,2,j,1]],color = "red")# plot the s semi axis

                    plt_obst1[:remove]()
                    obst1 = patches.Ellipse([obstacle.xy_vector[i,1,j,1],obstacle.xy_vector[i,2,j,1]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[i,j,1]*180/pi)
                    plt_obst1 = ax10[:add_patch](obst1)

                end
                if obstacle.n_obstacle >=2
                    # obsttraj_plot2[1][:remove]()
                    obstacle_plot2[1][:remove]()
                    y_obst_plot2[1][:remove]()
                    s_obst_plot2[1][:remove]()
                    # obsttraj_plot2 = ax10[:plot](obstacle.xy_vector[1:i,1,j,2], obstacle.xy_vector[1:i,2,j,2], color = "red", linestyle= ":")
                    obstacle_plot2 = ax10[:plot](obstacle.xy_vector[i,1,j,2], obstacle.xy_vector[i,2,j,2], color = "red", marker="o")  
                    y_obst_plot2 = ax10[:plot]([obstacle.axis_y_up[i,1,j,2],obstacle.axis_y_down[i,1,j,2]],[obstacle.axis_y_up[i,2,j,2],obstacle.axis_y_down[i,2,j,2]],color = "red")#plot the y semi axis
                    s_obst_plot2 = ax10[:plot]([obstacle.axis_s_up[i,1,j,2],obstacle.axis_s_down[i,1,j,2]],[obstacle.axis_s_up[i,2,j,2],obstacle.axis_s_down[i,2,j,2]],color = "red")# plot the s semi axis

                    plt_obst2[:remove]()
                    obst2 = patches.Ellipse([obstacle.xy_vector[i,1,j,2],obstacle.xy_vector[i,2,j,2]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[i,j,2]*180/pi)
                    plt_obst2 = ax10[:add_patch](obst2)
                end
                if obstacle.n_obstacle >=3
                    # obsttraj_plot3[1][:remove]()
                    obstacle_plot3[1][:remove]()
                    y_obst_plot3[1][:remove]()
                    s_obst_plot3[1][:remove]()
                    # obsttraj_plot2 = ax10[:plot](obstacle.xy_vector[1:i,1,j,3], obstacle.xy_vector[1:i,2,j,3], color = "red", linestyle= ":")
                    obstacle_plot3 = ax10[:plot](obstacle.xy_vector[i,1,j,3], obstacle.xy_vector[i,2,j,3], color = "red", marker="o")  
                    y_obst_plot3 = ax10[:plot]([obstacle.axis_y_up[i,1,j,3],obstacle.axis_y_down[i,1,j,3]],[obstacle.axis_y_up[i,2,j,3],obstacle.axis_y_down[i,2,j,3]],color = "red")#plot the y semi axis
                    s_obst_plot3 = ax10[:plot]([obstacle.axis_s_up[i,1,j,3],obstacle.axis_s_down[i,1,j,3]],[obstacle.axis_s_up[i,2,j,3],obstacle.axis_s_down[i,2,j,3]],color = "red")# plot the s semi axis

                    plt_obst3[:remove]()
                    obst3 = patches.Ellipse([obstacle.xy_vector[i,1,j,3],obstacle.xy_vector[i,2,j,3]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[i,j,3]*180/pi)
                    plt_obst3 = ax10[:add_patch](obst3)
                end
                if obstacle.n_obstacle >=4
                    # obsttraj_plot4[1][:remove]()
                    obstacle_plot4[1][:remove]()
                    y_obst_plot4[1][:remove]()
                    s_obst_plot4[1][:remove]()
                    # obsttraj_plot4 = ax10[:plot](obstacle.xy_vector[1:i,1,j,4], obstacle.xy_vector[1:i,2,j,4], color = "red", linestyle= ":")
                    obstacle_plot4 = ax10[:plot](obstacle.xy_vector[i,1,j,4], obstacle.xy_vector[i,2,j,4], color = "red", marker="o")  
                    y_obst_plot4 = ax10[:plot]([obstacle.axis_y_up[i,1,j,4],obstacle.axis_y_down[i,1,j,4]],[obstacle.axis_y_up[i,2,j,4],obstacle.axis_y_down[i,2,j,4]],color = "red")#plot the y semi axis
                    s_obst_plot4 = ax10[:plot]([obstacle.axis_s_up[i,1,j,4],obstacle.axis_s_down[i,1,j,4]],[obstacle.axis_s_up[i,2,j,4],obstacle.axis_s_down[i,2,j,4]],color = "red")# plot the s semi axis

                    plt_obst4[:remove]()
                    obst4 = patches.Ellipse([obstacle.xy_vector[i,1,j,4],obstacle.xy_vector[i,2,j,4]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[i,j,4]*180/pi)
                    plt_obst4 = ax10[:add_patch](obst4)
                end
                if obstacle.n_obstacle >=5
                    # obsttraj_plot5[1][:remove]()
                    obstacle_plot5[1][:remove]()
                    y_obst_plot5[1][:remove]()
                    s_obst_plot5[1][:remove]()
                    # obsttraj_plot5 = ax10[:plot](obstacle.xy_vector[1:i,1,j,5], obstacle.xy_vector[1:i,2,j,5], color = "red", linestyle= ":")
                    obstacle_plot5 = ax10[:plot](obstacle.xy_vector[i,1,j,5], obstacle.xy_vector[i,2,j,5], color = "red", marker="o")  
                    y_obst_plot5 = ax10[:plot]([obstacle.axis_y_up[i,1,j,5],obstacle.axis_y_down[i,1,j,5]],[obstacle.axis_y_up[i,2,j,5],obstacle.axis_y_down[i,2,j,5]],color = "red")#plot the y semi axis
                    s_obst_plot5 = ax10[:plot]([obstacle.axis_s_up[i,1,j,5],obstacle.axis_s_down[i,1,j,5]],[obstacle.axis_s_up[i,2,j,5],obstacle.axis_s_down[i,2,j,5]],color = "red")# plot the s semi axis

                    plt_obst5[:remove]()
                    obst5 = patches.Ellipse([obstacle.xy_vector[i,1,j,5],obstacle.xy_vector[i,2,j,5]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[i,j,5]*180/pi)
                    plt_obst5 = ax10[:add_patch](obst5)
                end
                if obstacle.n_obstacle >=6
                    # obsttraj_plot6[1][:remove]()
                    obstacle_plot6[1][:remove]()
                    y_obst_plot6[1][:remove]()
                    s_obst_plot6[1][:remove]()
                    # obsttraj_plot6 = ax10[:plot](obstacle.xy_vector[1:i,1,j,6], obstacle.xy_vector[1:i,2,j,6], color = "red", linestyle= ":")
                    obstacle_plot6 = ax10[:plot](obstacle.xy_vector[i,1,j,6], obstacle.xy_vector[i,2,j,6], color = "red", marker="o")  
                    y_obst_plot6 = ax10[:plot]([obstacle.axis_y_up[i,1,j,6],obstacle.axis_y_down[i,1,j,6]],[obstacle.axis_y_up[i,2,j,6],obstacle.axis_y_down[i,2,j,6]],color = "red")#plot the y semi axis
                    s_obst_plot6 = ax10[:plot]([obstacle.axis_s_up[i,1,j,6],obstacle.axis_s_down[i,1,j,6]],[obstacle.axis_s_up[i,2,j,6],obstacle.axis_s_down[i,2,j,6]],color = "red")# plot the s semi axis

                    plt_obst6[:remove]()
                    obst6 = patches.Ellipse([obstacle.xy_vector[i,1,j,6],obstacle.xy_vector[i,2,j,6]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[i,j,6]*180/pi)
                    plt_obst6 = ax10[:add_patch](obst6)
                end
                if obstacle.n_obstacle >=7
                    # obsttraj_plot7[1][:remove]()
                    obstacle_plot7[1][:remove]()
                    y_obst_plot7[1][:remove]()
                    s_obst_plot7[1][:remove]()
                    # obsttraj_plot7 = ax10[:plot](obstacle.xy_vector[1:i,1,j,7], obstacle.xy_vector[1:i,2,j,7], color = "red", linestyle= ":")
                    obstacle_plot7 = ax10[:plot](obstacle.xy_vector[i,1,j,7], obstacle.xy_vector[i,2,j,7], color = "red", marker="o")  
                    y_obst_plot7 = ax10[:plot]([obstacle.axis_y_up[i,1,j,7],obstacle.axis_y_down[i,1,j,7]],[obstacle.axis_y_up[i,2,j,7],obstacle.axis_y_down[i,2,j,7]],color = "red")#plot the y semi axis
                    s_obst_plot7 = ax10[:plot]([obstacle.axis_s_up[i,1,j,7],obstacle.axis_s_down[i,1,j,7]],[obstacle.axis_s_up[i,2,j,7],obstacle.axis_s_down[i,2,j,7]],color = "red")# plot the s semi axis

                    plt_obst7[:remove]()
                    obst7 = patches.Ellipse([obstacle.xy_vector[i,1,j,7],obstacle.xy_vector[i,2,j,7]],2*obstacle.rs,2*obstacle.ry,color="red",alpha=1.0,fill=false, angle = obstOrientation[i,j,7]*180/pi)
                    plt_obst7 = ax10[:add_patch](obst7)
                end
            end
            ##plot a line in the cost function that always show the current s
            if plot_costs ==1
                act_s_plot[1][:remove]()
                act_s_plot= ax9[:plot]([oldTraj.oldTraj[i,1,j],oldTraj.oldTraj[i,1,j]],[-3,oldTraj.oldNIter[j]], color ="black")  
            end
            if plot_lambda ==1
                act_lambda_plot[1][:remove]()
                act_lambda_plot= ax11[:plot]([oldTraj.oldTraj[i,1,j],oldTraj.oldTraj[i,1,j]],[0,1],label="_nolegend_", color ="black")  
            end

            # #plot inputs and their predictions
            if plot_inputs ==1 
                inta_plot_a[1][:remove]()
                inta_plot_df[1][:remove]()
                inta_plot_a = ax_Inp1[:plot](oldTraj.z_pred_sol[1:end-1,1,i,j], oldTraj.u_pred_sol[:,1,i,j] ,marker="o", color = "yellow")
                inta_plot_df = ax_Inp2[:plot](oldTraj.z_pred_sol[1:end-1,1,i,j], oldTraj.u_pred_sol[:,2,i,j] ,marker="o", color = "yellow")   
            end

            println("Press Enter for next plot step")
            println("Press c to cancel plotting of this round")
            a = ' '
            a = readline()
            if a == "c\r\n"
            stop_interactive = true 
                 break
            end
        end # end for loop interactive plots
        end #end if statement interactive plots
        if stop_interactive == true
            println("Do you want to cancel the whole plotting process? (c)")
            b = ' '
            b = readline()
            if b == "c\r\n"
                break
            end
        end
    end
    #return mpcParams, oldTraj
#end
