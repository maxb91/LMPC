using JLD
using PyPlot
function plots(j::Int64 = 1)
    
    include("helper/classes.jl")
    Data = load("2016-11-15-Data.jld")
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

    include("helper/calculateObstacleXY.jl")


xy_track  = [x_track; y_track]
t   = collect(0:dt:(buffersize-1)*dt)
xy_pred = zeros(mpcParams.N+1,2,length(t),n_rounds)

println("Number of simulated rounds in data: $n_rounds")

for k = 1:n_rounds
    for i = 1:i_final[j]
      
        xy_pred[:,:,i,j] = calculatePredictedXY(z_pred_log[:,:,:,k], mpcParams, trackCoeff, xy_track, convert(Int64,i))

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
figure(1)
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
# plot(t[1:i_final[j]-1],cost[1:i_final[j]-1,1],"r",t[1:i_final[j]-1],cost[1:i_final[j]-1,3],"b",t[1:i_final[j]-1],cost[1:i_final[j]-1,4],"y",t[1:i_final[j]-1],cost[1:i_final[j]-1,5],"m",t[1:i_final[j]-1],cost[1:i_final[j]-1,6],"c", t[1:i_final[j]-1], cost[1:i_final[j]-1,7])
# grid(1)
# title("Cost distribution")
# legend(["z","z_Term_const","deriv","control","lane", "Obstacle"], bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)


####plot of Cost in seperated axis
figure(2)
clf()
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



 # ####T
if j == 1
   
    # figure(4)
    # for i = 1:i_final[j]
    #     clf()
    #     ax1=subplot(311)
    #     plot(zCurr_s[1:i_final[j],1],zCurr_s[1:i_final[j],4])
    #     ax1=subplot(312)
    #     plot(zCurr_s[1:i_final[j],1],zCurr_s[1:i_final[j],2])
    #     ax1=subplot(313)
    #     plot(zCurr_s[1:i_final[j],1],zCurr_s[1:i_final[j],3])
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
car_plot = ax10[:plot](xStates_log[:,1,j], xStates_log[:,2,j], color = "blue")
pred_plot = ax10[:plot](xy_pred[:,1,1,j],xy_pred[:,2,1,j],color = "yellow", marker="o") 
obstacle_plot = ax10[:plot](obstacle.xy_vector[1,1,j], obstacle.xy_vector[1,2,j], color = "red",marker="o")
y_obst_plot = ax10[:plot]([obstacle.axis_y_up[1,1,j],obstacle.axis_y_down[1,1,j]],[obstacle.axis_y_up[1,2,j],obstacle.axis_y_down[1,2,j]],color = "red")#plot the y semi axis
s_obst_plot = ax10[:plot]([obstacle.axis_s_up[1,1,j],obstacle.axis_s_down[1,1,j]],[obstacle.axis_s_up[1,2,j],obstacle.axis_s_down[1,2,j]],color = "red")# plot the s semi axis
ax10[:plot](boundary_up[1,:], boundary_up[2,:],color="green",linestyle="--")
ax10[:plot](boundary_down[1,:], boundary_down[2,:],color="green",linestyle="--")
gca()[:set_aspect]("equal", adjustable="box")
grid()     
##########

#not everything necessary if we do the live plot
#plot of trajectory of car and obstacle
# ax10[:plot](zCurr_x[1:i_final[j],1], zCurr_x[1:i_final[j],2], color = "blue", linewidth = 1)# plot the trajectory of the car
# ax10[:plot](obstacle.xy_vector[1:i_final[j],1], obstacle.xy_vector[1:i_final[j],2], color = "red")
# car_plot[1][:remove]()
# obstacle_plot[1][:remove]()
# car_plot = ax10[:plot](zCurr_x[i_final[j],1], zCurr_x[i_final[j],2], color = "blue", marker="o")
# obstacle_plot = ax10[:plot](obstacle.xy_vector[i_final[j],1], obstacle.xy_vector[i_final[j],2], color = "red", marker="o")  

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

if j >= 2
    
    if j >= 3
        figure(7)
        clf()
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

    figure(4)   
    clf()
    ax4 = subplot(311)
    plot(sStates_log[1:i_final[j],1,1], sStates_log[1:i_final[j],4,1], color= "blue")
    plot(sStates_log[1:i_final[j],1,2], sStates_log[1:i_final[j],4,2], color= "yellow")
    plot(sStates_log[1:i_final[j],1,3], sStates_log[1:i_final[j],4,3], color= "green")
    grid()
    xlabel("s in [m]")
    ylabel("v in [m/s]")
    legend(["v current round","v old round", "3rd last round"])
    p1 = plot(1,1)

    ax5 = subplot(312, sharex=ax4)
    plot(sStates_log[1:i_final[j],1,1], sStates_log[1:i_final[j],2,1], color= "blue")
    plot(sStates_log[1:i_final[j],1,2], sStates_log[1:i_final[j],2,2], color= "yellow")
    plot(sStates_log[1:i_final[j],1,3], sStates_log[1:i_final[j],2,3], color= "green")
    grid()
    xlabel("s in [m]")
    ylabel("e_y in [m]")
    legend(["e_y current round","e_y old round", "3rd last round"])
    p2 = ax5[:plot](1,1)

    ax6 = subplot(313, sharex=ax4)
    hold(true)
    plot(sStates_log[1:i_final[j],1,1], sStates_log[1:i_final[j],3,1], color= "blue")
    plot(sStates_log[1:i_final[j],1,2], sStates_log[1:i_final[j],3,2], color= "yellow")
    plot(sStates_log[1:i_final[j],1,3], sStates_log[1:i_final[j],3,3], color= "green")
    grid()
    xlabel("s in [m]")
    ylabel("e_psi in [rad]")
    legend(["e_psi current round","e_psi old round", "3rd last round"])
    p3 = ax6[:plot](1,1)

    
    
    
    for i = 1:i_final[j]


        #plot predicted states over s
      
        p1[1][:remove]()
        p1 = ax4[:plot](z_pred_log[:,1,i,j], z_pred_log[:,4,i,j] ,marker="o") #plot predicted velocity

        p2[1][:remove]()
        p2 = ax5[:plot](z_pred_log[:,1,i,j], z_pred_log[:,2,i,j] ,marker="o") #plot predicted e_y
     
        p3[1][:remove]()
        p3 = ax6[:plot](z_pred_log[:,1,i,j], z_pred_log[:,3,i,j] ,marker="o")  #plot predicted e_psi
        
       
       #x-y plot
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


        # act_s_plot[1][:remove]()
        # act_s_plot= ax9[:plot]([oldTraj.oldTraj[i,1,1],oldTraj.oldTraj[i,1,1]],[0,30])  
        # #plot inputs
        # figure(5)
        # clf()
        # subplot(211) 
        # plot(oldTraj.oldTraj[1:i_final[j]-1,1,1], oldTraj.oldInput[1:i_final[j]-1,1,1], color= "green")
        # plot(z_log[1:end-1,1,i], u_log[:,1,i] ,marker="o")
        # xlabel("s in [m]")
        # ylabel("acceleration in [m/s^2]")
        # legend(["applied acceleration","predicted acceleration"])
        # grid()
        # subplot(212)
        #  plot(oldTraj.oldTraj[1:i_final[j]-1,1,1], oldTraj.oldInput[1:i_final[j]-1,2,1], color= "green")  
        # plot(z_log[1:end-1,1,i], u_log[:,2,i] ,marker="o")  
        #  xlabel("s in [m]")
        # ylabel("steering angle in [rad]")
        # legend(["applied steering angle","predicted steering angle"])
        # grid() 


    
        # figure(6)
        # plot(zCurr_s[1:i_final[j],1], zCurr_x[1:i_final[j],3])
        # xlabel("s in [m]")
        # ylabel("absolute angle psi in [rad]")
        # legend(["absolute angle psi plotted over s"]) 

        i = i+20



        println("Press Enter for next plot step")
        println("Press c to cancel plot")
        a = ' '
        a = readline()
        if a == "c\r\n" 
             break
        end
    end
end
end