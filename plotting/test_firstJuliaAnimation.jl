# Minimal-Example on how to do Julia animations with subplots and some static content using matplotlib
using PyCall
using PyPlot
using JLD
n_oldTrajPlots = 1
global j  =7
include("../helper/classes.jl")
file = "data/2017-02-16-02-21-Data.jld"
close("all")

    ####load data from file
    Data = load(file)
global l  =1
   
    x_track = Data["x_track"]
    y_track = Data["y_track"]
    trackCoeff = Data["trackCoeff"]
    obstacle = Data["obstacle"]
    modelParams = Data["modelParams"]
    mpcParams = Data["mpcParams"]
    buffersize = Data["buffersize"]
    # curv_approx = Data["curv_approx"]
    oldTraj     = Data["oldTraj"]
    #include("calculateObstacleXY.jl")
    #include("colorModule.jl")
    include("calculateObstacleXY.jl")
    include("colorModule.jl")
    #end of data loading

    #####create additional data for plotting
    dt = modelParams.dt
    xy_track  = [x_track; y_track]
    t   = collect(0:dt:(buffersize-1)*dt)
    xy_pred = zeros(mpcParams.N+1,2,length(t),oldTraj.n_oldTraj)

    println("Number of simulated rounds in data: $(oldTraj.n_oldTraj)")
    println("Load File located in: $file")
    obstOrientation = zeros(size(obstacle.xy_vector)[1],size(obstacle.xy_vector)[3],size(obstacle.xy_vector)[4])
    for k = 1:oldTraj.n_oldTraj
        for i = 1:oldTraj.oldNIter[j]
            # caluclate the predicted XY postion of the car from the s-ey values
            xy_pred[:,:,i,k] = calculatePredictedXY(oldTraj.z_pred_sol[:,:,:,k], mpcParams, trackCoeff, xy_track, convert(Int64,i),k)
            #calculate the obstacle postion from the s-ey values
            calculateObstacleXY!(obstacle, trackCoeff, xy_track,i,k,obstOrientation) #this funciton computes values for row i
        end
    end

    ##this part is to calculate the tracks boundaries and plot them later
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







@pyimport matplotlib.animation as animation# First set up the figure, the axis, and the plot element we want to animate
#@pyimport matplotlib.pyplot as plt
#plt[:rcParams]["animation.ffmpeg_path"] ="C:/ffmpeg-3.2-win64-static/bin/"
# create dummy data


# define axis limits
xmin = 0
xmax = 20
ymin = 0
ymax = 15

# create figure to plot
fig = figure(figsize=(15,8))
# gs = matplotlib[:gridspec][:GridSpec](3, 2) #width_ratios=[3, 1])  
# ax1 = fig[:add_subplot](gs[1,1:2])
# add two subplots to the figure and define axis limits
#ax1 = fig[:add_subplot](1, 2, 1)
 ax1= subplot2grid([3,2],[0,0], rowspan=3)
ax1[:set_xlim](xmin, xmax)
ax1[:set_ylim](ymin, ymax)
# set labels
ax1[:set_xlabel](L"x",fontsize=20)
ax1[:set_ylabel](L"y",fontsize=20)
ax1[:set_aspect]("equal", adjustable="box")
#ax1[:canvas][:set_window_title]("Track and cars in XY plane")
grid()
ax2= subplot2grid([3,2],[0,1])
ax2[:set_ylabel](L"v",fontsize=20)
grid()#
ax3= subplot2grid([3,2],[1,1], sharex = ax2)
ax3[:set_ylabel](L"e_y",fontsize=20)
grid()
ax4= subplot2grid([3,2],[2,1], sharex = ax2)
ax4[:set_xlabel](L"s",fontsize=20)

grid()

###########################
#x-y plot
###########################
colorObjectXY= colorModule.ColorManager()
for k= 1:n_oldTrajPlots
    colorXYold = colorModule.getColor(colorObjectXY)
    ax1[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+k],1,j+k], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+k],2,j+k],linestyle="--", color = colorXYold, label= "$k old Traj")[1]
end
 #plot s length in 5m steps over track
for l=1:8
    ax1[:plot]([boundary_down[1,l*50+1],boundary_up[1,l*50+1]],[boundary_down[2,l*50+1],boundary_up[2,l*50+1]], color = "black", linestyle = ":", linewidth = 0.5)
    ax1[:text](boundary_down[1,l*50+1]+1,boundary_down[2,l*50+1],"s = $(l*5)",fontsize=8,clip_on = true)
end
# plot trajectory less visible to overwrite it with interactiv simulation
ax1[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],2,j], color = "black" ,linestyle=":", label="current Traj")

# define the lines for the subplots. Their data will later be changed creating the animation
 #plot the boundary lines
ax1[:plot](boundary_up[1,:], boundary_up[2,:],color="green",linestyle=":")
ax1[:plot](boundary_down[1,:], boundary_down[2,:],color="green",linestyle=":")
middle = ax1[:plot]([],[],linestyle="--", color = "#FF9900", linewidth = 0.5, label="_nolegend_")[1]
pred_plot = ax1[:plot]([],[],color = "yellow", marker="o")[1]
car_plot = ax1[:plot]([],[], color = "black")[1]
car_plot_a = ax1[:plot]([],[], color = "black", marker = "o")[1]

obstacle_plot = ax1[:plot]([],[], color = "red",linestyle=":", label ="_nolegend_")[1]
obstacle_plot_a = ax1[:plot]([],[], color = "red",marker="o", label = "obstacle Traj")[1]
y_obst_plot = ax1[:plot]([],[],color = "red")[1]#plot the y semi axis
s_obst_plot = ax1[:plot]([],[],color = "red")[1]# plot the s semi axis

####################
#v-plot
####################


k = j+1
colorObjectV= colorModule.ColorManager()
while k<=j+n_oldTrajPlots
    colorV = colorModule.getColor(colorObjectV)
    ax2[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[k],1,k], oldTraj.oldTraj[1:oldTraj.oldNIter[k],4,k], color = colorV)
    k  = k+1
end
ax2[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTraj[1:oldTraj.oldNIter[j],4,j], color = "black")
v_plot  = ax2[:plot]([],[],marker="o", color = "yellow")[1]
ax2[:set_ylim](0, 2.1)
ax2[:set_ylabel](L"v",fontsize=20)

####################
#lambda-plot
####################
# colorObjectLambda = colorModule.ColorManager()
# while k<= j+n_oldTrajPlots
#     colorLambda = colorModule.getColor(colorObjectLambda)
#     ax2[:scatter](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.lambda_sol[k,1:oldTraj.oldNIter[j]-1,j].*oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorLambda, marker = "x")
#     k = k+1
# end
# act_lambda_plot = ax2[:plot]([],[],color = "black", label="_nolegend_")[1]
# ax2[:set_ylim](-0.01, 1.01)
# ax2[:set_ylabel](L"\lambda",fontsize=20)
####################
#e_y-plot
####################

k = j+1
colorObject_eY= colorModule.ColorManager()
while k<=j+n_oldTrajPlots
    color_eY = colorModule.getColor(colorObject_eY)
    ax3[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[k],1,k], oldTraj.oldTraj[1:oldTraj.oldNIter[k],2,k], color = color_eY)
    k  = k+1
end
ax3[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTraj[1:oldTraj.oldNIter[j],2,j], color = "black")
e_y_plot  = ax3[:plot]([],[],marker="o", color = "yellow")[1]

####################
#obstacle cost plot
####################
ax4[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[2,1:oldTraj.oldNIter[j]-1,j])  
s_cost_plot= ax4[:plot]([],[], color ="black")[1]
ax4[:set_ylabel](L"cost_{term}",fontsize=20)

####################
#Init function
####################
# Function is called by animator to create the non-animated "background" of the animation
# important: return the line and a nothing at the end
function initPlot()
    ###########################
    #x-y plot
    ###########################
    #x-y Plot of the racetrack and obstacle and car
    middle[:set_data](x_track',y_track')
    pred_plot[:set_data](xy_pred[:,1,1,j],xy_pred[:,2,1,j])
    car_plot[:set_data](oldTraj.oldTrajXY[1,1,j], oldTraj.oldTrajXY[1,2,j])
    car_plot_a[:set_data](oldTraj.oldTrajXY[1,1,j], oldTraj.oldTrajXY[1,2,j])
    obstacle_plot[:set_data](obstacle.xy_vector[1,1,j], obstacle.xy_vector[1,2,j])
    obstacle_plot_a[:set_data](obstacle.xy_vector[1,1,j], obstacle.xy_vector[1,2,j])  
    y_obst_plot[:set_data]([obstacle.axis_y_up[1,1,j],obstacle.axis_y_down[1,1,j]],[obstacle.axis_y_up[1,2,j],obstacle.axis_y_down[1,2,j]])#plot the y semi axis
    s_obst_plot[:set_data]([obstacle.axis_s_up[1,1,j],obstacle.axis_s_down[1,1,j]],[obstacle.axis_s_up[1,2,j],obstacle.axis_s_down[1,2,j]])
    #######################
    #v plot
    #######################
    v_plot[:set_data](oldTraj.z_pred_sol[:,1,1,j], oldTraj.z_pred_sol[:,4,1,j]) #plot predicted velocity

    #######################
    #lambda plot
    #######################
    # act_lambda_plot[:set_data]([oldTraj.oldTraj[1,1,j],oldTraj.oldTraj[1,1,j]],[0,1])
    #######################
    #e_y plot
    #######################
    e_y_plot[:set_data](oldTraj.z_pred_sol[:,1,1,j], oldTraj.z_pred_sol[:,2,1,j]) #plot predicted e_y    
    #######################
    #obstacle cost plot
    #######################
    s_cost_plot[:set_data]([oldTraj.oldTraj[1,1,j],oldTraj.oldTraj[1,1,j]],[-3,oldTraj.oldNIter[j]])
       
        # gca()[:set_aspect]("equal", adjustable="box")
        # grid() 
        #legend(bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)    
    return(middle,pred_plot,car_plot,car_plot_a,obstacle_plot,obstacle_plot_a,y_obst_plot,s_obst_plot,e_y_plot,s_cost_plot, nothing)
end

# the function is called by the animator to change the plot elements that are returned. I basically change the
# x and y data of lines 1 to 3.
# important: return the lines + a nothing
# k is automatically incremented everytime the animate function fires. k=k+1 accounts for the fact that python starts with 0 where julia starts with 1
start_offset = -10
function animate(i)
    i = i+5
    
    # k = i-start_offset
    k = i
    if  i >= oldTraj.oldNIter[j]
        global j =6
       k  =l 
       global l = l+1
    end
#     l = i-20
# if i <= 20
#     ax1[:set_xlim](xmin, xmax)
#     ax1[:set_ylim](ymin, ymax)
# elseif i>=20 && i <= 40
#     ax1[:set_xlim](xmin, xmax-0.75*l)
#     ax1[:set_ylim](ymin, ymax-0.51*l)
# else
    ##################
    #x-y plot
    ##################
    #update the position of the car and the prediction
    pred_plot[:set_data](xy_pred[:,1,k,j],xy_pred[:,2,k,j])
    car_plot[:set_data](oldTraj.oldTrajXY[1:k,1,j], oldTraj.oldTrajXY[1:k,2,j])
    car_plot_a[:set_data](oldTraj.oldTrajXY[k,1,j], oldTraj.oldTrajXY[k,2,j]) #plot current position of car with a marker

    #update the obstacle position
    obstacle_plot[:set_data](obstacle.xy_vector[1:k,1,j], obstacle.xy_vector[1:k,2,j])
    obstacle_plot_a[:set_data](obstacle.xy_vector[k,1,j], obstacle.xy_vector[k,2,j])  
    y_obst_plot[:set_data]([obstacle.axis_y_up[k,1,j],obstacle.axis_y_down[k,1,j]],[obstacle.axis_y_up[k,2,j],obstacle.axis_y_down[k,2,j]])#plot the y semi axis
    s_obst_plot[:set_data]([obstacle.axis_s_up[k,1,j],obstacle.axis_s_down[k,1,j]],[obstacle.axis_s_up[k,2,j],obstacle.axis_s_down[k,2,j]])

    #update axis of x-y plot every time step
    ax1[:set_xlim](xy_pred[6,1,k,j]-2, xy_pred[6,1,k,j]+2)
    ax1[:set_ylim](xy_pred[6,2,k,j]-2, xy_pred[6,2,k,j]+2)
    #######################
    #v plot
    #######################
    v_plot[:set_data](oldTraj.z_pred_sol[:,1,k,j], oldTraj.z_pred_sol[:,4,k,j]) #plot predicted velocity
    #######################
    #lambda plot
    #######################
    # act_lambda_plot[:set_data]([oldTraj.oldTraj[k,1,j],oldTraj.oldTraj[k,1,j]],[0,1])
    #######################
    #e_y plot
    #######################
    e_y_plot[:set_data](oldTraj.z_pred_sol[:,1,k,j], oldTraj.z_pred_sol[:,2,k,j]) #plot predicted e_y   
    #######################
    #obstacle cost plot
    #######################
    s_cost_plot[:set_data]([oldTraj.oldTraj[k,1,j],oldTraj.oldTraj[k,1,j]],[-3,oldTraj.oldNIter[j]])
# end
    return (pred_plot,car_plot,car_plot_a,obstacle_plot,obstacle_plot_a,y_obst_plot,s_obst_plot,e_y_plot,s_cost_plot,nothing)
end
# animator object runs the animation
# frames: defines how often animate function is called and interval how often in milliseconds

anim = animation.FuncAnimation(fig, animate, init_func=initPlot,frames=oldTraj.oldNIter[j]+oldTraj.oldNIter[j-1]+start_offset-1, interval=1,repeat=false )

# to save video: you need encoder and writer, for ubuntu it should work if you do: sudo apt-get install ffmpeg
FFwriter = animation.FFMpegWriter( fps=10, bitrate=1500, extra_args=["-vcodec", "libx264"])   #(fps=10,bitrate=3000,extra_args=["-vcodec", "libx264","-pix_fmt","yuv420p"])
anim[:save]("basic_animation.mp4", writer = FFwriter)
