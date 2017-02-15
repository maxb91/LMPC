# Minimal-Example on how to do Julia animations with subplots and some static content using matplotlib
using PyCall
using PyPlot
using JLD
close("all")
@pyimport matplotlib as mpl
@pyimport matplotlib.patches as patches
@pyimport matplotlib.animation as animation# First set up the figure, the axis, and the plot element we want to animate
matplotlib[:style][:use]("classic") # somehow my julia version changed plotting style 

include("helper/classes.jl")
include("helper/plot_functions.jl")
include("helper/calculateObstacleXY.jl")
include("helper/colorModule.jl")


file = "data/2017-02-14-16-34-Data.jld"
copied_plot = 1
j=1

obstacle_color = "red"
boundary_color = "black"
ego_color = "green"



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

println("Number of simulated rounds in data: $(oldTraj.n_oldTraj)")
println("Load File located in: $file")

obstacle.rs = 0.2
obstacle.ry = 0.1

ds = trackCoeff.ds
dt = modelParams.dt
xy_track  = [x_track; y_track]
t   = collect(0:dt:(buffersize-1)*dt)
trackLength = (size(x_track)[2]-1)*trackCoeff.ds



boundary_up, boundary_down, trackL = calc_track_boundaryXY(xy_track, trackCoeff)

################################
##calculate predcted position of the car in XY plane
# ################################
xy_pred = zeros(mpcParams.N+1,2,size(t)[1],oldTraj.n_oldTraj)
obstOrientation = zeros(size(obstacle.xy_vector)[1],size(obstacle.xy_vector)[3],size(obstacle.xy_vector)[4])
for k = 1:oldTraj.n_oldTraj
    for i = 1:oldTraj.oldNIter[k]
        # caluclate the predicted XY postion of the car from the s-ey values
        xy_pred[:,:,i,k] = calculatePredictedXY(oldTraj.z_pred_sol[:,:,:,k], mpcParams, trackCoeff, xy_track, convert(Int64,i),k)
        #calculate the obstacle postion from the s-ey values
        obstOrientation = calculateObstacleXY!(obstacle, trackCoeff, xy_track,i,k,obstOrientation ) #this funciton computes values for row i
    end
end

#oval limits
# xmin = -15
# xmax = 15
# ymin = -15
# ymax = 15

xmin = -10
xmax = 13
ymin = -13
ymax = 10
# create figure to plot
fig = figure(figsize=(15,8))
ax1= subplot2grid([2,3],[0,0], colspan =2, rowspan=2)
ax1[:set_xlim](xmin, xmax)
ax1[:set_ylim](ymin, ymax)
# set labels
ax1[:set_xlabel](L"x",fontsize=20)
ax1[:set_ylabel](L"y",fontsize=20)
ax1[:tick_params]( axis="x", which="both", bottom="off", top="off", labelbottom="off")
ax1[:tick_params]( axis="y", which="both", left="off", right="off", labelleft="off")
ax1[:set_aspect]("equal", adjustable="box")
#ax1[:canvas][:set_window_title]("Track and cars in XY plane")
# grid()
ax2= subplot2grid([2,3],[0,2])
ax2[:set_ylabel](L"v \> [m/s]",fontsize=20)
# ax2[:set_xlim](findmin(oldTraj.copyInfo[:,3,j])[1], findmax(oldTraj.copyInfo[:,3,j])[1])
# ax2[:set_ylim](findmin(oldTraj.copyInfo[:,2,j])[1], findmax(oldTraj.copyInfo[:,2,j])[1])
ax2[:set_xlim](0, trackLength)
ax2[:set_ylim](0, modelParams.v_max+0.1)
grid()#


ax3= subplot2grid([2,3],[1,2])
ax3[:set_xlabel](L"s \> [m]",fontsize=20)
ax3[:set_xlim](0, trackLength)
if copied_plot ==1
ax3[:set_ylabel](L"s \> [m]",fontsize=20)
ax3[:set_ylim](0, trackLength)
else 
    ax3[:set_ylabel](L"obstacle \> cost",fontsize=20)
end
# ax3[:set_ylim](0, 35)
grid()
# fig[:subplots_adjust](hspace=.5)

###########################
#x-y plot
###########################

# plot the boundary lines
ax1[:plot](x_track',y_track', linestyle = (0, (4.0, 8.0)), color = color=boundary_color, linewidth = 0.4, label="_nolegend_")#plot the racetrack
ax1[:plot](boundary_up[1,:], boundary_up[2,:],color=boundary_color, linewidth = 0.7)#,linestyle=":")
ax1[:plot](boundary_down[1,:], boundary_down[2,:],color=boundary_color, linewidth = 0.7)#,linestyle=":")
for l=1:convert(Int64,trunc(trackL/51))
    ax1[:plot]([boundary_down[1,l*50+1],boundary_up[1,l*50+1]],[boundary_down[2,l*50+1],boundary_up[2,l*50+1]], color = "black", linestyle = ":", linewidth = 0.5)
    boundvec = [boundary_up[1,l*50+1]-boundary_down[1,l*50+1];boundary_up[2,l*50+1]-boundary_down[2,l*50+1]]
    ax1[:text](boundary_down[1,l*50+1]+1.85*boundvec[1],boundary_down[2,l*50+1]+1.85*boundvec[2],"$(convert(Int64,l*50*ds))",fontsize=20,clip_on = true)
end


#obstacle and car
obst_patch =Array{PyCall.PyObject}(obstacle.n_obstacle)
plt_obst =Array{PyCall.PyObject}(obstacle.n_obstacle)
for ii = 1:obstacle.n_obstacle
    obst_patch[ii] = patches.Ellipse([obstacle.xy_vector[1,1,j,ii],obstacle.xy_vector[1,2,j,ii]],2*obstacle.rs,2*obstacle.ry,color=obstacle_color,alpha=1.0,fill=false, angle = obstOrientation[1,j,ii]*180/pi)
    plt_obst[ii] = ax1[:add_patch](obst_patch[ii])
end

carParts=drawCar(ax1,[oldTraj.oldTrajXY[2,1,j],oldTraj.oldTrajXY[2,2,j],oldTraj.oldTrajXY[2,5,j]*180/pi], ego_color)
for p in carParts
                ax1[:add_patch](p)
end
car_plot = ax1[:plot]([],[], color = ego_color)[1]
pred_plot = ax1[:plot](xy_pred[:,1,2,j],xy_pred[:,2,2,j],color = "yellow", marker="o")[1]#(oldTraj.oldTrajXY[1,1,j],oldTraj.oldTrajXY[1,2,j],color = "yellow", marker="o")


# last trajectory in xy plot
# lastj_plot = ax1[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+1],1,j+1], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+1],2,j+1], color = ego_color, linewidth = 0.7, linestyle = ":")
####################
#copied-plot
####################
if copied_plot == 1
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
    for i=1:oldTraj.oldNIter[j]-1
        if oldTraj.copyInfo[i,1,j]>0.0 # if lambda of copied traj is greater 0.1 -> if traj s used for solving.
            ax3[:scatter](oldTraj.copyInfo[i,3,j],oldTraj.copyInfo[i,2,j], color = "#DCDCDC",edgecolor="black", lw =0.5)
        end
    end
    for i=1:oldTraj.oldNIter[j]-1
        if oldTraj.copyInfo[i,4,j]>0.6 # if lambda of copied traj is greater 0.1 -> if traj s used for solving.
            ax3[:scatter](oldTraj.copyInfo[i,3,j],oldTraj.copyInfo[i,2,j], color = "blue",edgecolor="black",lw =0.5) #color = colordefs[convert(Int64,oldTraj.copyInfo[i,1,j])]
        end
    end
    s_curr_copied = ax3[:plot]([],[], color ="black")[1]
else
####################
#cost plot
####################
    ax3[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[8,1:oldTraj.oldNIter[j]-1,j])
    s_curr_copied = ax3[:plot]([],[], color ="black")[1]
    max_cost = findmax(oldTraj.costs[8,1:oldTraj.oldNIter[j]-1,j])[1]
end


####################
#v-plot
####################
ax2[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTraj[1:oldTraj.oldNIter[j],4,j], color = ego_color)
v_plot  = ax2[:plot](oldTraj.z_pred_sol[:,1,1,j], oldTraj.z_pred_sol[:,4,1,j],marker="o", color = "yellow")[1]
# s_curr_v = ax2[:plot]([],[], color ="black")[1]


####################
#Init function
####################
# Function is called by animator to create the non-animated "background" of the animation
# important: return the line and a nothing at the en

# the function is called by the animator to change the plot elements that are returned. I basically change the
# x and y data of lines 1 to 3.
# important: return the lines + a nothing
# k is automatically incremented everytime the animate function fires. k=k+1 accounts for the fact that python starts with 0 where julia starts with 1
        # tight_layout()
start_offset = 50
function animate(frame)
    frame+=1
    
    k = frame-start_offset
    # if  i >= oldTraj.oldNIter[j]
    #     global j =6
    #    k  =l 
    #    global l = l+1
    # end
    l = frame-10
    if frame <= 10
        ax1[:set_xlim](xmin, xmax)
        ax1[:set_ylim](ymin, ymax)
    elseif frame>=11 && frame <= 51
        #oval
        # ax1[:set_xlim](xmin+13/40*l, xmax-11/40*l)
        # ax1[:set_ylim](ymin+12/40*l, ymax-12/40*l)

        ax1[:set_xlim](xmin+8/40*l, xmax-9/40*l)
        ax1[:set_ylim](ymin+10/40*l, ymax-7/40*l)
    else
        ##################
        #x-y plot
        ##################
        #update the position of the car and the prediction
        for ii=1:obstacle.n_obstacle
            plt_obst[ii][:remove]()
            obst_patch[ii] = patches.Ellipse([obstacle.xy_vector[k,1,j,ii],obstacle.xy_vector[k,2,j,ii]],2*obstacle.rs,2*obstacle.ry,color=obstacle_color,alpha=1.0,fill=false, angle = obstOrientation[k,j,ii]*180/pi)
            plt_obst[ii] = ax1[:add_patch](obst_patch[ii])
        end

    #     #update axis of x-y plot every time step
        ax1[:set_xlim](xy_pred[6,1,k,j]-3, xy_pred[6,1,k,j]+3)
        ax1[:set_ylim](xy_pred[6,2,k,j]-3, xy_pred[6,2,k,j]+3)


        car_plot[:set_data](oldTraj.oldTrajXY[1:k,1,j], oldTraj.oldTrajXY[1:k,2,j])

        #plot patches for current car position
        pos =[oldTraj.oldTrajXY[k,1,j],oldTraj.oldTrajXY[k,2,j], oldTraj.oldTrajXY[k,5,j]*180/pi]
        updateCarParts(ax1,carParts,pos)

        pred_plot[:set_data](xy_pred[:,1,k,j],xy_pred[:,2,k,j])
        # plot(oldTraj.oldTrajXY[k,1,j],oldTraj.oldTrajXY[k,2,j])
        if copied_plot ==1
    #     #######################
    #     #copied plot
    #     #######################
            s_curr_copied[:set_data]([oldTraj.oldTraj[k,1,j],oldTraj.oldTraj[k,1,j]],[0,trackLength])
        else
            s_curr_copied[:set_data]([oldTraj.oldTraj[k,1,j],oldTraj.oldTraj[k,1,j]],[0,max_cost])
        end
    #     #######################
    #     #v plot
    #     #######################
        v_plot[:set_data](oldTraj.z_pred_sol[:,1,k,j], oldTraj.z_pred_sol[:,4,k,j])
    end
    return carParts, nothing
end
# animator object runs the animation
# frames: defines how often animate function is called and interval how often in milliseconds

# anim = animation.FuncAnimation(fig, animate,frames=oldTraj.oldNIter[j]+oldTraj.oldNIter[j-1]+start_offset-1, interval=1,repeat=false )
# oldTraj.oldNIter[j]
# anim = animation.FuncAnimation(fig, animate,frames=oldTraj.oldNIter[j], interval=10,repeat=false)
anim = animation.FuncAnimation(fig, animate,frames=oldTraj.oldNIter[j]+start_offset-1, interval=10,repeat=false)#, blit=true)
# to save video: you need encoder and writer, for ubuntu it should work if you do: sudo apt-get install ffmpeg
FFwriter = animation.FFMpegWriter( fps=10, bitrate=-1, extra_args=["-vcodec", "libx264","-pix_fmt","yuv420p"])   #(fps=10,bitrate=3000,extra_args=["-vcodec", "libx264","-pix_fmt","yuv420p"])
anim[:save]("basic_animation.mp4", writer = FFwriter)
