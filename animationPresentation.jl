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


file = "data/2017-02-08-16-24-Data.jld"
global j  = 1
global l  = 1
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

ds = trackCoeff.ds
dt = modelParams.dt
xy_track  = [x_track; y_track]
t   = collect(0:dt:(buffersize-1)*dt)





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


xmin = -7
xmax = 12
ymin = -10
ymax = 2
# create figure to plot
fig = figure(figsize=(15,8))
ax1= subplot2grid([4,1],[0,0], rowspan=3)
ax1[:set_xlim](xmin, xmax)
ax1[:set_ylim](ymin, ymax)
# set labels
ax1[:set_xlabel](L"x",fontsize=20)
ax1[:set_ylabel](L"y",fontsize=20)
ax1[:set_aspect]("equal", adjustable="box")
#ax1[:canvas][:set_window_title]("Track and cars in XY plane")
grid()
ax2= subplot2grid([4,1],[3,0])
ax2[:set_ylabel](L"v",fontsize=20)
grid()#


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

# obst_patch =Array{PyCall.PyObject}(obstacle.n_obstacle)
# plt_obst =Array{PyCall.PyObject}(obstacle.n_obstacle)

carParts=drawCar(ax1,[0.0,0.0,0.0], "green")
for p in carParts
                ax1[:add_patch](p)
end


# cartraj_plot = ax1[:plot](1,1)  #just for initialization
# car_plot = ax1[:plot](1,1)#(oldTraj.oldTrajXY[1,1,j], oldTraj.oldTrajXY[1,2,j], color = ego_color) # just dummy to use remove func later
# pred_plot = ax1[:plot](1,1)#(oldTraj.oldTrajXY[1,1,j],oldTraj.oldTrajXY[1,2,j],color = "yellow", marker="o")
####################
#copied-plot
####################


####################
#e_y-plot
####################

# k
# e_y_plot  = ax3[:plot]([],[],marker="o", color = "yellow")[1]

####################
#obstacle cost plot
####################
# ax4[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[2,1:oldTraj.oldNIter[j]-1,j])  
# s_cost_plot= ax4[:plot]([],[], color ="black")[1]
# ax4[:set_ylabel](L"cost_{term}",fontsize=20)

####################
#Init function
####################
# Function is called by animator to create the non-animated "background" of the animation
# important: return the line and a nothing at the en

# the function is called by the animator to change the plot elements that are returned. I basically change the
# x and y data of lines 1 to 3.
# important: return the lines + a nothing
# k is automatically incremented everytime the animate function fires. k=k+1 accounts for the fact that python starts with 0 where julia starts with 1

function animate(frame)
    frame+=1
    
    # k = frame-start_offset
    k = frame
    # if  i >= oldTraj.oldNIter[j]
    #     global j =6
    #    k  =l 
    #    global l = l+1
    # end
#     l = frame-20
# if frame <= 20
#     ax1[:set_xlim](xmin, xmax)
#     ax1[:set_ylim](ymin, ymax)
# elseif frame>=20 && i <= 40
#     ax1[:set_xlim](xmin, xmax-0.75*l)
#     ax1[:set_ylim](ymin, ymax-0.51*l)
# else
    ##################
    #x-y plot
    ##################
    #update the position of the car and the prediction

    #update the obstacle position
#     obstacle_plot[:set_data](obstacle.xy_vector[1:k,1,j], obstacle.xy_vector[1:k,2,j])
# cle.axis_s_down[k,1,j]],[obstacle.axis_s_up[k,2,j],obstacle.axis_s_down[k,2,j]])
#     #update axis of x-y plot every time step
#     ax1[:set_xlim](xy_pred[6,1,k,j]-2, xy_pred[6,1,k,j]+2)
#     ax1[:set_ylim](xy_pred[6,2,k,j]-2, xy_pred[6,2,k,j]+2)
    pos =[oldTraj.oldTrajXY[k,1,j],oldTraj.oldTrajXY[k,2,j], oldTraj.oldTrajXY[k,5,j]*180/pi]
    updateCarParts(ax1,carParts,pos)
    # plot(oldTraj.oldTrajXY[k,1,j],oldTraj.oldTrajXY[k,2,j])
#     #######################
#     #copied plot
#     #######################
#     v_plot[:set_data](oldTraj.z_pred_sol[:,1,k,j], oldTraj.z_pred_sol[:,4,k,j]) #plot predicted velocity

#     #######################
#     #obstacle cost plot
#     #######################
#     s_cost_plot[:set_data]([oldTraj.oldTraj[k,1,j],oldTraj.oldTraj[k,1,j]],[-3,oldTraj.oldNIter[j]])
# end
    return carParts, nothing
end
# animator object runs the animation
# frames: defines how often animate function is called and interval how often in milliseconds

# anim = animation.FuncAnimation(fig, animate,frames=oldTraj.oldNIter[j]+oldTraj.oldNIter[j-1]+start_offset-1, interval=1,repeat=false )
anim = animation.FuncAnimation(fig, animate,frames=oldTraj.oldNIter[j], interval=1,repeat=false)#, blit=true)
# to save video: you need encoder and writer, for ubuntu it should work if you do: sudo apt-get install ffmpeg
FFwriter = animation.FFMpegWriter( fps=8, bitrate=-1, extra_args=["-vcodec", "libx264","-pix_fmt","yuv420p"])   #(fps=10,bitrate=3000,extra_args=["-vcodec", "libx264","-pix_fmt","yuv420p"])
anim[:save]("basic_animation.mp4", writer = FFwriter)
