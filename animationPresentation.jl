# Minimal-Example on how to do Julia animations with subplots and some static content using matplotlib
using PyCall
using PyPlot
using JLD
n_oldTrajPlots = 1
global j  =7
global l  =1

@pyimport matplotlib.animation as animation
@pyimport matplotlib.patches as patches
@pyimport matplotlib.animation as animation# First set up the figure, the axis, and the plot element we want to animate
xmin = 0
xmax = 20
ymin = 0
ymax = 15

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


####################
#copied-plot
####################




####################
#e_y-plot
####################

k
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
    pred_plot[:set_data](xy_pred[:,1,1,j],xy_pred[:,2,1,j])

    car_plot_a[:set_data](oldTraj.oldTrajXY[1,1,j], oldTraj.oldTrajXY[1,2,j])

    #######################
    #copied plot
    #######################
    v_plot[:set_data](oldTraj.z_pred_sol[:,1,1,j], oldTraj.z_pred_sol[:,4,1,j]) 

  
    return(, nothing)
end

# the function is called by the animator to change the plot elements that are returned. I basically change the
# x and y data of lines 1 to 3.
# important: return the lines + a nothing
# k is automatically incremented everytime the animate function fires. k=k+1 accounts for the fact that python starts with 0 where julia starts with 1
start_offset = -10
function animate(i)
    i = i+1
    
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


    #update the obstacle position
    obstacle_plot[:set_data](obstacle.xy_vector[1:k,1,j], obstacle.xy_vector[1:k,2,j])
cle.axis_s_down[k,1,j]],[obstacle.axis_s_up[k,2,j],obstacle.axis_s_down[k,2,j]])

    #update axis of x-y plot every time step
    ax1[:set_xlim](xy_pred[6,1,k,j]-2, xy_pred[6,1,k,j]+2)
    ax1[:set_ylim](xy_pred[6,2,k,j]-2, xy_pred[6,2,k,j]+2)




    #######################
    #copied plot
    #######################
    v_plot[:set_data](oldTraj.z_pred_sol[:,1,k,j], oldTraj.z_pred_sol[:,4,k,j]) #plot predicted velocity

    #######################
    #obstacle cost plot
    #######################
    s_cost_plot[:set_data]([oldTraj.oldTraj[k,1,j],oldTraj.oldTraj[k,1,j]],[-3,oldTraj.oldNIter[j]])
# end
    return (,nothing)
end
# animator object runs the animation
# frames: defines how often animate function is called and interval how often in milliseconds

anim = animation.FuncAnimation(fig, animate, init_func=initPlot,frames=oldTraj.oldNIter[j]+oldTraj.oldNIter[j-1]+start_offset-1, interval=1,repeat=false )

# to save video: you need encoder and writer, for ubuntu it should work if you do: sudo apt-get install ffmpeg
FFwriter = animation.FFMpegWriter( fps=8, bitrate=1500, extra_args=["-vcodec", "libx264"])   #(fps=10,bitrate=3000,extra_args=["-vcodec", "libx264","-pix_fmt","yuv420p"])
anim[:save]("basic_animation.mp4", writer = FFwriter)
