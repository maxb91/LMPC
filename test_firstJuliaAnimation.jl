# Minimal-Example on how to do Julia animations with subplots and some static content using matplotlib
n_oldTrajPlots = 6
j =1
file = "data/2016-12-05-23-40-Data.jld"
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
    #include("calculateObstacleXY.jl")
    #include("colorModule.jl")
    include("helper/calculateObstacleXY.jl")
    include("helper/colorModule.jl")
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






using PyCall
using PyPlot
@pyimport matplotlib.animation as animation# First set up the figure, the axis, and the plot element we want to animate
#@pyimport matplotlib.pyplot as plt
#plt[:rcParams]["animation.ffmpeg_path"] ="C:/ffmpeg-3.2-win64-static/bin/"
# create dummy data
xdata = 0:1:300
ydata = 2*xdata
zdata = 3*xdata

# define axis limits
xmin = 0
xmax = 20
ymin = 0
ymax = 15

# create figure to plot
fig = figure(figsize=(15,8))

# add two subplots to the figure and define axis limits
ax10 = fig[:add_subplot](1, 2, 1)
 ax10[:set_xlim](xmin, xmax)
 ax10[:set_ylim](ymin, ymax)
ax2 = fig[:add_subplot](1,2,2)
ax2[:set_xlim](xmin, xmax)
ax2[:set_ylim](ymin, ymax)

# set labels
# ax1[:set_xlabel]("xdata",fontsize=20)
# ax1[:set_ylabel]("ydata",fontsize=20)
ax2[:set_xlabel]("xdata",fontsize=20)
ax2[:set_ylabel]("ydata",fontsize=20)

# define the lines for the subplots. Their data will later be changed creating the animation
middle = ax10[:plot]([],[],"r-+")[1]
 car_plot = ax10[:plot]([],[],"b--")[1]
# line3 = ax2[:plot]([],[],"g-+")[1]
# constLine = ax2[:plot]([],[],"k--")[1]


# Function is called by animator to create the non-animated "background" of the animation
# important: return the line and a nothing at the end
function initPlot()
  #x-y Plot of the racetrack and obstacle and car
        #f_xy_plot= figure(3)
        #f_xy_plot[:canvas][:set_window_title]("Track and cars in XY plane")
        #ax10[:plot](x_track',y_track', linestyle="--", color = "#FF9900", linewidth = 0.5, label="_nolegend_")#plot the racetrack
        middle[:set_data](x_track',y_track')
        #plot older trajectory
        # colorObjectXY= colorModule.ColorManager()
        # for k= 1:n_oldTrajPlots
        #     colorXYold = colorModule.getColor(colorObjectXY)
        #     ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+k],1,j+k], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j+k],2,j+k],linestyle="--", color = colorXYold, label= "$k old Traj")
        # end
        
        # if interactive_plot == 1
        #     car_plot = ax10[:plot](oldTraj.oldTrajXY[1,1,j], oldTraj.oldTrajXY[1,2,j], color = "black") # just dummy to use remove func later
        #     ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],2,j], color = "black" ,linestyle=":", label="current Traj")# plot trajectory less visible to overwrite it with interactiv simulation
        #     pred_plot = ax10[:plot](xy_pred[:,1,1,j],xy_pred[:,2,1,j],color = "yellow", marker="o") 
        #     obstacle_plot = ax10[:plot](obstacle.xy_vector[1,1,j], obstacle.xy_vector[1,2,j], color = "red",marker="o", label = "obstacle Traj")
        #     y_obst_plot = ax10[:plot]([obstacle.axis_y_up[1,1,j],obstacle.axis_y_down[1,1,j]],[obstacle.axis_y_up[1,2,j],obstacle.axis_y_down[1,2,j]],color = "red")#plot the y semi axis
        #     s_obst_plot = ax10[:plot]([obstacle.axis_s_up[1,1,j],obstacle.axis_s_down[1,1,j]],[obstacle.axis_s_up[1,2,j],obstacle.axis_s_down[1,2,j]],color = "red")# plot the s semi axis
        # else
        #     car_plot = ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],2,j], color = "black", label="current Traj")
        #     ax10[:plot](obstacle.xy_vector[1:oldTraj.oldNIter[j],1,j], obstacle.xy_vector[1:oldTraj.oldNIter[j],2,j], color = "red")
        #     ax10[:plot](obstacle.xy_vector[oldTraj.oldNIter[j],1,j], obstacle.xy_vector[oldTraj.oldNIter[j],2,j], color = "red", marker = "o")
        #     y_obst_plot = ax10[:plot]([obstacle.axis_y_up[oldTraj.oldNIter[j],1,j],obstacle.axis_y_down[oldTraj.oldNIter[j],1,j]],[obstacle.axis_y_up[oldTraj.oldNIter[j],2,j],obstacle.axis_y_down[oldTraj.oldNIter[j],2,j]],color = "red")#plot the y semi axis
        #     s_obst_plot = ax10[:plot]([obstacle.axis_s_up[oldTraj.oldNIter[j],1,j],obstacle.axis_s_down[oldTraj.oldNIter[j],1,j]],[obstacle.axis_s_up[oldTraj.oldNIter[j],2,j],obstacle.axis_s_down[oldTraj.oldNIter[j],2,j]],color = "red")# plot the s semi axis
        # end
        # #plot the boundary lines
        # ax10[:plot](boundary_up[1,:], boundary_up[2,:],color="green",linestyle=":")
        # ax10[:plot](boundary_down[1,:], boundary_down[2,:],color="green",linestyle=":")
        # for l=1:20
        #     ax10[:plot]([boundary_down[1,l*50+1],boundary_up[1,l*50+1]],[boundary_down[2,l*50+1],boundary_up[2,l*50+1]], color = "black", linestyle = ":", linewidth = 0.5)
        #     text(boundary_down[1,l*50+1]+1,boundary_down[2,l*50+1],"s = $(l*5)",fontsize=8)
        # end
        # gca()[:set_aspect]("equal", adjustable="box")
        # grid() 
        # legend(bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    #end    
  return(middle, nothing)
end

# the function is called by the animator to change the plot elements that are returned. I basically change the
# x and y data of lines 1 to 3.
# important: return the lines + a nothing
# k is automatically incremented everytime the animate function fires. k=k+1 accounts for the fact that python starts with 0 where julia starts with 1
function animate(k)
  k = k+1

    car_plot[:set_data](oldTraj.oldTrajXY[1:k,1,j], oldTraj.oldTrajXY[1:k,2,j])
    # line2[:set_data]([xdata[1:k]],[zdata[1:k]])
    # line3[:set_data]([xdata[1:k]],[zdata[1:k]])

    return (car_plot,nothing)
end

# animator object runs the animation
# frames: defines how often animate function is called and interval how often in milliseconds
anim = animation.FuncAnimation(fig, animate, init_func=initPlot,frames=180, interval=50,repeat=false)

# to save video: you need encoder and writer, for ubuntu it should work if you do: sudo apt-get install ffmpeg
# FFwriter = animation.FFMpegWriter( fps=30, extra_args=["-vcodec", "libx264"])
# anim[:save]("basic_animation.mp4", writer = FFwriter)
