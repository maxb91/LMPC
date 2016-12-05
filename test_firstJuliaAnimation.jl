# Minimal-Example on how to do Julia animations with subplots and some static content using matplotlib

using PyCall
using PyPlot
@pyimport matplotlib.animation as animation# First set up the figure, the axis, and the plot element we want to animate
#@pyimport matplotlib.pyplot as plt
plt[:rcParams]["animation.ffmpeg_path"] ="C:/ffmpeg-3.2-win64-static/bin/"
# create dummy data
xdata = 0:1:300
ydata = 2*xdata
zdata = 3*xdata

# define axis limits
xmin = 0
xmax = 300
ymin = 0
ymax = 600

# create figure to plot
fig = figure(figsize=(15,8))

# add two subplots to the figure and define axis limits
ax1 = fig[:add_subplot](1, 2, 1)
ax1[:set_xlim](xmin, xmax)
ax1[:set_ylim](ymin, ymax)
ax2 = fig[:add_subplot](1,2,2)
ax2[:set_xlim](xmin, xmax)
ax2[:set_ylim](ymin, ymax)

# set labels
ax1[:set_xlabel]("xdata",fontsize=20)
ax1[:set_ylabel]("ydata",fontsize=20)
ax2[:set_xlabel]("xdata",fontsize=20)
ax2[:set_ylabel]("ydata",fontsize=20)

# define the lines for the subplots. Their data will later be changed creating the animation
line1 = ax1[:plot]([],[],"r-+")[1]
line2 = ax1[:plot]([],[],"b--")[1]
line3 = ax2[:plot]([],[],"g-+")[1]
constLine = ax2[:plot]([],[],"k--")[1]


# Function is called by animator to create the non-animated "background" of the animation
# important: return the line and a nothing at the end
function initPlot()
  constLine[:set_data]([xdata],[ydata])
  return(constLine, nothing)
end

# the function is called by the animator to change the plot elements that are returned. I basically change the
# x and y data of lines 1 to 3.
# important: return the lines + a nothing
# k is automatically incremented everytime the animate function fires. k=k+1 accounts for the fact that python starts with 0 where julia starts with 1
function animate(k)
  k = k+1

    line1[:set_data]([xdata[1:k]],[ydata[1:k]])
    line2[:set_data]([xdata[1:k]],[zdata[1:k]])
    line3[:set_data]([xdata[1:k]],[zdata[1:k]])

    return (line1,line2,line3,nothing)
end

# animator object runs the animation
# frames: defines how often animate function is called and interval how often in milliseconds
anim = animation.FuncAnimation(fig, animate, init_func=initPlot,frames=100, interval=50,repeat=false)

# to save video: you need encoder and writer, for ubuntu it should work if you do: sudo apt-get install ffmpeg
FFwriter = animation.FFMpegWriter( fps=30, extra_args=["-vcodec", "libx264"])
anim[:save]("basic_animation.mp4", writer = FFwriter)
