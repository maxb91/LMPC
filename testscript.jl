using PyCall
using PyPlot
using JLD
close("all")
@pyimport matplotlib as mpl
@pyimport matplotlib.patches as patches
@pyimport matplotlib.animation as animation#
fig = figure()
a = patches.Rectangle([0.2,0.0],6.0,5.0)
ax5 = gca()
ax5[:add_patch](a)

# function initPlot()

#    ax5[:add_patch](a)    
#    # return a
# end


function animate(frame)
    frame+=1
    a[:set_xy]([0.0,0.01*frame])

    return a
end



anim = animation.FuncAnimation(fig, animate, frames=200,interval=100)# init_func=initPlot, )#,repeat=false)#, blit = true )
# to save video: you need encoder and writer, for ubuntu it should work if you do: sudo apt-get install ffmpeg
# FFwriter = animation.FFMpegWriter( fps=8, bitrate=1500, extra_args=["-vcodec", "libx264"])   #(fps=10,bitrate=3000,extra_args=["-vcodec", "libx264","-pix_fmt","yuv420p"])
# anim[:save]("basic_animation.mp4", writer = FFwriter)


# fig = figure()
# ax1 = subplot(1,1,1)
# drawPatch=drawCar(ax,[0.0,0.0,0.0], "green")

# for p in drawPatch
#                 ax1[:add_patch](p)
# end

# updateCarParts(ax,patches,pos)
