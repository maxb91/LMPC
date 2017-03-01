
using PyPlot
using JLD
using PyCall
@pyimport matplotlib as mpl
@pyimport matplotlib.patches as patches
matplotlib[:style][:use]("classic")
include("../helper/classes.jl")
include("plot_functions.jl")
include("calculateObstacleXY.jl")
include("colorModule.jl")

close("all")
# file = "LMPC/trackdata.jld"
# file = "tracks/oval.jld"
# file = "LMPC/small_oval.jld"
# file = "LMPC/adv_track.jld"
# file = "LMPC/optm_adv.jld"
# file = "LMPC/sim_curv.jld"
file = "data/2017-02-10-16-00-Data.jld"

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


ds = trackCoeff.ds
dt = modelParams.dt
xy_track  = [x_track; y_track]
t   = collect(0:dt:(buffersize-1)*dt)

xy_track  = [x_track; y_track]
trackL = size(xy_track,2)
boundary_up = zeros(2,trackL)
boundary_down = zeros(2,trackL)
boundary_color = "black"
ego_color = "green"
obstacle_color = "white"
width=0.6
width = width+0.2
for kkk = 1:1:trackL
    if 1< kkk < trackL 
        xt_secant = x_track[kkk+1] - x_track[kkk-1]
        yt_secant = y_track[kkk+1] - y_track[kkk-1]
        normVec = [-yt_secant; xt_secant]/norm([-yt_secant; xt_secant])
        boundary_up[:,kkk] = xy_track[:,kkk] + normVec * width /2
        boundary_down[:,kkk] = xy_track[:,kkk] - normVec * width /2
    elseif kkk == 1
        xt_secant = x_track[kkk+1] - x_track[kkk]
        yt_secant = y_track[kkk+1] - y_track[kkk]
        normVec = [-yt_secant; xt_secant]/norm([-yt_secant; xt_secant])
        boundary_up[:,kkk] = xy_track[:,kkk] + normVec * width /2
        boundary_down[:,kkk] = xy_track[:,kkk] - normVec * width /2
    else
        xt_secant = x_track[kkk] - x_track[kkk-1]
        yt_secant = y_track[kkk] - y_track[kkk-1]
        normVec = [-yt_secant; xt_secant]/norm([-yt_secant; xt_secant])
        boundary_up[:,kkk] = xy_track[:,kkk] + normVec * width /2
        boundary_down[:,kkk] = xy_track[:,kkk] - normVec * width /2
    end
end


# xy_pred = zeros(mpcParams.N+1,2,size(t)[1],oldTraj.n_oldTraj)
# obstOrientation = zeros(size(obstacle.xy_vector)[1],size(obstacle.xy_vector)[3],size(obstacle.xy_vector)[4])
# for k = 1:oldTraj.n_oldTraj
#     for i = 1:oldTraj.oldNIter[k]
#         # caluclate the predicted XY postion of the car from the s-ey values
#         xy_pred[:,:,i,k] = calculatePredictedXY(oldTraj.z_pred_sol[:,:,:,k], mpcParams, trackCoeff, xy_track, convert(Int64,i),k)
#         #calculate the obstacle postion from the s-ey values
#         obstOrientation = calculateObstacleXY!(obstacle, trackCoeff, xy_track,i,k,obstOrientation ) #this funciton computes values for row i
#     end
# end

f_xy_plot= figure(1)
f_xy_plot[:canvas][:set_window_title]("Track and cars in XY plane")
ax10 = subplot(1,1,1)
# plot the boundary lines
ax10[:plot](x_track',y_track', linestyle = (0, (4.0, 8.0)), color =boundary_color, linewidth = 0.4, label="_nolegend_")#plot the racetrack
ax10[:plot](boundary_up[1,:], boundary_up[2,:],color=boundary_color, linewidth = 0.7)#,linestyle=":")
ax10[:plot](boundary_down[1,:], boundary_down[2,:],color=boundary_color, linewidth = 0.7)#,linestyle=":")


ax10[:plot]([boundary_down[1,1],boundary_up[1,1]],[boundary_down[2,1],boundary_up[2,1]], color = "black", linewidth = 2)
ax10[:plot]([boundary_down[1,trackL],boundary_up[1,trackL]],[boundary_down[2,trackL],boundary_up[2,trackL]], color = "black", linewidth = 2)
for l=1:convert(Int64,trunc(trackL/51))
    ax10[:plot]([boundary_down[1,l*50+1],boundary_up[1,l*50+1]],[boundary_down[2,l*50+1],boundary_up[2,l*50+1]], color = "black", linestyle = ":", linewidth = 0.5)
    boundvec = [boundary_up[1,l*50+1]-boundary_down[1,l*50+1];boundary_up[2,l*50+1]-boundary_down[2,l*50+1]]
    # ax10[:text](boundary_down[1,l*50+1]+1.85*boundvec[1],boundary_down[2,l*50+1]+1.85*boundvec[2],"$(convert(Int64,l*50*ds))",fontsize=22,clip_on = true)
end

# limits fuer sim_curv
# ax10[:set_xlim](-8,11)
# ax10[:set_ylim](-11.5,3)
# ax10[:set_xticks]([ -5, 0, 5, 10])
# ax10[:set_yticks]([-10, -5, 0])

# limits fuer optm_adv
# ax10[:set_xlim](-10,6)
# ax10[:set_ylim](-11.5,3)
# ax10[:set_xticks]([ -5, 0, 5])
# ax10[:set_yticks]([-10, -5, 0])

# limits fuer oval
# ax10[:set_xlim](-15,15)
# ax10[:set_ylim](-11.5,3)
# ax10[:set_xticks]([-15, -10 ,-5, 0, 5, 10, 15])
# ax10[:set_yticks]([-10, -5, 0])

ax10[:tick_params](axis="both", which="major", labelsize=24)
# ax10[:set_xlabel](L"X \>\>\>[m] " , fontsize = 24)
# ax10[:set_ylabel](L"Y \>\>\> [m] ", fontsize = 24)
ax10[:set_aspect]("equal", adjustable="box")
ax10[:tick_params]( axis="x", which="both", bottom="off", top="off", labelbottom="off")
ax10[:tick_params]( axis="y", which="both", left="off", right="off", labelleft="off")
# ax10[:grid]()


j = 1
# car_plot = ax10[:plot](oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],2,j], color = ego_color, label="current Traj")
cm = mpl.cm[:get_cmap]("bwr")#("RdYlBu")
car_plot = scatter(oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],1,j], oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],2,j], c =  sqrt(oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],3,j].^2 +oldTraj.oldTrajXY[1:oldTraj.oldNIter[j],4,j].^2),vmin=1.3, vmax=2.1,s= 80)#, cmap = cm)
cbar = f_xy_plot[:colorbar](car_plot,ax=ax10, shrink=1.0, pad=.01)#, aspect=10) shrink = 0.77 for oval
cbar[:ax][:tick_params](labelsize=20)
cbar[:set_label]("velocity [m/s]", rotation=90, fontsize =20)


# obst_patch =Array{PyCall.PyObject}(obstacle.n_obstacle)
# plt_obst =Array{PyCall.PyObject}(obstacle.n_obstacle)
# for ii = 1:obstacle.n_obstacle
#     obst_patch[ii] = patches.Ellipse([obstacle.xy_vector[1,1,j,ii],obstacle.xy_vector[1,2,j,ii]],2*obstacle.rs,2*obstacle.ry,color=obstacle_color,alpha=1.0,fill=false, angle = obstOrientation[1,j,ii]*180/pi, linewidth=2)
#     plt_obst[ii] = ax10[:add_patch](obst_patch[ii])
# end