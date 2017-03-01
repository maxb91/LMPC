function calc_track_boundaryXY(xy_track::Array{Float64}, trackCoeff)
    trackL = size(xy_track,2)
    boundary_up = zeros(2,trackL)
    boundary_down = zeros(2,trackL)
    trackCoeff.width = trackCoeff.width+0.2 # for plotting add real boundary
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
    return boundary_up, boundary_down, trackL
end


function drawCar(ax,pos::Array{Float64},mycolor::String="black")

    x = pos[1]
    y = pos[2]
    psi = pos[3]

    pos = [x;y]

    # define car geometry
    length = 40 / 100.0
    width = 18.0 / 100.0
    btmLeft_vortex = pos + [-length/2,-width/2]

    # define window geometry
    window_length = 10.0/100.0
    window_width = 10.0 / 100.0
    window_distance = 6.0 / 100.0 #from c.g. 
    windowLeft_vortex = pos + [window_distance,-window_width/2]

    # define tire geometry
    tire_width = 1.0/100.0
    tire_length = 5.0/100.0
    tire_distance = 8.0/100.0

    # define non-rotated patches that make up the car's top-view
    car = patches.Rectangle(btmLeft_vortex,length,width,facecolor=mycolor, linewidth=3.0,edgecolor="black",zorder=4)
    window = patches.Rectangle(windowLeft_vortex,window_length,window_width,fill=false, linewidth=3.0,edgecolor="black",zorder=4)
    tl = patches.Rectangle(pos + [-tire_distance-tire_length,width/2],tire_length,tire_width,color="black",zorder=4)
    tr = patches.Rectangle(pos + [tire_distance,width/2],tire_length,tire_width,color="black",zorder=4)
    bl = patches.Rectangle(pos + [-tire_distance-tire_length,-width/2-tire_width],tire_length,tire_width,color="black",zorder=4)
    br = patches.Rectangle(pos + [tire_distance,-width/2-tire_width],tire_length,tire_width,color="black",zorder=4)

    car_patches = [car,window,tl,tr,bl,br]
    # define a rotation transformation and a transformation from data c.s. to display c.s.
    t1 = mpl.transforms[:Affine2D]()
    t1[:rotate_deg_around](pos[1],pos[2], psi)
    t2 = ax[:transData]
    # combine transformations
    t3 = t1[:__add__](t2)

    # apply rotation transformation and add rotated patch to axis
    for p in car_patches
        p[:set_transform](t3)
    end

    return car_patches    
end

function updateCarParts(ax,patches,pos::Array{Float64})
    # [car,window,tl,tr,bl,br]

    car = patches[1]
    window = patches[2]
    tl = patches[3]
    tr = patches[4]
    bl = patches[5]
    br = patches[6]

    x = pos[1]
    y = pos[2]
    psi = pos[3]
    pos = [x;y]


    # define car geometry
    length = 36.5 / 100.0
    width = 18.0 / 100.0
    btmLeft_vortex = pos + [-length/2,-width/2]

    # define window geometry
    window_length = 10.0/100.0
    window_width = 10.0 / 100.0
    window_distance = 6.0 / 100.0 #from c.g. 
    windowLeft_vortex = pos + [window_distance,-window_width/2]

    # define tire geometry
    tire_width = 1.0/100.0
    tire_length = 5.0/100.0
    tire_distance = 8.0/100.0

    # set 
    car[:set_xy](btmLeft_vortex)
    window[:set_xy](windowLeft_vortex)
    tl[:set_xy](pos + [-tire_distance-tire_length,width/2])
    tr[:set_xy](pos + [tire_distance,width/2])
    bl[:set_xy](pos + [-tire_distance-tire_length,-width/2-tire_width])
    br[:set_xy](pos + [tire_distance,-width/2-tire_width])
    # define a rotation transformation and a transformation from data c.s. to display c.s.
    t1 = mpl.transforms[:Affine2D]()
    t1[:rotate_deg_around](x,y, psi)
    t2 = ax[:transData]
    # combine transformations
    t3 = t1[:__add__](t2)

    # apply rotation transformation and add rotated patch to axis
    for p in patches
        p[:set_transform](t3)
    end
    return patches
end






function plotfct_states_over_t(oldTraj, t, j)
    fig_1 = figure(1)
    fig_1[:canvas][:set_window_title]("States and Inputs over t")
    ax1=subplot(211)
    plot(t[1:oldTraj.oldNIter[j]],oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],"y",t[1:oldTraj.oldNIter[j]],oldTraj.oldTraj[1:oldTraj.oldNIter[j],2,j],"r",t[1:oldTraj.oldNIter[j]],oldTraj.oldTraj[1:oldTraj.oldNIter[j],3,j],"g",t[1:oldTraj.oldNIter[j]],oldTraj.oldTraj[1:oldTraj.oldNIter[j],4,j],"b")
    grid(1)
    legend(["s","eY","ePsi","v"], bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    title("States")
    ax2=subplot(212,sharex=ax1)
    plot(t[1:oldTraj.oldNIter[j]-1],oldTraj.oldInput[1:oldTraj.oldNIter[j]-1,1,j],"r",t[1:oldTraj.oldNIter[j]-1],oldTraj.oldInput[1:oldTraj.oldNIter[j]-1,2,j],"g")
    grid(1)
    title("Control input")
    legend(["a","d_f"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    # ax3=subplot(313,sharex=ax1)
    # plot(t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,1,j],"r",t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,3,j],"b",t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,4,j],"y",t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,5,j],"m",t[1:oldTraj.oldNIter[j]-1],oldTraj.costs[1:oldTraj.oldNIter[j]-1,6,j],"c", t[1:oldTraj.oldNIter[j]-1], oldTraj.costs[1:oldTraj.oldNIter[j]-1,7,j])
    # grid(1)
    # title("Cost distribution")
    # legend(["z","z_Term_const","deriv","control","lane", "Obstacle"], bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    nothing
end







function plotfct_costs(oldTraj,j)
    plt_cost = figure(2)
    plt_cost[:canvas][:set_window_title]("Costs over s")
    ax9 = subplot(3,2,1)
    ax9[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[2,1:oldTraj.oldNIter[j]-1,j])  
    grid()
    xlabel("s in [m]")
    ylabel("Terminal cost ") 

    ax7= subplot(3,2,2,sharex=ax9)
    ax7[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[3,1:oldTraj.oldNIter[j]-1,j])  
    grid()
    xlabel("s in [m]")
    ylabel("cost constraint")    
    ax8= subplot(3,2,3,sharex=ax9)
    ax8[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[8,1:oldTraj.oldNIter[j]-1,j])  
    grid()
    xlabel("s in [m]")
    ylabel("cost Obstacle")     
    ax12= subplot(3,2,4,sharex=ax9)
    ax12[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[4,1:oldTraj.oldNIter[j]-1,j])  
    grid()
    xlabel("s in [m]")
    ylabel("deriv cost")    
    ax13= subplot(3,2,5,sharex=ax9)
    ax13[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[6,1:oldTraj.oldNIter[j]-1,j])  
    grid()
    xlabel("s in [m]")
    ylabel("lane Cost")     
    ax14= subplot(3,2,6,sharex=ax9)
    ax14[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.costs[1,1:oldTraj.oldNIter[j]-1,j])  
    ax14[:grid]()
    xlabel("s in [m]")
    ylabel("z Cost")
    return ax9
end









function plotfct_lambda(oldTraj,j)
    f_lambda =figure(4)
    f_lambda[:canvas][:set_window_title]("Lambda and ssOn values over t")
    ax11= subplot(2,2,1)
    colorObjectLambda = colorModule.ColorManager()
    for k=1: oldTraj.n_oldTraj//4
        k = convert(Int64,k)
        colorLambda = colorModule.getColor(colorObjectLambda)
        scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.lambda_sol[k,1:oldTraj.oldNIter[j]-1,j].*oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorLambda, marker = "x", label = string(L"\lambda","k"))
    end
    xlabel("s in [m]")
    ylabel("lambda")
    grid()
    #ax11[:set_ylim]([-0.1,1.1])
    legend()
    #legend([L"\lambda 1",L"\lambda 2",L"\lambda 3",L"\lambda 4",L"\lambda 5", L"\lambda 6",L"\lambda 7",L"\lambda 8",L"\lambda 9",L"\lambda 10"],bbox_to_anchor=(-0.201, 1), loc=2, borderaxespad=0.)
    
    axlambda2= subplot(2,2,2,sharex= ax11, sharey =ax11)
    for k=oldTraj.n_oldTraj//4+1: oldTraj.n_oldTraj//2
        k = convert(Int64,k)
        colorLambda = colorModule.getColor(colorObjectLambda)
        scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.lambda_sol[k,1:oldTraj.oldNIter[j]-1,j].*oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorLambda, marker = "x")
    end
    xlabel("s in [m]")
    ylabel("lambda")
    grid()
    legend([L"\lambda /4+1",L"\lambda /4+2",L"\lambda 3",L"\lambda 4",L"\lambda 5", L"\lambda 6"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)

    axlambda3= subplot(2,2,3,sharex= ax11, sharey =ax11)
    for k=oldTraj.n_oldTraj//2+1: oldTraj.n_oldTraj//4*3
        k = convert(Int64,k)
        colorLambda = colorModule.getColor(colorObjectLambda)
        scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.lambda_sol[k,1:oldTraj.oldNIter[j]-1,j].*oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorLambda, marker = "x")
    end
    xlabel("s in [m]")
    ylabel("lambda")
    grid()
    legend([L"\lambda /2+1",L"\lambda /2+2",L"\lambda 3",L"\lambda ",L"\lambda 5"],bbox_to_anchor=(-0.201, 1), loc=2, borderaxespad=0.)

    axlambda4= subplot(2,2,4,sharex= ax11, sharey =ax11)
    for k=oldTraj.n_oldTraj*3//4+1: oldTraj.n_oldTraj
        k = convert(Int64,k)
        colorLambda = colorModule.getColor(colorObjectLambda)
        scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.lambda_sol[k,1:oldTraj.oldNIter[j]-1,j].*oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorLambda, marker = "x")
    end
    xlabel("s in [m]")
    ylabel("lambda")
    grid()
    legend([L"\lambda *3/4+1",L"\lambda *3/4+2",L"\lambda 3",L"\lambda 4",L"\lambda 5"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    return ax11
end
    







function plotfct_ssOn(oldTraj,j, ax11)   
    f_ssOn =figure(10)
    f_ssOn[:canvas][:set_window_title](" ssOn values over t")
    axssOn = subplot(2,2,1,sharex= ax11)
    colorObjectSafeSet= colorModule.ColorManager()
    for k=1: oldTraj.n_oldTraj//4
        k = convert(Int64,k)
        colorSafeSet= colorModule.getColor(colorObjectSafeSet)
        scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j], color = colorSafeSet, marker = "x")
    end
    xlabel("s in [m]")
    ylabel("ssOn")
    grid()
    legend(["ssOn1","ssOn2","ssOn3","ssOn4","ssOn5", "ssOn6","ssOn7","ssOn8","ssOn9","ssOn10"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    
    axssOn2 = subplot(2,2,2,sharex= ax11)
    for k=oldTraj.n_oldTraj//4+1: oldTraj.n_oldTraj//2
        k = convert(Int64,k)
        colorSafeSet= colorModule.getColor(colorObjectSafeSet)
        scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j],color = colorSafeSet, marker = "x")
    end
    xlabel("s in [m]")
    ylabel("ssOn")
    grid()
    legend(["ssOn/4+1","ssOn/4+2","ssOn3","ssOn4","ssOn5", "ssOn6","ssOn7","ssOn8","ssOn9","ssOn10"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    axssOn3 = subplot(2,2,3,sharex= ax11)
    for k=oldTraj.n_oldTraj//2+1: oldTraj.n_oldTraj*3//4
        k = convert(Int64,k)
        colorSafeSet= colorModule.getColor(colorObjectSafeSet)
        scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j],color = colorSafeSet, marker = "x")
    end
    xlabel("s in [m]")
    ylabel("ssOn")
    grid()
    legend(["ssOn/2+1","ssOn/2+2","ssOn3","ssOn4","ssOn5", "ssOn6","ssOn7","ssOn8","ssOn9","ssOn10"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
    axssOn4 = subplot(2,2,4,sharex= ax11)
    for k=oldTraj.n_oldTraj*3//4+1: oldTraj.n_oldTraj
        k = convert(Int64,k)
        colorSafeSet= colorModule.getColor(colorObjectSafeSet)
        scatter(oldTraj.oldTraj[1:oldTraj.oldNIter[j]-1,1,j],oldTraj.ssInfOn_sol[k,1:oldTraj.oldNIter[j]-1,j],color = colorSafeSet, marker = "x")
    end
    xlabel("s in [m]")
    ylabel("ssOn")
    grid()
    legend(["ssOn*3/4+1","ssOn/2+2","ssOn3","ssOn4","ssOn5", "ssOn6","ssOn7","ssOn8","ssOn9","ssOn10"],bbox_to_anchor=(1.001, 1), loc=2, borderaxespad=0.)
end










function plotfct_epsilon(oldTraj, j)
    f_eps = figure(8)
    f_eps[:canvas][:set_window_title]("values of Eps over s") 

    plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],oldTraj.eps[1,1:oldTraj.oldNIter[j],j])
    plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],oldTraj.eps[2,1:oldTraj.oldNIter[j],j])
    plot(oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],oldTraj.eps[3,1:oldTraj.oldNIter[j],j])
    xlabel("s in [m]")
    ylabel("value of epsilons")
    legend(["epsilon left boundary","epsilon right boundary", "epsilon velocity"])
    grid()
end








function plotfct_curvature(oldTraj,j)
    set_fontsize = 30
    f_curv_app = figure(7)
    f_curv_app[:canvas][:set_window_title]("Curvature approximation over s")  
    axCurv = subplot(1,1,1)
    axCurv[:plot](oldTraj.oldTraj[1:oldTraj.oldNIter[j],1,j],oldTraj.curvature[1:oldTraj.oldNIter[j],j])
    axCurv[:tick_params](axis="both", which="major", labelsize=set_fontsize)
    axCurv[:set_xlabel](L"s \>\>\>[m] " , fontsize = set_fontsize)
    axCurv[:set_ylabel](L"\kappa \>\>\> [1/m] ", fontsize = set_fontsize)
    # axCurv[:set_xticks]([-15, -10 ,-5, 0, 5, 10, 15])
    # axCurv[:set_yticks]([-0.3, -0.2,-0.1, 0])
    axCurv[:set_xlim](-0,52.7)
    # axCurv[:set_ylim](-1.3,1.2)
    axCurv[:grid]()
end







function plotfct_copied(oldTraj,j)
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
    f_copied_plot= figure(9)
    f_copied_plot[:canvas][:set_window_title]("Copied s over current s")
    axCopied = subplot(1,1,1)
    axCopied[:clear]()

    for i=1:oldTraj.oldNIter[j]-1
        if oldTraj.copyInfo[i,1,j]>0.0 # if lambda of copied traj is greater 0.1 -> if traj s used for solving.
            scatter(oldTraj.copyInfo[i,3,j],oldTraj.copyInfo[i,2,j], color = "#DCDCDC",s=100)
        end
    end
    for i=1:oldTraj.oldNIter[j]-1
        if oldTraj.copyInfo[i,4,j]>0.6 # if lambda of copied traj is greater 0.1 -> if traj s used for solving.
            # scatter(oldTraj.copyInfo[i,3,j],oldTraj.copyInfo[i,2,j], color = colordefs[convert(Int64,oldTraj.copyInfo[i,1,j])])
            scatter(oldTraj.copyInfo[i,3,j],oldTraj.copyInfo[i,2,j], color = "blue", s=100)
        end
    end
    axCopied[:set_aspect]("equal", adjustable="box")
    # axCopied[:set_ylim](0,20)
    # axCopied[:set_xlim](22,48)
    axCopied[:set_xlabel](L"current \> s \> [m]", fontsize =30)
    axCopied[:grid]()
    axCopied[:set_ylabel](L"old \> s \>  [m]", fontsize =30)
    axCopied[:tick_params]( axis="both", which="major", labelsize=30)
end