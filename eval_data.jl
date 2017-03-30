# This script is used to evaluate simulated data

using PyPlot
using JLD

function initPlot()
    linewidth = 0.5
    rc("axes", linewidth=linewidth)
    rc("lines", linewidth=linewidth, markersize=2)
    #rc("font", family="")
    rc("axes", titlesize="small", labelsize="small")        # can be set in Latex
    rc("xtick", labelsize="x-small")
    rc("xtick.major", width=linewidth/2)
    rc("ytick", labelsize="x-small")
    rc("ytick.major", width=linewidth/2)
    rc("text", usetex="true")
    rc("legend", fontsize="small")
    rc("font",family="serif")
    rc("font",size=10)
    rc("figure",figsize=[5.0,3])
    #rc("text.latex", preamble = """\\usepackage[utf8x]{inputenc}\\usepackage[T1]{fontenc}\\usepackage{lmodern}""")               # <- tricky! -- gotta actually tell tex to use!
    #rc("pgf", texsystem="pdflatex",preamble=L"""\usepackage[utf8x]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}""")
end

function smooth(x,n)
    y = zeros(size(x))
    for i=1:size(x,1)
        start = max(1,i-n)
        fin = min(size(x,1),start + 2*n)
        y[i,:] = mean(x[start:fin,:],1)
    end
    return y
end

initPlot()

include("barc_lib/classes.jl")
include("barc_lib/LMPC/functions.jl")

log_path_sim = "sim_dyn.jld"
d_sim = load(log_path_sim)

z    = d_sim["z"]
x    = d_sim["x"]
t    = d_sim["t"]
u    = d_sim["u"]
cost = d_sim["totalCost"]

z[z.==0] = NaN
x[x.==0] = NaN

# Create Track
s_track = 0.01:.01:50.87
c_track = zeros(5087)
c_track[201:400] = linspace(0,-pi/4,200)
c_track[401:600] = linspace(-pi/4,0,200)
c_track[701:900] = linspace(0,-pi/4,200)
c_track[901:1100] = linspace(-pi/4,0,200)
c_track[1101:1200] = linspace(0,-pi/4,100)
c_track[1201:1300] = linspace(-pi/4,0,100)
c_track[1601:2100] = linspace(0,10*pi/4/10,500)
c_track[2101:2600] = linspace(10*pi/4/10,0,500)
c_track[2601:2900] = linspace(0,-pi/3,300)
c_track[2901:3200] = linspace(-pi/3,0,300)
c_track[3501:3700] = linspace(0,-2*pi/2/4,200)
c_track[3701:3900] = linspace(-2*pi/2/4,0,200)
c_track[4041:4340] = linspace(0,-2*pi/2/6,300)
c_track[4341:4640] = linspace(-2*pi/2/6,0,300)

path_x,xl,xr = s_to_x(s_track,c_track)

# Plot and save e_Y
lapn = [1,10,20,30,40]      # specify which laps should be plotted
figure(1)
for i=1:size(lapn,1)
    plot(z[:,6,lapn[i]],z[:,5,lapn[i]],label="Lap $(lapn[i])")
end
lgd0=legend(loc="center left", bbox_to_anchor=(1, 0.5))
grid("on")
title("Comparison of lateral deviation")
xlabel("\$s\$ [m]")
ylabel("\$e_Y\$ [m]")
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Dyn_eY.pgf"
savefig(path_to_file, bbox_extra_artists=(lgd0,), bbox_inches="tight")


# Plot and save x,y
figure(2)
plot(path_x[:,1],path_x[:,2],"b--",xl[:,1],xl[:,2],"b-",xr[:,1],xr[:,2],"b-")
axis("equal")
arrow(0.3,-0.7,1,0,head_width=0.1)

for i=1:size(lapn,1)
    plot(x[:,1,lapn[i]],x[:,2,lapn[i]],label="Lap $(lapn[i])")
end
lgd4=legend(loc="center left", bbox_to_anchor=(1, 0.5))
grid("on")
title("x-y-view")
xlabel("\$x\$ [m]")
ylabel("\$y\$ [m]")
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Dyn_xy.pgf"
savefig(path_to_file, bbox_extra_artists=(lgd4,), bbox_inches="tight")

# Plot and save v
figure(3)
for i=1:size(lapn,1)
    subplot(2,1,1)
    plot(z[:,6,lapn[i]],z[:,1,lapn[i]],label="Lap $(lapn[i])")
    subplot(2,1,2)
    plot(z[:,6,lapn[i]],z[:,2,lapn[i]],label="Lap $(lapn[i])")
end
subplot(2,1,1)
lgd1=legend(loc="center left", bbox_to_anchor=(1, 0.5))
grid("on")
xlabel("\$s\$ [m]")
ylabel("\$v_x\$ \$\\left[\\frac{m}{s}\\right]\$")
subplot(2,1,2)
lgd2=legend(loc="center left", bbox_to_anchor=(1, 0.5))
grid("on")
xlabel("\$s\$ [m]")
ylabel("\$v_y\$ \$\\left[\\frac{m}{s}\\right]\$")
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Dyn_v.pgf"
savefig(path_to_file, bbox_extra_artists=(lgd1,lgd2), bbox_inches="tight")

# Plot and save cost
figure(4)
plot(cost[1:40]/50,"-o")
grid("on")
title("Total cost")
xlabel("Iteration")
ylabel("Iteration cost (\$t\$ [\$s\$])")
ylim([0,maximum(cost)/50*1.1])
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Dyn_cost.pgf"
savefig(path_to_file)

lapn = [2,40]
# Plot friction circle/accelerations
figure()
pstyle = ("--","-")
for i=1:size(lapn,1)
    a_x = diff(z[:,1,lapn[i]])./diff(t)
    a_y = diff(z[:,2,lapn[i]])./diff(t)
    #plot(a_x,a_y)
    plot(smooth(a_x,2),smooth(a_y,1),pstyle[i],label="Lap $(lapn[i])")
end
grid("on")
legend()
axis("equal")
xlabel("\$a_x \\ \\left[\\frac{m}{s^2}\\right]\$")
ylabel("\$a_y \\ \\left[\\frac{m}{s^2}\\right]\$")
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Dyn_fcircle.pgf"
savefig(path_to_file)

# Plot velocity over 2D-track
i = 40
v = sqrt(z[:,1,i].^2+z[:,2,i].^2)
figure()
plot(path_x[:,1],path_x[:,2],"b--",xl[:,1],xl[:,2],"b-",xr[:,1],xr[:,2],"b-")
scatter(x[:,1,i],x[:,2,i],c=v,cmap=ColorMap("jet"),edgecolors="face",vmin=minimum(v),vmax=maximum(v))
grid("on")
axis("equal")
xlabel("x [m]")
ylabel("y [m]")
title("x-y-view in iteration $i")
cb = colorbar()
cb[:set_label]("Velocity \$\\left[\\frac{m}{s}\\right]\$")
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Dyn_v_over_xy.pgf"
savefig(path_to_file)

it_cost = zeros(40)
for i=1:40
    it_cost[i] = size(z[z[:,6,i].>0,:,i],1)
end
figure()
plot(it_cost/50,"-x")
plot(cost/50/2,"-")
legend(["Single","Double"])
grid("on")
xlabel("Iteration")
ylabel("Iteration cost")
xlim([0,20])
ylim([4,10])
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Dyn_double_it_cost.pgf"
savefig(path_to_file)

