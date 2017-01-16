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
    rc("figure",figsize=[5.5,3.0])
    #rc("text.latex", preamble = """\\usepackage[utf8x]{inputenc}\\usepackage[T1]{fontenc}\\usepackage{lmodern}""")               # <- tricky! -- gotta actually tell tex to use!
    #rc("pgf", texsystem="pdflatex",preamble=L"""\usepackage[utf8x]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}""")
end

initPlot()

include("barc_lib/classes.jl")
include("barc_lib/LMPC/functions.jl")

log_path_sim = "sim_kin.jld"
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
c_track[1:200] = 0
c_track[201:400] = linspace(0,-pi/4,200)
c_track[401:600] = linspace(-pi/4,0,200)
c_track[601:700] = 0
c_track[701:900] = linspace(0,-pi/4,200)
c_track[901:1100] = linspace(-pi/4,0,200)
c_track[1101:1200] = linspace(0,-pi/4,100)
c_track[1201:1300] = linspace(-pi/4,0,100)
c_track[1301:1600] = 0
c_track[1601:2100] = linspace(0,10*pi/4/10,500)
c_track[2101:2600] = linspace(10*pi/4/10,0,500)
c_track[2601:2900] = linspace(0,-pi/3,300)
c_track[2901:3200] = linspace(-pi/3,0,300)
c_track[3201:3500] = 0
c_track[3501:3700] = linspace(0,-2*pi/2/4,200)
c_track[3701:3900] = linspace(-2*pi/2/4,0,200)
c_track[3901:3902] = 0
c_track[4041:4340] = linspace(0,-2*pi/2/6,300)
c_track[4341:4640] = linspace(-2*pi/2/6,0,300)

path_x,xl,xr = s_to_x(s_track,c_track)
#xl[end,:] = xl[1,:]
#xr[end,:] = xr[1,:]
# Plot and save e_Y
lapn = [1,2,10]      # specify which laps should be plotted
figure(1)
for i=1:size(lapn,1)
    plot(z[:,1,lapn[i]],z[:,2,lapn[i]],label="Lap $(lapn[i])")
    #readline()
end
legend()
grid("on")
title("Comparison of lateral deviation")
xlabel("\$s\$ [m]")
ylabel("\$e_Y\$ [m]")
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Kin_eY.pgf"
savefig(path_to_file)


# Plot and save x,y
figure(2)
plot(path_x[:,1],path_x[:,2],"b--",xl[:,1],xl[:,2],"b-",xr[:,1],xr[:,2],"b-")
axis("equal")
for i=1:size(lapn,1)
    plot(x[:,1,lapn[i]],x[:,2,lapn[i]],label="Lap $(lapn[i])")
    #readline()
end
legend()
grid("on")
title("x-y-view")
xlabel("\$x\$ [m]")
ylabel("\$y\$ [m]")
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Kin_xy.pgf"
savefig(path_to_file)

# Plot and save v
figure(3)
for i=1:size(lapn,1)
    plot(z[:,1,lapn[i]],z[:,4,lapn[i]],label="Lap $(lapn[i])")
    #readline()
end
legend()
grid("on")
title("Velocities")
xlabel("\$s\$ [m]")
ylabel("\$v\$ \$\\left[\\frac{m}{s}\\right]\$")
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Kin_v.pgf"
savefig(path_to_file)

# Plot and save cost
figure(4)
plot(cost[1:10]/10,"-o")
grid("on")
title("Total cost")
xlabel("Iteration")
ylabel("Iteration cost")
ylim([0,60])
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Kin_cost.pgf"
savefig(path_to_file)

figure(5)
plot(s_track, c_track)
grid("on")
xlabel("\$s\$ [\$ m\$]")
ylabel("\$ c(s) \$ \$\\left[ \\frac{1}{m} \\right] \$")
tight_layout()
path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/6. Final report/Figures/Simulation/Kin_curv.pgf"
savefig(path_to_file)


