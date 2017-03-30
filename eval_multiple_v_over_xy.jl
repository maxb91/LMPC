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
    rc("figure",figsize=[6,4])
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

log_path_sim = "sim_dyn_good.jld"
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

# Plot velocity over 2D-track
for i=1:40
    # i = 40
    v = sqrt(z[:,1,i].^2+z[:,2,i].^2)
    figure()
    plot(path_x[:,1],path_x[:,2],"b--",xl[:,1],xl[:,2],"b-",xr[:,1],xr[:,2],"b-")
    #scatter(x[:,1,i],x[:,2,i],c=v,cmap=ColorMap("jet"),edgecolors="face",vmin=minimum(v),vmax=maximum(v))
    scatter(x[:,1,i],x[:,2,i],c=v,cmap=ColorMap("jet"),edgecolors="face",s=1,vmin=2,vmax=3.5)
    xlim([-9,5])
    grid("on")
    #axis("equal")
    gca()[:set_aspect]("equal")
    xlabel("x [m]")
    ylabel("y [m]")
    title("x-y-view in iteration $i")
    cb = colorbar()
    cb[:set_label]("Velocity \$\\left[\\frac{m}{s}\\right]\$")
    tight_layout()
    path_to_file = "/Users/Maximilian/Documents/ETH/Master/Master thesis/7. Presentation/Zurich/Figures/sim/lap$(i).pdf"
    savefig(path_to_file)
end
