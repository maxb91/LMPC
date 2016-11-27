using PyPlot
using JLD

include("barc_lib/classes.jl")
include("barc_lib/LMPC/functions.jl")

log_path_sim = "sim.jld"
d_sim = load(log_path_sim)

z    = d_sim["z"]
x    = d_sim["x"]
t    = d_sim["t"]
u    = d_sim["u"]

z[z.==0] = NaN
x[x.==0] = NaN

# Create Track
s_track = 0.01:.01:50.49
c_track = zeros(5049)
c_track[1:300] = 0
c_track[301:400] = linspace(0,-pi/2,100)
c_track[401:500] = linspace(-pi/2,0,100)
c_track[501:900] = 0
c_track[901:1000] = linspace(0,-pi/2,100)
c_track[1001:1100] = linspace(-pi/2,0,100)
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
c_track[3901:4102] = 0
c_track[4103:4402] = linspace(0,-2*pi/2/6,300)
c_track[4403:4702] = linspace(-2*pi/2/6,0,300)

figure(1)
for i=1:15
    plot(z[:,6,i],z[:,1:5,i])
    #readline()
end
grid("on")
title("Comparison of 15 laps")
xlabel("s [m]")



path_x,xl,xr = s_to_x(s_track,c_track)

figure(2)
plot(path_x[:,1],path_x[:,2],"b--",xl[:,1],xl[:,2],"b-",xr[:,1],xr[:,2],"b-")
axis("equal")
for i=1:14:15
    plot(x[:,1,i],x[:,2,i])
    readline()
end
grid("on")
title("Comparison of 15 laps")
xlabel("s [m]")

