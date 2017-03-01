# include("PrintFig.jl")
# using PyPlot, PrintFig
# printfig(fig,filename="test.tex")


# using PyPlot
close("all")
   m = 1.98#1.98 # kg
    mu  = 0.85
    g = 9.81 # m/s^2
    I_z = 0.03 # kg * m^2
    B = 6.0#1.0  8.0
    C = 1.6    #1.35
 FMax = mu*m*g / 2.0 

alpha_f = vec(-0.7:0.005:0.7)

F_yf = FMax * sin(C*atan(B*alpha_f))

F_yf2 = 2*FMax *alpha_f
# F_yf = -10*FMax * alpha_f
fig = figure(1)
ax1 = subplot(111)
ax1[:plot](alpha_f,F_yf)
ax1[:plot](alpha_f,F_yf2)
ax1[:grid]()
ax1[:set_xlabel](L"\alpha \> [rad] " , fontsize = 24)
ax1[:set_ylabel](L"F_y \> [N] ", fontsize = 24)
ax1[:set_xlim](-0.7,0.7)