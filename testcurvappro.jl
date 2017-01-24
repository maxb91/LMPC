module testAppro

using PyPlot

dt = 0.1
t  = collect(0:dt:40)
zCurr_s = zeros(length(t)+1,4)  
zCurr_x = zeros(length(t)+1,4)   


include("helper/loadTestMap.jl")
x_track, y_track = loadTestMap() #load a track to test controller


type TrackCoeff         # coefficients of track
    coeffAngle::Array{Float64,1}
    coeffCurvature::Array{Float64,1}
    nPolyCurvature::Int64      # order of the interpolation polynom
    nPolyXY::Int64              # order of the interpolation polynom of the x y coordinates
    width::Float64               # lane width -> is used in cost function as soft constraints (to stay on track)
    TrackCoeff(coeffAngle=Float64[], coeffCurvature=Float64[], nPolyCurvature=4, nPolyXY = 6, width=1.0) = new(coeffAngle,coeffCurvature,nPolyCurvature,nPolyXY)
end
trackCoeff                  = TrackCoeff()   
trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]         # polynomial coefficients for curvature approximation (zeros for straight line)
trackCoeff.nPolyCurvature = 4 # has to be 4 cannot be changed freely at the moment orders are still hardcoded in some parts of localizeVehicleCurvAbslizeVehicleCurvAbs
trackCoeff.nPolyXY = 6 


include("helper/localizeVehicleCurvAbs.jl")

zCurr_x[1,1] = 140 #140 #20
zCurr_x[1,2] = 126 #126 #25
zCurr_x[1,3] = 0.1


i=1
zCurr_s[i,:], trackCoeff.coeffCurvature, s_start = localizeVehicleCurvAbs(zCurr_x[i,:],x_track,y_track,trackCoeff)
end


using PyPlot
close("all")
   m = 1.98#1.98 # kg
    mu  = 0.85
    g = 9.81 # m/s^2
    I_z = 0.03 # kg * m^2
    B = 8.0#1.0  8.0
    C = 1.35    #1.35
 FMax = mu*m*g / 2.0 

alpha_f = vec(-0.7:0.005:0.7)

F_yf = -FMax * sin(C*atan(B*alpha_f))
plot(alpha_f,F_yf)
grid()
