
using JuMP
using Ipopt
include("ModelFcn.jl")

# define model parameters
L_a     = 0.125         # distance from CoG to front axel
L_b     = 0.125         # distance from CoG to rear axel
dt      = 0.1           # time step of system
coeffCurvature   = [0,0,0,0,0]

# preview horizon
N       = 10

# define targets [generic values]
v_ref   = 0.2

c0      = [0.5431, 1.2767, 2.1516, -2.4169]

p = modelParams(L_a,L_b,dt,coeffCurvature,N,v_ref,c0)

# define decision defVars 
# states: position (x,y), yaw angle, and velocity
# inputs: acceleration, steering angle 
println("Creating kinematic bicycle model ....")

InitModel()


# status update
println("initial solve ...")

status = solve(mdl)

println("finished initial solve!")
println("Solve Status: ")
println(status)
println("Changing objective...\n")

changeObjective()

println("Solve again...")
status = solve(mdl)
print("Solved, Status: ")
println(status)