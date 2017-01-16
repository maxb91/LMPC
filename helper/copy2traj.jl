using JLD

include("classes.jl")
file = "data/2017-01-14-16-19-Data.jld"
Data = load(file)
obstacle = Data["obstacle"]
oldTraj     = Data["oldTraj"]



k=1

v=obstacle.v[:,k]
sY = obstacle.sy_obstacle[:,k]
olTra=oldTraj.oldTraj[:,:,k]
curvature=oldTraj.curvature[:,k]
cost=oldTraj.cost2Target[:,k]
distanc=oldTraj.distance2obst[:,k]

file = "data/2017-01-14-16-36-Data.jld"
Data = load(file)
obstacle = Data["obstacle"]
oldTraj = Data["oldTraj"]

k=10
obstacle.v[:,k]=v
obstacle.sy_obstacle[:,k] =sY
oldTraj.oldTraj[:,2:4,k]=olTra[:,2:4]
oldTraj.oldTraj[:,1,k]=olTra[:,1]
oldTraj.curvature[:,k]=curvature
oldTraj.cost2Target[:,k]=cost
oldTraj.distance2obst[:,k]=distanc


filenameS = string("data/"string(Dates.today()),"-",Dates.format(now(), "HH-MM"),"-SafeSet.jld")
if isfile(filenameS)
        filenameS = string("data/"string(Dates.today()),"-",Dates.format(now(), "HH-MM"),"-SafeSet-2.jld")
        warn("SafeSet file already exists. Added extension \"-2\" ")
end
println("Save SafeSet to $filenameS .......")
jldopen(filenameS, "w") do file
    #addrequire(file, classes) #ensures that custom data types are working when loaded
    JLD.write(file, "oldTraj", oldTraj)
    JLD.write(file,"obstacle", obstacle)
end
println("finished")

