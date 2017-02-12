function loadTestMap()
 #x = readdlm("C:/Users/Felix/Desktop/Project/x.txt",'  ')
 #y = readdlm("C:/Users/Felix/Desktop/Project/y.txt",'  ')


# x = Array(0:0.1:69.9)'
# y = zeros(700)'

# x = readcsv("C:/Users/Felix/AppData/Local/Julia-0.5.0/LMPC/helper/x_dm.dat")
# y = readcsv("C:/Users/Felix/AppData/Local/Julia-0.5.0/LMPC/helper/y_dm.dat")

# file = "tracks/trackdata.jld"
file = "tracks/oval.jld"
# file = "tracks/small_oval.jld"
# file = "tracks/adv_track.jld"
# file = "tracks/optm_adv.jld"
# file = "tracks/sim_curv.jld" #similar curvature in wo curves
Data = load(file)
x = Data["x"]'
y = Data["y"]'
return x, y
end