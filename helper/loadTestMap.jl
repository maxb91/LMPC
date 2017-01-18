function loadTestMap()
 #x = readdlm("C:/Users/Felix/Desktop/Project/x.txt",'  ')
 #y = readdlm("C:/Users/Felix/Desktop/Project/y.txt",'  ')


# x = Array(0:0.1:69.9)'
# y = zeros(700)'

# x = readcsv("C:/Users/Felix/AppData/Local/Julia-0.5.0/LMPC/helper/x_dm.dat")
# y = readcsv("C:/Users/Felix/AppData/Local/Julia-0.5.0/LMPC/helper/y_dm.dat")

# file = "trackdata.jld"
# file = "oval.jld"
file = "small_oval.jld"
Data = load(file)
x = Data["x"]'
y = Data["y"]'
return x, y
end