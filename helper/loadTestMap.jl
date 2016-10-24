function loadTestMap()
 #x = readdlm("C:/Users/Felix/Desktop/Project/x.txt",'  ')
 #y = readdlm("C:/Users/Felix/Desktop/Project/y.txt",'  ')
x = readcsv("C:/Users/Felix/AppData/Local/Julia-0.5.0/LMPC/helper/x_m.dat")
y = readcsv("C:/Users/Felix/AppData/Local/Julia-0.5.0/LMPC/helper/y_m.dat")
return x, y
end