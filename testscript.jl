fig = figure()
ax1 = subplot(1,1,1)
drawPatch=drawCar(fig,[0.0,0.0,0.0], "green")

for p in drawPatch
                ax1[:add_patch](p)
end

