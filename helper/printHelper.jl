function printPrediction(mpcSol::MpcSol)
    figure(101)
    for i=1:5
        subplot(6,1,i)
        plot(mpcSol.z[:,6],mpcSol.z[:,i],"-o")
        grid("on")
        legend(["$i"])
    end
    subplot(6,1,6)
    plot(mpcSol.z[1:end-1,6],mpcSol.u[:,1],"-o")
    grid("on")
    legend(["a"])
    readline()
end