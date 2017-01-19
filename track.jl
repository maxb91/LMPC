using JLD
using PyPlot

function create_track(w = 0.3)
    close("all")
    x = [0.0]           # starting point
    y = [0.0]
    x_l = [0.0]           # starting point
    y_l = [w]
    x_r = [0.0]           # starting point
    y_r = [-w]
    ds = 0.1

    theta = [0.0]

    # SOPHISTICATED TRACK
    # add_curve(theta,30,0.0)
    # add_curve(theta,60,-2*pi/3)
    # add_curve(theta,90,pi)
    # add_curve(theta,80,-5*pi/6)
    # add_curve(theta,10,0.0)
    # add_curve(theta,50,-pi/2)
    # add_curve(theta,50,0.0)
    # add_curve(theta,40,-pi/4)
    # add_curve(theta,30,pi/4)
    # add_curve(theta,20,0.0)
    # add_curve(theta,50,-pi/2)
    # add_curve(theta,25,0.0)
    # add_curve(theta,50,-pi/2)
    # add_curve(theta,28,0.0)


     # # racetarck
    add_curve(theta,20,0.0)
    # add_curve(theta,40,-pi*3/4)
    add_curve(theta,34,-pi/2)
    add_curve(theta,10,-pi/5)
    add_curve(theta,68,0.0)
    add_curve(theta,10,-pi/4)
    add_curve(theta,62,-pi)
    add_curve(theta,10,pi/5)
    add_curve(theta,11,0)
    add_curve(theta,32,pi*3/4)
    add_curve(theta,30,0.0)
    add_curve(theta,94,-pi)
    add_curve(theta,12,-pi*2/5*0.9366)
    add_curve(theta,12,pi*2/5*0.9366)
    add_curve(theta,17,0.0)


    # # SIMPLE oval
    # add_curve(theta,7,0)
    # add_curve(theta,200,-pi)
    # add_curve(theta,14,0)
    # add_curve(theta,200,-pi)
    # add_curve(theta,5,0)

    # GOGGLE TRACK
    # add_curve(theta,30,0)
    # add_curve(theta,40,-pi/2)
    # add_curve(theta,40,-pi/2)
    # add_curve(theta,20,-pi/6)
    # add_curve(theta,30,pi/3)
    # add_curve(theta,20,-pi/6)
    # add_curve(theta,40,-pi/2)
    # add_curve(theta,40,-pi/2)
    # add_curve(theta,35,0)

    # SIMPLE GOGGLE TRACK
    # add_curve(theta,30,0)
    # add_curve(theta,40,-pi/2)
    # add_curve(theta,10,0)
    # add_curve(theta,40,-pi/2)
    # add_curve(theta,20,pi/10)
    # add_curve(theta,30,-pi/5)
    # add_curve(theta,20,pi/10)
    # add_curve(theta,40,-pi/2)
    # add_curve(theta,10,0)
    # add_curve(theta,40,-pi/2)
    # add_curve(theta,35,0)

    #  # SHORT SIMPLE track
    # add_curve(theta,10,0)
    # add_curve(theta,80,-pi)
    # add_curve(theta,20,0)
    # add_curve(theta,80,-pi)
    # add_curve(theta,9,0)

    for i=1:length(theta)
            push!(x, x[end] + cos(theta[i])*ds)
            push!(y, y[end] + sin(theta[i])*ds)
            push!(x_l, x[end-1] + cos(theta[i]+pi/2)*w)
            push!(y_l, y[end-1] + sin(theta[i]+pi/2)*w)
            push!(x_r, x[end-1] + cos(theta[i]-pi/2)*w)
            push!(y_r, y[end-1] + sin(theta[i]-pi/2)*w)
    end
    # track = cat(2, x, y, x_l, y_l, x_r, y_r)
    track =[x'; y']
    plot(track[1,:],track[2,:])
    grid()
    return track
    # plot(x,y,x_l,y_l,x_r,y_r)
end

function add_curve(theta::Array{Float64}, length::Int64, angle)
    d_theta = 0
    curve = 2*sum(1:length/2)+length/2
    for i=0:length-1
        if i < length/2+1
            d_theta = d_theta + angle / curve
        else
            d_theta = d_theta - angle / curve
        end
        push!(theta, theta[end] + d_theta)
    end
end

function save_track(filenameS)
    jldopen(filenameS, "w") do file
        #addrequire(file, classes) #ensures that custom data types are working when loaded
        JLD.write(file, "x", track[1,:])
        JLD.write(file, "y", track[2,:])
    end
end