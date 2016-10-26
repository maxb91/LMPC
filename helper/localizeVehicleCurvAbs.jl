function localizeVehicleCurvAbs(states_x,x_track,y_track,TrackCoeff, itercount)
    # Outputs: zCurr_s, coeffCurv 
    # zCurr_s = [s, ey, epsi, states_x[4]]
    # itercount is solely used for debugging purposes plots

    OrderXY = TrackCoeff.nPolyXY
    OrderThetaCurv = TrackCoeff.nPolyCurvature

    # grab current states of the vehicle
    x       = states_x[1]
    y       = states_x[2]
    psi     = states_x[3]

    nodes_center      = [x_track; y_track]
     # if (x == 0) && (y == 0) && (psi == 0)
     #     x = nodes_center[1,31]
     #     y = nodes_center[2,31]
     # end


    # Select order of interpolation (Now hard coded)
    nPoints = 60
    N_nodes_poly_back = 30 
    N_nodes_poly_front = nPoints-N_nodes_poly_back



    #--- HERE I am assuming that the distance between two points in the map is
    #--- 1M <----------- IMPORTANT ASSUMPTION --------------------------------
    Matrix = zeros(nPoints+1,OrderXY+1)
    for i = 0:nPoints
        for k = 0:OrderXY
            Matrix[i+1,OrderXY+1-k] = i^k
        end
    end

    Matrix4th = zeros(nPoints+1,OrderThetaCurv+1) #generate a matrix of 4th order to approximate Theta and the curvature
    for i = 0:nPoints
        for k = 0:OrderThetaCurv
            Matrix4th[i+1,OrderThetaCurv+1-k] = i^k
        end
    end


    #--- NEED to check if in nodes_lane there are no point repeated. I had
    #this problem when I used the map that Jason gave me
        
    # distances from current point to all relevant points
    dist_vec        = sqrt((x-x_track).^2 + (y - y_track).^2)
        
    # find minimum distance 
    dist_lane, idx_min = findmin(dist_vec)

    # i assumed that idx_min is the index of the closest node  

    N_nodes_center  = size(nodes_center,2) # how many nodes define the track 

    if idx_min <= N_nodes_poly_back #if the nearest point is less meters away from the start of the track then we need to change the structure and interpolate just from 1 to the future
        ind_start   = 1     #this has to be tested at the moment we just start at s= 30 to prevetn this part
        ind_end     = nPoints+1
        error("too close to start") #because this does not work yet we set an error msg
    else 
        ind_start   = idx_min-N_nodes_poly_back #max(1,idx_min-N_nodes_poly_back)
        ind_end     = idx_min+N_nodes_poly_front #min(N_nodes_center,idx_min+N_nodes_poly_front)
    end
    nodes_near = nodes_center[:,ind_start:ind_end]

    # Select node for X Y
    nodes_near_X = vec(nodes_near[1,:])
    nodes_near_Y = vec(nodes_near[2,:])
    
    # xyVectorModule=((x-nodes_near_X(ind_start))^2+(y-nodes_near_Y(ind_start))^2)^(1/2);
    # compute least squares coefficients that approximate the road
    coeffY = Matrix\nodes_near_Y
    coeffX = Matrix\nodes_near_X  

    # now compute (s,y) needed for the MPC
    Counter = 1
    Discretization = 0.01
    # Initializate the values to find exact s look 1 meter in front and 1 meter in the back of nearest point no track
    S_Value = zeros(0:Discretization:2) #create an vector for the elements 1 meter before and after the nearest point of the track 
    DistanceNew = zeros(0:Discretization:2)

    #T those were just to test
    # XCurve_t = zeros(0:Discretization:2)
    # YCurve_t = zeros(0:Discretization:2)

    # Evaluate all points to find the current s cloesest to vehicle
    j_vec = zeros(OrderXY+1)                     

    #!! this is only correct if the car is at a as greater 30  we need to adjust a if statement
    if idx_min <= N_nodes_poly_back #if the nearest point is less meters away from the start of the track then we need to change the structure and interpolate just from 1 to the future
    #use point from end of track or make some points up
    error("not implemeted")
    else   
        for j=(N_nodes_poly_back-1):Discretization:(N_nodes_poly_back+1) #?? j between 29 and 31 should be between idx_min+-1?
            # j does not stand for the id of the node but for the length of s in meters so j = 29 means node number 30
            for i = 1:OrderXY+1
                j_vec[i] =j^(OrderXY+1-i)
            end
            XCurve = dot(coeffX, j_vec)
            YCurve = dot(coeffY, j_vec)


            #T just test
            # XCurve_t[Counter] = dot(coeffX, j_vec)
            # YCurve_t[Counter] = dot(coeffY, j_vec)


            S_Value[Counter] = j
            DistanceNew[Counter] = sqrt((x-XCurve).^2+(y-YCurve).^2) #distance of vehicle to every interpolated node
            Counter = Counter + 1
        end
    end
    # T
    # plot(XCurve_t,YCurve_t , linewidth=2.0)
    # scatter(XCurve_t[1],YCurve_t[1])
    # scatter(XCurve_t[end],YCurve_t[end])
    # scatter(x,y)
    # plot(nodes_near_X, nodes_near_Y ) 
    #scatter(x_track[idx_min],y_track[idx_min])
  

    # from the evaluated points get the best as [s, y] ---> use the as Feedback
    eyabs, idx_min_Dist = findmin(DistanceNew)

    # Extract the current s
    s=S_Value[idx_min_Dist] #s is always between 29 and 31 have to ad s start for real position

    # Find the sign of ey #?? how does this work?
    s0 = s-0.01 
    s0_vec = zeros(OrderXY+1,1)
    s_vec = zeros(OrderXY+1,1)
    sdot_vec = zeros(OrderXY+1,1)::Array{Float64,2}

    for i = 1:OrderXY+1
            s_vec[i] =s^(OrderXY+1-i)
            s0_vec[i] =s0^(OrderXY+1-i)
            sdot_vec[i] = (OrderXY+1-i)*s^(OrderXY-i)
    end

    XCurve0 = coeffX'*s0_vec
    YCurve0 = coeffY'*s0_vec

    XCurve = coeffX'*s_vec
    YCurve = coeffY'*s_vec
    dX = coeffX'*sdot_vec #[6*s^5 5*s^4 4*s^3 3*s^2 2*s^1 1 0]' comment can be deleted if sdot_vec is verified for all polynomials
    dY = coeffY'*sdot_vec      

    xyVectorAngle = atan2(y-YCurve0,x-XCurve0)
    xyPathAngle = atan2(dY,dX)
    #?? calculate epsi here?

   @show ey = eyabs*sign(sin(xyVectorAngle-xyPathAngle))

    # Calcuate the error due to the conversion in the curvilinear abscissa #?? what does it do and how
    yBack = YCurve + ey*cos(xyPathAngle)
    xBack = XCurve - ey*sin(xyPathAngle)
    @show Error = sqrt((y-yBack)^2 + (x-xBack)^2)

    # now compute the angle and the curvature needed for the interpolation

    b_theta_vector = zeros(nPoints+1,1)
    b_curvature_vector = zeros(nPoints+1,1)
    angle = 0.0


    Counter = 1
    jd_vec = zeros(OrderXY+1,1)::Array{Float64,2}
    jdd_vec = zeros(OrderXY+1,1)::Array{Float64,2}
    #we go from s = 0 because s is 0 at s_start and then we interpolate over the interval we defined above 
    for j=0.0:nPoints #j must be 0.0 to be initialized as a float to be able to do j^-1 in for loop
    
       
        #this generic approach did not work because last elements become NaN
        #for i = 1:OrderXY+1
        #   jd_vec[i] =(OrderXY+1-i)*j^(OrderXY-i)
        #    jdd_vec[i] =(OrderXY+1-i)*(OrderXY-i)*j^(OrderXY-1-i)
        #end
        dX = dot(coeffX,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0]) 
        dY = dot(coeffY,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0])
        ddX = dot(coeffX,[30*j^4, 20*j^3, 12*j^2, 6*j, 2, 0, 0])
        ddY = dot(coeffY,[30*j^4, 20*j^3, 12*j^2, 6*j, 2, 0, 0])
            
        angle = atan2(dY,dX);

        psi
        if Counter > 1 #what do we do here? do we normalize angle when we drive against "Lane direction"
            DummyVar = angle-(b_theta_vector[Counter-1,1]) #?? angle(1)
            if (DummyVar > pi)
                angle = angle - 2*pi
            elseif DummyVar < -pi
                angle = angle + 2*pi
            end
        else
            if (angle-psi) > pi #?? angle(1)
                angle = angle - 2*pi #?? angle(1)
            elseif (angle-psi) < -pi #?? angle(1)
                angle = angle + 2*pi #?? angle(1)
            end
        end
            
        b_theta_vector[Counter,1] = angle #?? angle(1)
            
        curvature = (dX*ddY-dY*ddX)/(dX^2+dY^2)^(3/2) #standard curvature formula
        b_curvature_vector[Counter,1] = curvature #?? curvature(1)
           
        Counter = Counter + 1
    end
  
    # compute coeff for theta and curvature
    coeffTheta=Matrix4th\b_theta_vector
    coeffCurv=Matrix4th\b_curvature_vector
    coeffCurv = vec(coeffCurv)
   
   #####T
   #test the approximation of the curvature
   # if itercount%20 == 0 
        

   #      #EXTRACT the caluclated curvatue form drivatives formula
   #      curv_from_formula = b_curvature_vector[1:nPoints+1,1]
   #      s_curv=collect(0.0:nPoints)
        
   #      #calculate the curvature back from the coefficients
   #      Count = 1
   #      polyt = zeros(201)
   #      for s_t = 10:0.2:50
   #      polyt[Count] = dot(coeffCurv, [s_t^4 s_t^3 s_t^2 s_t 1])
   #      Count+=1
   #      end
   #      s_t = collect(10:0.2:50)
   #      # # s_t = zeros(25:0.1:35)
   #      # # for i = 1:size(s_t)[1]
   #      # # s_t[i] =25+0.1*(i-1)
   #      # # end

   #      #plot
   #      close(3)
   #      figure(3)
   #      plot(s_curv, curv_from_formula)
   #      plot(s_t,polyt)
   #      legend(["curv_from_formula","curv from back calculation"])
   #  end
    #####endT


    XConvertedBackS=coeffX'*[s^6 s^5 s^4 s^3 s^2 s 1]' #?? what means convertedbackS?
    YConvertedBackS=coeffY'*[s^6 s^5 s^4 s^3 s^2 s 1]'   
        
    angle=coeffTheta'*[s^4 s^3 s^2 s 1]'; #Make the sketch to understand why #??

        

    XConvertedBack = XConvertedBackS - sin(angle)*ey #??
    YConvertedBack = YConvertedBackS + cos(angle)*ey
    ErrorConvertedBack = sqrt((x-XConvertedBack[1])^2+((y-YConvertedBack[1])^2))
        

    # Finally compute epsi
    j = s
    dX = dot(coeffX,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0])
    dY = dot(coeffY,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0])
    angle=atan2(dY,dX) # gives a value between -pi < x <= pi

    if angle < 0
        angle = 2*pi + angle
  
    end

    epsi = psi - angle

    if abs(epsi)>(pi/2)#?? what do we normalize?
        if epsi<(pi/2)
            epsi = psi - (angle - 2*pi)  
        else
            epsi = (psi - 2*pi) - angle          
        end
    end

    #T
    #this was yjust to test correcntess of teh dfrisr derivative
    # posTrackx = x_track[idx_min]
    # posTracky = y_track[idx_min]
    # drawtangentx =posTrackx + cos(angle)*3
    # drawtangenty =posTracky+ sin(angle)*3
    # xs= [posTrackx; drawtangentx]
    # ys= [posTracky; drawtangenty]
    
    #changeMetersforBARC
   # @show s_start = ind_start/10-0.1 #--- HERE assuming all point are 1m distant, # if we are at the 4th point we are 3m away from the start
   # s= s/10
    s_start = ind_start - 1

    #return s_start, s, ey, coeffX,coeffY, coeffTheta, coeffCurv, epsi
    zCurr_s = zeros(4)
    zCurr_s = [s+s_start ey epsi states_x[4]]
    return zCurr_s, coeffCurv, s_start
end
