function localizeVehicleCurvAbs(states_x::Array{Float64},x_track::Array{Float64},y_track::Array{Float64},trackCoeff::classes.TrackCoeff, itercount::Int64, N::Int64, dt::Float64 , Pcurvature)
    #Pcurvature just for plotting bugfixing
    # Outputs: zCurr_s, coeffCurv 
    # zCurr_s = [s, ey, epsi, states_x[4]]
    # itercount is solely used for debugging purposes plots

    OrderXY = trackCoeff.nPolyXY
    OrderThetaCurv = trackCoeff.nPolyCurvature
    ds = trackCoeff.ds

    # grab current states of the vehicle
    x       = states_x[1]
    y       = states_x[2]
    psi     = states_x[5]
    v_x     = states_x[3]
    v_y     = states_x[4]
    psi_dot     = states_x[6]
    v_abs = sqrt(states_x[3].^2 + states_x[4].^2)

    nodes_center      = [x_track; y_track]
     # if (x == 0) && (y == 0) && (psi == 0)
     #     x = nodes_center[1,31]
     #     y = nodes_center[2,31]
     # end


    # Select order of interpolation (Now hard coded)

    # N_nodes_poly_back = 10
    # N_nodes_poly_front = 30
    # nPoints = N_nodes_poly_back + N_nodes_poly_front
    
    N_nodes_poly_back = convert(Int64,ceil(0.5*v_abs*N*dt/trackCoeff.ds))
    N_nodes_poly_front = convert(Int64,ceil(2.0*v_abs*N*dt/trackCoeff.ds))
    nPoints = N_nodes_poly_back + N_nodes_poly_front


    #--- HERE I am assuming that the distance between two points in the map is
    #--- 1M <----------- IMPORTANT ASSUMPTION --------------------------------


    #--- NEED to check if in nodes_lane there are no point repeated. I had
    #this problem when I used the map that Jason gave me
        
    # distances from current point to all relevant points
    dist_vec        = sqrt((x-x_track).^2 + (y - y_track).^2)
        
    # find minimum distance 
    dist_lane, idx_min = findmin(dist_vec)

    # i assumed that idx_min is the index of the closest node  

    N_nodes_center  = size(nodes_center,2) # how many nodes define the track 
    ind_start   = idx_min-N_nodes_poly_back #max(1,idx_min-N_nodes_poly_back)
    ind_end     = idx_min+N_nodes_poly_front #min(N_nodes_center,idx_min+N_nodes_poly_front)
    
    if idx_min <= N_nodes_poly_back #if the nearest point is less meters away from the start of the track then we need to change the structure and interpolate just from 1 to the future
          ###not closed used this:
        # ind_start   = 1
        # ind_end     = nPoints+1
        # nodes_near = nodes_center[:,ind_start:ind_end]
        #if closed use this:
        nodes_near= hcat(nodes_center[:,N_nodes_center+ind_start:N_nodes_center],nodes_center[:,1:ind_end]) 
    elseif idx_min+N_nodes_poly_front >= N_nodes_center
       ###not closed used this:
        # ind_start   = N_nodes_center-nPoints 
        # ind_end     = N_nodes_center
        # nodes_near = nodes_center[:,ind_start:ind_end]
        #if closed use this:stack the end and beginning of the lap together
       
        # ind_start   = idx_min-N_nodes_poly_back 
        # ind_end     = idx_min+N_nodes_poly_front
        # nodes_near = append!(nodes_center[:,ind_start:N_nodes_center],nodes_center[:,0:idx_min+N_nodes_poly_front+1-N_nodes_center])
         nodes_near= hcat(nodes_center[:,ind_start:N_nodes_center],nodes_center[:,1:ind_end-N_nodes_center]) 
        # error("stacking nodes from beginning of track not yet implemented")
    else 
        
        nodes_near = nodes_center[:,ind_start:ind_end]
    end
    

    s_interp_start = (ind_start-1)*ds # in the first index s is equal 0. (always one less)
    #s_interp_end = (ind_end-1)*ds
    s_nearest = (idx_min-1)*ds
    # Select node for X Y
    nodes_near_X = vec(nodes_near[1,:])
    nodes_near_Y = vec(nodes_near[2,:])

    index = 0 
#this just works because s = 1 m between points
        Matrix = zeros(nPoints+1,OrderXY+1)
    # for i = s_interp_start:ds:s_interp_end
    for i = s_interp_start:ds:s_interp_start+nPoints*ds
        for k = 0:OrderXY
            index = convert(Int,(i-s_interp_start)/ds+1)
            Matrix[index,OrderXY+1-k] = i^k
        end
    end

    Matrix4th = zeros(nPoints+1,OrderThetaCurv+1) #generate a matrix of 4th order to approximate Theta and the curvature
    for i = s_interp_start:ds:s_interp_start+nPoints*ds
        for k = 0:OrderThetaCurv
            index = convert(Int,(i-s_interp_start)/ds+1)
            Matrix4th[index,OrderThetaCurv+1-k] = i^k
        end
    end

    
    # xyVectorModule=((x-nodes_near_X(ind_start))^2+(y-nodes_near_Y(ind_start))^2)^(1/2);
    # compute least squares coefficients that approximate the road
    coeffY = Matrix\nodes_near_Y
    coeffX = Matrix\nodes_near_X  

#######test
    # @show g=(idx_min-1)*ds
    # testdX = dot(coeffX,[6*g^5, 5*g^4, 4*g^3, 3*g^2, 2*g, 1, 0]) 
    #  @show testX = dot(coeffX,[g^6, g^5, g^4, g^3, g^2, g, 1]) 
    #######
    # now compute (s,y) needed for the MPC
    Counter = 1
    Discretization = 0.01*ds
    j_vec = zeros(OrderXY+1) 
    # Initializate the values to find exact s look 1 meter in front and 1 meter in the back of nearest point no track
   
    #T those were just to test
    # XCurve_t = zeros(0:Discretization:2)
    # YCurve_t = zeros(0:Discretization:2)

    # Evaluate all points to find the current s cloesest to vehicle

   
    if s_nearest >= ds
        S_Value = zeros(0:Discretization:2*ds) #create an vector for the elements 1 meter before and after the nearest point of the track 
        DistanceNew = zeros(0:Discretization:2*ds)

        for j=(s_nearest-ds):Discretization:(s_nearest+ds) 
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
    else # just active while the car is nearest to the starting point as we cannot evaluate functon for negative s
        S_Value = zeros(0:Discretization:1*ds) #create an vector for the elements 1 after the nearest point of the track 
        DistanceNew = zeros(0:Discretization:1*ds)

        #T those were just to test
        # XCurve_t = zeros(0:Discretization:2)
        # YCurve_t = zeros(0:Discretization:2)
        for j=s_nearest:Discretization:(s_nearest+ds) 
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
    # scatter(XCurve_t[1:100:end],YCurve_t[1:100:end] , color= "yellow")
    # scatter(XCurve_t[1],YCurve_t[1],color= "blue")
    # scatter(XCurve_t[end],YCurve_t[end],color= "blue")
    # scatter(x,y,color ="black")
    # plot(nodes_near_X, nodes_near_Y ) 
    # scatter(x_track[idx_min],y_track[idx_min], color = "green")
  

    # from the evaluated points get the best as [s, y] ---> use the as Feedback
    #@show 
    eyabs, idx_min_Dist = findmin(DistanceNew)

    # Extract the current s
    s=S_Value[idx_min_Dist] #s is always between 29 and 31 have to ad s start for real position

    # Find the sign of ey
    if s >= 0.01*ds
        s0 = s-0.01*ds
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
    else
        s0 = s+0.01*ds
        s0_vec = zeros(OrderXY+1,1)
        s_vec = zeros(OrderXY+1,1)
        sdot_vec = zeros(OrderXY+1,1)::Array{Float64,2}

        for i = 1:OrderXY+1
                s_vec[i] =s^(OrderXY+1-i)
                s0_vec[i] =s0^(OrderXY+1-i)
               #sdot_vec = (OrderXY+1-i)*s^(OrderXY-i) #this calucaltion does not work: problems with NaN last i  gives 0 *Inf =NaN
        end
        sdot_vec = [6*s^5 5*s^4 4*s^3 3*s^2 2*s 1 0]'

        XCurve0 = coeffX'*s0_vec
        YCurve0 = coeffY'*s0_vec

        XCurve = coeffX'*s_vec
        YCurve = coeffY'*s_vec
           
        dX = dot(coeffX,sdot_vec) #[6*s^5 5*s^4 4*s^3 3*s^2 2*s^1 1 0]' comment can be deleted if sdot_vec is verified for all polynomials
        dY = dot(coeffY,sdot_vec) 
        xyVectorAngle = atan2(y-YCurve0,x-XCurve0)
        xyPathAngle = atan2(dY,dX)# gives a value between -pi < x <= pi  
    end
    ey = eyabs*sign(sin(xyVectorAngle-xyPathAngle))

  #compute epsi

    # s_vec = zeros(OrderXY+1,1)
    # sdot_vec = zeros(OrderXY+1,1)::Array{Float64,2}

    # for i = 1:OrderXY+1
    #         s_vec[i] =s^(OrderXY+1-i)
    #         sdot_vec[i] = (OrderXY+1-i)*s^(OrderXY-i)
    # end
    # sdot_vec = [6*s^5 5*s^4 4*s^3 3*s^2 2*s 1 0]'


    # dX = dot(coeffX,sdot_vec) #[6*s^5 5*s^4 4*s^3 3*s^2 2*s^1 1 0]' comment can be deleted if sdot_vec is verified for all polynomials
    # # dY = dot(coeffY,sdot_vec) 
    # xyPathAngle = atan2(dY,dX)
    epsi = mod((psi+pi),(2*pi))-pi-xyPathAngle
    epsi = mod((epsi+pi),(2*pi))-pi

     # ey = eyabs*sign(epsi)


  # # Finally compute epsi
  #   j = s
  #   dX = dot(coeffX,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0])
  #   dY = dot(coeffY,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0])
  #   angle=atan2(dY,dX) # gives a value between -pi < x <= pi

  #   # epsi = psi - angle


  #   epsi = mod((psi+pi),(2*pi))-pi-angle
  #   epsi = mod((epsi+pi),(2*pi))-pi


    # XCurve= dot(coeffX,s_vec)
    # YCurve= dot(coeffY,s_vec)
    #T Calcuate the error due to the conversion in the curvilinear abscissa
    yBack = YCurve + ey*cos(xyPathAngle)
    xBack = XCurve - ey*sin(xyPathAngle)
    Error = sqrt((y-yBack)^2 + (x-xBack)^2)
    if Error[1] >= 0.001
        warn("problem with approximation of x and y pos. Error: $Error, i: $itercount")
    end
    #endT

    # now compute the angle and the curvature needed for the interpolation

    b_theta_vector = zeros(nPoints+1,1)
    b_curvature_vector = zeros(nPoints+1)
    angle = 0.0


    Counter = 1
    jd_vec = zeros(OrderXY+1,1)::Array{Float64,2}
    jdd_vec = zeros(OrderXY+1,1)::Array{Float64,2}
    #we go from s = 0 because s is 0 at s_start and then we interpolate over the interval we defined above 
    #for j=s_interp_start:ds:s_interp_end #j must be 0.0 to be initialized as a float to be able to do j^-1 in for loop
    for j = s_interp_start:ds:s_interp_start+nPoints*ds
       
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

       
        if Counter > 1 #what do we do here? do we normalize angle when we drive against "Lane direction"
            DummyVar = angle-(b_theta_vector[Counter-1,1]) 
            if (DummyVar > pi)
                angle = angle - 2*pi
            elseif DummyVar < -pi
                angle = angle + 2*pi
            end
        else
            if (angle-psi) > pi
                angle = angle - 2*pi
            elseif (angle-psi) < -pi 
                angle = angle + 2*pi 
            end
        end
            
        b_theta_vector[Counter,1] = angle 
            
        curvature = (dX*ddY-dY*ddX)/(dX^2+dY^2)^(3/2) #standard curvature formula
        b_curvature_vector[Counter] = curvature #there might still be problesm with the accuracy of the curvature calculation
           
        Counter = Counter + 1
    end
  
    # compute coeff for theta and curvature
    coeffTheta = Matrix4th\b_theta_vector
    coeffCurv  = Matrix4th\b_curvature_vector
   #####T
   #test the approximation of the curvature
   # if itercount%20 == 0 
        

   #      #EXTRACT the caluclated curvatue form drivatives formula
   #      curv_from_formula = b_curvature_vector[1:nPoints+1,1]
   #      s_curv=collect(0.0:nPoints)
        
   #      #calculate the curvature back from the coefficients
        # Count = 1
        # s_p = collect(s-N_nodes_poly_back*ds:0.01:s+N_nodes_poly_front*ds)
        # polyt = zeros(size(s_p))
        # for s_t = s-N_nodes_poly_back*ds:0.01:s+N_nodes_poly_front*ds
        #     polyt[Count] = dot(coeffCurv, [s_t^4 s_t^3 s_t^2 s_t 1])
        #     Count+=1
        # end
        # plot(s_p,polyt, color = "blue")
        # readline()


        j = s
        dX = dot(coeffX,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0]) 
        dY = dot(coeffY,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0])
        ddX = dot(coeffX,[30*j^4, 20*j^3, 12*j^2, 6*j, 2, 0, 0])
        ddY = dot(coeffY,[30*j^4, 20*j^3, 12*j^2, 6*j, 2, 0, 0])
        Pcurvature[itercount,1] = s
        Pcurvature[itercount,2] = (dX*ddY-dY*ddX)/(dX^2+dY^2)^(3/2)
        
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

 
        
    # angle=coeffTheta'*[s^4 s^3 s^2 s 1]'; #Make the sketch to understand why 

    
        

  

    #T
    #this was yjust to test correcntess of teh dfrisr derivative
    # posTrackx = x_track[idx_min]
    # posTracky = y_track[idx_min]
    # drawtangentx =posTrackx + cos(angle)*3
    # drawtangenty =posTracky+ sin(angle)*3
    # xs= [posTrackx; drawtangentx]
    # ys= [posTracky; drawtangenty]
    
    

    #return s_start, s, ey, coeffX,coeffY, coeffTheta, coeffCurv, epsi
    zCurr_s = zeros(6)
    zCurr_s = [v_x v_y psi_dot epsi ey s] 
    return zCurr_s, coeffCurv
end
