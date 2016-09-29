function [s_start, s, ey, coeffX,coeffY, coeffTeta, coeffCurv, epsi] = localizeVehicleCurvAbs(states,x_track,y_track,mpcParams)
# this function will properly localize the vehicle in all lanes 

# Outputs: 
# currLane - the current lane that the vehicle is closest to
# coeffs - 4xn array of the coefficients for each lane 
# sPositions - nx1 array of s positions of the vehicle in each lane

OrderXY = mpcParams.OrderXY;
OrderTetaCurv = mpcParams.OrderTetaCurv;

# grab current states of the vehicle
x       = states(1);
y       = states(2);
psi     = states(3);

nodes_center      = [x_track; y_track];

if (x == 0) && (y == 0) && (psi == 0)
    x = nodes_center[1,50];
    y = nodes_center[2,50];
end
# Select order of interpolation (Now hard coded)
nPoints = 60; #50
N_nodes_poly_back = +30; #20
N_nodes_poly_front = nPoints-N_nodes_poly_back; #--- Need 36m of curvilinear abscissa



#--- HERE I am assuming that the distance between two points in the map is
#--- 1M <----------- IMPORTANT ASSUMPTION --------------------------------
Matrix = zeros(nPoints+1,OrderXY+1);
for i = 0:nPoints
    for k = 0:OrderXY
        Matrix[i+1,OrderXY+1-k] = i^k;
    end
end

Matrix3rd = zeros(nPoints+1,OrderTetaCurv+1);
for i = 0:nPoints
    for k = 0:OrderTetaCurv
        Matrix3rd[i+1,OrderTetaCurv+1-k] = i^k;
    end
end


#--- NEED to check if in nodes_lane there are no point repeated. I had
#this problem when I used the map that Jason gave me
    
# distances from current point to all relevant points
dist_vec        = sqrt((x-x_track).^2 + (y - y-track).^2);
    
# find minimum distance 
dist_lane, idx_min = findmin(dist_vec);

# i assumed that idx_min is the index of the closest node  

N_nodes_center  = size(nodes_center,2); # how many nodes define the track 


if idx_min==1 # why disgtinguish just for the first? < polyback ?
    ind_start   = max(1,idx_min); #?? can idx_min be sth else?
    ind_end     = min(N_nodes_center,idx_min+N_nodes_poly_back+N_nodes_poly_front);
else 
    ind_start   = max(1,idx_min-N_nodes_poly_back);
    ind_end     = min(N_nodes_center,idx_min+N_nodes_poly_front);
end
nodes_near = nodes_center[:,ind_start:ind_end];

# Select node for X Y
nodes_near_X = nodes_near[1,:]';
nodes_near_Y = nodes_near[2,:]';

# xyVectorModule=((x-nodes_near_X(ind_start))^2+(y-nodes_near_Y(ind_start))^2)^(1/2);

# compute least squares coefficients that approximate the road
coeffY= Matrix\nodes_near_Y;
coeffX= Matrix\nodes_near_X;    

# now compute (s,y) needed for the MPC
Counter= 1;
Discretization=0.01;

# Initializate the values #?? 2m?
S_Value = 0*[0:Discretization:2];
DistanceNew = zeros(0:Discretization:2);

# Evaluate all points to find the current s
for j=[(N_nodes_poly_back-1):Discretization:(N_nodes_poly_back+1)]
    XCurve = coeffX'*[j^6 j^5 j^4 j^3 j^2 j^1 1]';
    YCurve = coeffY'*[j^6 j^5 j^4 j^3 j^2 j^1 1]';

    S_Value(Counter) = j;
    DistanceNew(Counter) = sqrt(((x-XCurve)^2+(y-YCurve)^2));
    Counter = Counter + 1;
end

# from the evaluated points get the best as [s, y] ---> use the as Feedback
[eyabs, idx_min_Dist] = min(DistanceNew);

# Exstract the current s
s=S_Value(idx_min_Dist);

# Find the sign of ey
s0 = s-0.01;
XCurve0 = coeffX'*[s0^6 s0^5 s0^4 s0^3 s0^2 s0^1 1]';
YCurve0 = coeffY'*[s0^6 s0^5 s0^4 s0^3 s0^2 s0^1 1]';

XCurve = coeffX'*[s^6 s^5 s^4 s^3 s^2 s^1 1]';
YCurve = coeffY'*[s^6 s^5 s^4 s^3 s^2 s^1 1]';
dX = coeffX'*[6*s^5 5*s^4 4*s^3 3*s^2 2*s^1 1 0]';
dY = coeffY'*[6*s^5 5*s^4 4*s^3 3*s^2 2*s^1 1 0]';        

xyVectorAngle = atan2(y-YCurve0,x-XCurve0);
xyPathAngle=atan2(dY,dX);

ey = eyabs*sign(sin(xyVectorAngle-xyPathAngle));

# Calcuate the error due to the conversion in the curvilinear abscissa
yBack = YCurve + ey*cos(xyPathAngle);
xBack = XCurve - ey*sin(xyPathAngle);
Error = sqrt((y-yBack)^2 + (x-xBack)^2);

# now compute the angle and the curvature needed for the interpolation
Counter= 1;
b_teta_vector = zeros(nPoints+1,1);
b_curvature_vector = zeros(nPoints+1,1);
angle = 0;

for j=[0:nPoints]
    dX = coeffX'*[6*j^5 5*j^4 4*j^3 3*j^2 2*j^1 1 0]';
    dY = coeffY'*[6*j^5 5*j^4 4*j^3 3*j^2 2*j^1 1 0]';
    ddX = coeffX'*[30*j^4 20*j^3 12*j^2 6*j^1 2 0 0]';
    ddY = coeffY'*[30*j^4 20*j^3 12*j^2 6*j^1 2 0 0]';
        
    angle=atan2(dY,dX);
    
    if Counter > 1
        DummyVar = angle(1)-(b_teta_vector(Counter-1,1));
        if (DummyVar > pi)
            angle = angle - 2*pi;
        elseif (DummyVar) < -pi
            angle = angle + 2*pi;
        end
    else
        if (angle(1)-psi) > pi
            angle = angle(1) - 2*pi;
        elseif (angle(1)-psi) < -pi
            angle = angle(1) + 2*pi;
        end
    end
        
    b_teta_vector(Counter,1) = angle(1);
        
    curvature = (dX*ddY-dY*ddX)/(dX^2+dY^2)^(3/2);
    b_curvature_vector(Counter,1) = curvature(1);
       
    Counter = Counter + 1;
end
    
# compute coeff for teta and curvature
coeffTeta=Matrix3rd\b_teta_vector; 
coeffCurv=Matrix3rd\b_curvature_vector; 

XConvertedBackS=coeffX'*[s^6 s^5 s^4 s^3 s^2 s 1]';
YConvertedBackS=coeffY'*[s^6 s^5 s^4 s^3 s^2 s 1]';   
    
angle=coeffTeta'*[s^3 s^2 s 1]'; #Make the sketch to understand why
    
XConvertedBack = XConvertedBackS - sin(angle)*(ey);
YConvertedBack = YConvertedBackS + cos(angle)*(ey);
ErrorConvertedBack = sqrt((x-XConvertedBack(1))^2+((y-YConvertedBack(1))^2));
    

# Finally compute epsi
j = s;
dX = coeffX'*[6*j^5 5*j^4 4*j^3 3*j^2 2*j^1 1 0]';
dY = coeffY'*[6*j^5 5*j^4 4*j^3 3*j^2 2*j^1 1 0]';
angle=atan2(dY,dX);

if angle<0
    angle = 2*pi + angle;
end

epsi = psi - angle;

if abs(epsi)>(pi/2)
    if epsi<(pi/2)
        epsi = psi - (angle - 2*pi);    
    else
        epsi = (psi - 2*pi) - angle;            
    end
end


s_start = ind_start-1; #--- HERE assuming all point are 1m distant
