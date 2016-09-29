function [OutTime, OutTraj, OutInput, coeffCost, coeffConst, Coeff_Vy, Coeff_Psi_dot, Coeff_Vx] = coeffConstraintCost(InTime, OldTrajectory,OldInput, InTraj, InInput,State_In, Input_In, s_start, s, s_target, v_target, mpcParams,LapNumber,Q,R,OldTrajectory_ID,OldInput_ID, OldInput_ID_1, OldTrajectory_ID_1, OldInput_1, OldTrajectory_1)
% this computes the coeff for the cost and constraints

% Outputs: 
% coeffConst
% coeffCost
MultiLapSysID = 1;

if MultiLapSysID == 0
    ActivateSysIdMultiLap = 10000;
else
    ActivateSysIdMultiLap = 2;
end
% Parameters
N = mpcParams.Hp;
nz = mpcParams.nz;
Order = mpcParams.OrderCostCons;
coeffCost     = zeros(2*(Order+1),1);
coeffConst    = zeros(2*(Order+1)*(nz-2),1); % nz-2 beacuse no coeff for s and filter
%coeffCost     = zeros((Order+1),1);
%coeffConst    = zeros((Order+1)*(nz-1),1);
Coeff_Psi_dot = zeros(3,1);
Coeff_Vy      = zeros(4,1);
Coeff_Vx      = zeros(3,1);
RMatrix = diag(R);

% Select the old data
OldxDot   = OldTrajectory(1,:);
OldyDot   = OldTrajectory(2,:);
OldpsiDot = OldTrajectory(3,:);
Oldpsi    = OldTrajectory(4,:);
OldeY     = OldTrajectory(5,:);
OldS      = OldTrajectory(6,:);

OldxDot_1   = OldTrajectory_1(1,:);
OldyDot_1   = OldTrajectory_1(2,:);
OldpsiDot_1 = OldTrajectory_1(3,:);
Oldpsi_1    = OldTrajectory_1(4,:);
OldeY_1     = OldTrajectory_1(5,:);
OldS_1      = OldTrajectory_1(6,:);
% 
N_points_1 = length(OldS_1);

N_points = length(OldS);
DelayBeforeStaringID = 4;
N_CurrentPoints = 50 + DelayBeforeStaringID;
PointsBefore = 80;
PointsAfter  = 160;
    
NumberPoints_Curr = 0;


OutTraj  = [InTraj(1:6, 2:(N_CurrentPoints) ), State_In(1:6,1)];
OutInput = [InInput(1:2, 2:(N_CurrentPoints) ), [Input_In(1,1); State_In(7,1)]];
OutTime  = InTime + 1;

OldxDot_ID   = OldTrajectory_ID(1,:);
OldyDot_ID   = OldTrajectory_ID(2,:);
OldpsiDot_ID = OldTrajectory_ID(3,:);
OldeY_ID     = OldTrajectory_ID(5,:);

N_points_ID = length(OldxDot_ID);

OldxDot_ID_1   = OldTrajectory_ID_1(1,:);
OldyDot_ID_1   = OldTrajectory_ID_1(2,:);
OldpsiDot_ID_1 = OldTrajectory_ID_1(3,:);
OldeY_ID_1     = OldTrajectory_ID_1(5,:);

N_points_ID_1 = length(OldxDot_ID_1);

Vx_INT_Curr      = OutTraj(1,:);
Vy_INT_Curr      = OutTraj(2,:);
Psi_dot_INT_Curr = OutTraj(3,:);
% eY_INT_Curr      = OutTraj(5,:);
delta_IND_Curr   = OutInput(1,:);
a_IND_Curr       = OutInput(2,:);
    
Qfunction = zeros(N_points,1);
Qfunction_1 = zeros(N_points_1,1);

if LapNumber > 1
    deltaOld = OldInput(1,:);
    aOld     = OldInput(2,:);

	deltaOld_1 = OldInput_1(1,:);
    aOld_1     = OldInput_1(2,:);
    
    deltaOld_ID = OldInput_ID(1,:);
    deltaOld_ID_1 = OldInput_ID_1(1,:);

    aOld_ID = OldInput_ID(2,:);
    aOld_ID_1 = OldInput_ID_1(2,:);
    % Compute the total s
    s_total = s_start + s;

    % Compute the index
    DistS = ( s_total - OldS ).^2;
    [~, idx_s] = min(DistS);
    IndexBezierS=idx_s(1);
    
    DistS_1 = ( s_total - OldS_1 ).^2;
    [~, idx_s_1] = min(DistS_1);
    IndexBezierS_1=idx_s_1(1);
    
    % Create the vector used for the interpolation
    bxDot_Vector     = OldxDot(IndexBezierS:IndexBezierS+floor(2*N))';
    byDot_Vector     = OldyDot(IndexBezierS:IndexBezierS+floor(2*N))';
	bpsiDot_Vector   = OldpsiDot(IndexBezierS:IndexBezierS+floor(2*N))';
    bpsi_Vector    = Oldpsi(IndexBezierS:IndexBezierS+floor(2*N))';
    beY_Vector     = OldeY(IndexBezierS:IndexBezierS+floor(2*N))';

    
	bxDot_Vector_1     = OldxDot_1(IndexBezierS_1:IndexBezierS_1+floor(2*N))';
    byDot_Vector_1     = OldyDot_1(IndexBezierS_1:IndexBezierS_1+floor(2*N))';
	bpsiDot_Vector_1   = OldpsiDot_1(IndexBezierS_1:IndexBezierS_1+floor(2*N))';
    bpsi_Vector_1    = Oldpsi_1(IndexBezierS_1:IndexBezierS_1+floor(2*N))';
    beY_Vector_1     = OldeY_1(IndexBezierS_1:IndexBezierS_1+floor(2*N))';
    % The states are parametrized with resprect to the curvilinear abscissa,
    % so we select the point used for the interpolation. Need to subtract an
    % offset to be coherent with the MPC formulation
    s_forinterpy = OldS(IndexBezierS:IndexBezierS+floor(2*N)) - s_start*ones(1,floor(2*N)+1);
    
    s_forinterpy_1 = OldS_1(IndexBezierS_1:IndexBezierS_1+floor(2*N)) - s_start*ones(1,floor(2*N)+1);
    
    % Create the Matrix for the interpolation
    MatrixInterp = zeros(floor(2*N)+1,Order+1);
    for i = 1:(floor(2*N)+1)
        for k = 0:Order
            MatrixInterp(i,Order+1-k) = [s_forinterpy(i)^k];
        end
    end

	MatrixInterp_1 = zeros(floor(2*N)+1,Order+1);
    for i = 1:(floor(2*N)+1)
        for k = 0:Order
            MatrixInterp_1(i,Order+1-k) = [s_forinterpy_1(i)^k];
        end
    end
    
    % Compute the coefficients
    CoefficientsFor_xDot   = MatrixInterp\bxDot_Vector;
    CoefficientsFor_yDot   = MatrixInterp\byDot_Vector;
    CoefficientsFor_psiDot   = MatrixInterp\bpsiDot_Vector;
    CoefficientsFor_psi    = MatrixInterp\bpsi_Vector;
    CoefficientsFor_ey     = MatrixInterp\beY_Vector;

	CoefficientsFor_xDot_1   = MatrixInterp_1\bxDot_Vector_1;
    CoefficientsFor_yDot_1   = MatrixInterp_1\byDot_Vector_1;
    CoefficientsFor_psiDot_1   = MatrixInterp_1\bpsiDot_Vector_1;
    CoefficientsFor_psi_1    = MatrixInterp_1\bpsi_Vector_1;
    CoefficientsFor_ey_1     = MatrixInterp_1\beY_Vector_1;
    
    % Stack the coefficients
    coeffConst(1:(Order+1),1)                 = CoefficientsFor_xDot;
    coeffConst(Order+1+(1:(Order+1)),1)       = CoefficientsFor_yDot;
    coeffConst(2*(Order+1)+(1:(Order+1)),1)   = CoefficientsFor_psiDot;    
    coeffConst(3*(Order+1)+(1:(Order+1)),1)   = CoefficientsFor_psi;    
    coeffConst(4*(Order+1)+(1:(Order+1)),1)   = CoefficientsFor_ey; 
    
    coeffConst(5*(Order+1)+(1:(Order+1)),1)   = CoefficientsFor_xDot_1;    
	coeffConst(6*(Order+1)+(1:(Order+1)),1)   = CoefficientsFor_yDot_1;    
	coeffConst(7*(Order+1)+(1:(Order+1)),1)   = CoefficientsFor_psiDot_1;    
	coeffConst(8*(Order+1)+(1:(Order+1)),1)   = CoefficientsFor_psi_1;        
	coeffConst(9*(Order+1)+(1:(Order+1)),1)   = CoefficientsFor_ey_1; 
    
    % Now compute the final cost coefficients

    PowerIndex = 0;
    for k=1:N_points
        % Start from the end ---> reverse the index
        indx = N_points-k+1;

        % Here the if for the cost, minimize the distance to a straight line
        if OldS(indx) <= s_target
            x = [1;  0; 0;  0; 0; 0];
            QMatrix = diag([100, 0, 0, 0, 0, 0]);
        else
            x = [1;  0; 0;  0; 0; 0];
            QMatrix = diag([0, 0, 0, 0, 0, 0]);
        end

        % Here actually computing the cost
        if (indx >=  IndexBezierS)
            if k ==1
                % If last point --> No Input
                Qfunction(indx) = (x')*QMatrix*(x);
            else
                % Stack input
                u = [deltaOld(indx); aOld(indx)];
                u = [0; 0]; % This is just for numerical issues
                % Compute the cost Q(N-1) = Q(N) + x'Qx + u'Ru
                Qfunction(indx) = Qfunction(indx + 1) + x'*QMatrix*x + u'*RMatrix*u;
                if Qfunction(indx) > 10^PowerIndex
                    PowerIndex = PowerIndex + 1;
                end
            end
            ender
    end

    % Select the part needed for the interpolation
    bQfunction_Vector    = Qfunction(IndexBezierS:IndexBezierS+floor(2*N));
    % Just scaling
    bQfunction_Vector_Scaled    = bQfunction_Vector;
    % Compute coefficient for the cost
    coeffCost(1:Order+1,1)  = MatrixInterp\bQfunction_Vector_Scaled;

    % Now compute the final cost coefficients
    PowerIndex = 0;
    for k=1:N_points_1
        % Start from the end ---> reverse the index
        indx_1 = N_points_1-k+1;

        % Here the if for the cost, minimize the distance to a straight line
        if OldS_1(indx_1) <= s_target
            x = [1;  0; 0;  0; 0; 0];
            QMatrix = diag([100, 0, 0, 0, 0, 0]);
        else
            x = [1;  0; 0;  0; 0; 0];
            QMatrix = diag([0, 0, 0, 0, 0, 0]);
        end

        % Here actually computing the cost
        if (indx_1 >=  IndexBezierS_1)
            if k ==1
                % If last point --> No Input
                Qfunction_1(indx_1) = (x')*QMatrix*(x);
            else
                % Stak input
                u = [deltaOld_1(indx_1); aOld_1(indx_1)];
                u = [0; 0]; % This is just for numerical issues
                % Compute the cost Q(N-1) = Q(N) + x'Qx + u'Ru
                Qfunction_1(indx_1) = Qfunction_1(indx_1 + 1) + x'*QMatrix*x + u'*RMatrix*u;
                if Qfunction_1(indx_1) > 10^PowerIndex
                    PowerIndex = PowerIndex + 1;
                end
            end
        end
    end

    % Select the part needed for the interpolation
    bQfunction_Vector_1    = Qfunction_1(IndexBezierS_1:IndexBezierS_1+floor(2*N));
    % Just scaling
    bQfunction_Vector_Scaled_1    = bQfunction_Vector_1;
    % Compute coefficient for the cost
    coeffCost(Order+1+(1:(Order+1)),1)  = MatrixInterp_1\bQfunction_Vector_Scaled_1;
    
    % =================================================================== %
    % ========================== System ID ============================== %
    % =================================================================== %
    
    %Decide how many points are needed

    Before = IndexBezierS-PointsBefore;
    After  = IndexBezierS+PointsAfter;
    
    if LapNumber > ActivateSysIdMultiLap
        Before_1 = IndexBezierS_1 - PointsBefore;
        After_1  = IndexBezierS_1 + PointsAfter;
        
        if Before_1<=0
            After_1=After_1-Before_1;
            Before_1=1;
        end
        if After_1 > N_points_ID_1
            Before_1 = Before_1 - ( After_1 -N_points_ID_1 );
            After_1  = N_points_ID_1;
        end
    else
        Before_1 = Before;
        After_1  = After;
    end
    
    if Before<=0
        After=After-Before;
        Before=1;
    end
    if After > N_points_ID
        Before = Before - ( After -N_points_ID );
        After= N_points_ID;
    end

    % Select Points
    Psi_dot_INT = OldpsiDot_ID(Before:After);
    Vy_INT      = OldyDot_ID(Before:After);
    Vx_INT      = OldxDot_ID(Before:After);
    eY_INT      = OldeY_ID(Before:After);
    delta_IND   = deltaOld_ID(Before:After);
    a_IND       = aOld_ID(Before:After);

    if LapNumber > ActivateSysIdMultiLap
        Psi_dot_INT_1 = OldpsiDot_ID_1(Before_1:After_1);
        Vy_INT_1      = OldyDot_ID_1(Before_1:After_1);
        Vx_INT_1      = OldxDot_ID_1(Before_1:After_1);
        eY_INT_1      = OldeY_ID_1(Before_1:After_1);
        delta_IND_1   = deltaOld_ID_1(Before_1:After_1);
        a_IND_1       = aOld_ID_1(Before_1:After_1);
    else
        Psi_dot_INT_1 = Psi_dot_INT;
        Vy_INT_1      = Vy_INT;
        Vx_INT_1      = Vx_INT;
        eY_INT_1      = eY_INT;
        delta_IND_1   = delta_IND;
        a_IND_1       = a_IND;
    end
    
    % Interpolate the Yaw Rate
    NumberPoints = length(Vx_INT);
    NumberPoints_1 = length(Vx_INT_1);
    M_Psi_dot_INT = zeros(NumberPoints-2+1 + NumberPoints_1-2+1 + N_CurrentPoints-DelayBeforeStaringID+1, 3);
    y_Psi_dot_INT = zeros(NumberPoints-2+1 + NumberPoints_1-2+1 + N_CurrentPoints-DelayBeforeStaringID+1, 1);

    M_Vy_INT = zeros(NumberPoints-2+1 + NumberPoints_1-2+1 + N_CurrentPoints-DelayBeforeStaringID+1, 4);
    y_Vy_INT = zeros(NumberPoints-2+1 + NumberPoints_1-2+1 + N_CurrentPoints-DelayBeforeStaringID+1, 1);

    M_Vx_INT = zeros(NumberPoints-2+1 + NumberPoints_1-2+1 + N_CurrentPoints-DelayBeforeStaringID+1, 3);
    y_Vx_INT = zeros(NumberPoints-2+1 + NumberPoints_1-2+1 + N_CurrentPoints-DelayBeforeStaringID+1, 1);
    
    % Identify Psidot
    ii = 1;
    for j = 0:(NumberPoints-2)
        ind = NumberPoints - j;
        M_Psi_dot_INT(ii,1)   = [Psi_dot_INT(ind-1) / Vx_INT(ind-1)];
        M_Psi_dot_INT(ii,2)   = [Vy_INT(ind-1) / Vx_INT(ind-1)];
        M_Psi_dot_INT(ii,3)   = [delta_IND(ind-1)];
        
        y_Psi_dot_INT(ii,1) = Psi_dot_INT(ind) - Psi_dot_INT(ind-1);
        ii = ii + 1;
    end
    for j = 0:(NumberPoints_1-2)
        ind = NumberPoints_1 - j;
        M_Psi_dot_INT(ii,1)   = [Psi_dot_INT_1(ind-1) / Vx_INT_1(ind-1)];
        M_Psi_dot_INT(ii,2)   = [Vy_INT_1(ind-1) / Vx_INT_1(ind-1)];
        M_Psi_dot_INT(ii,3)   = [delta_IND_1(ind-1)];

        y_Psi_dot_INT(ii,1) = Psi_dot_INT_1(ind) - Psi_dot_INT_1(ind-1);
        ii = ii + 1;
    end
    
    if OutTime > DelayBeforeStaringID
        NumberPoints_Curr = (min(OutTime, N_CurrentPoints)-DelayBeforeStaringID);
        for j = 0:NumberPoints_Curr
            ind = N_CurrentPoints - j;
            M_Psi_dot_INT(ii,1)   = [Psi_dot_INT_Curr(ind-1) / Vx_INT_Curr(ind-1)];
            M_Psi_dot_INT(ii,2)   = [Vy_INT_Curr(ind-1) / Vx_INT_Curr(ind-1)];
            M_Psi_dot_INT(ii,3)   = [delta_IND_Curr(ind)];

            y_Psi_dot_INT(ii,1) = Psi_dot_INT_Curr(ind) - Psi_dot_INT_Curr(ind-1);
            ii = ii + 1;
        end
    end
    
    % Identify Vy
    ii = 1;
    for j = 0:(NumberPoints-2)
        ind = NumberPoints - j;
        M_Vy_INT(ii,1)   = [Vy_INT(ind-1) / Vx_INT(ind-1)];
        M_Vy_INT(ii,2:3) = [Psi_dot_INT(ind-1) * Vx_INT(ind-1)  Psi_dot_INT(ind-1)/Vx_INT(ind-1)];
        M_Vy_INT(ii,4)   = [delta_IND(ind-1)];

        y_Vy_INT(ii,1) = Vy_INT(ind) - Vy_INT(ind-1);
        ii = ii + 1;
    end
    for j = 0:(NumberPoints_1-2)
        ind = NumberPoints_1 - j;
        M_Vy_INT(ii,1)   = [Vy_INT_1(ind-1) / Vx_INT_1(ind-1)];
        M_Vy_INT(ii,2:3) = [Psi_dot_INT_1(ind-1) * Vx_INT_1(ind-1)  Psi_dot_INT_1(ind-1)/Vx_INT_1(ind-1)];
        M_Vy_INT(ii,4)   = [delta_IND_1(ind-1)];

        y_Vy_INT(ii,1) = Vy_INT_1(ind) - Vy_INT_1(ind-1);
        ii = ii + 1;
    end
    
    if OutTime > DelayBeforeStaringID
        for j = 0:NumberPoints_Curr
            ind = N_CurrentPoints - j;
            M_Vy_INT(ii,1)   = [Vy_INT_Curr(ind-1) / Vx_INT_Curr(ind-1)];
            M_Vy_INT(ii,2:3) = [Psi_dot_INT_Curr(ind-1) * Vx_INT_Curr(ind-1)  Psi_dot_INT_Curr(ind-1)/Vx_INT_Curr(ind-1)];
            M_Vy_INT(ii,4)   = [delta_IND_Curr(ind)];

            y_Vy_INT(ii,1) = Vy_INT_Curr(ind) - Vy_INT_Curr(ind-1);
            ii = ii + 1;
        end
    end

    % Identify Vx
    ii = 1;
    for j = 0:(NumberPoints-2)
        ind = NumberPoints - j;
        M_Vx_INT(ii,1)   = [Vy_INT(ind-1)];
        M_Vx_INT(ii,2)   = [Psi_dot_INT(ind-1)];
        M_Vx_INT(ii,3)   = [Vx_INT(ind-1)];

        y_Vx_INT(ii,1) = Vx_INT(ind) - a_IND(ind-1);
        ii = ii + 1;
    end
    for j = 0:(NumberPoints_1-2)
        ind = NumberPoints_1 - j;
        M_Vx_INT(ii,1)   = [Vy_INT_1(ind-1)];
        M_Vx_INT(ii,2)   = [Psi_dot_INT_1(ind-1)];
        M_Vx_INT(ii,3)   = [Vx_INT_1(ind-1)];

        y_Vx_INT(ii,1) = Vx_INT_1(ind) - a_IND_1(ind-1);
        ii = ii + 1;
    end
    
    if OutTime > DelayBeforeStaringID
        for j = 0:NumberPoints_Curr
            ind = N_CurrentPoints - j;
            M_Vx_INT(ii,1)   = [Vy_INT_Curr(ind-1)];
            M_Vx_INT(ii,2)   = [Psi_dot_INT_Curr(ind-1)];
            M_Vx_INT(ii,3)   = [Vx_INT_Curr(ind-1)];

            y_Vx_INT(ii,1) = Vx_INT_Curr(ind) - a_IND_Curr(ind);
            ii = ii + 1;
        end
    end
    
    Coeff_Psi_dot(:,1) = (transpose(M_Psi_dot_INT)*M_Psi_dot_INT) \ transpose(M_Psi_dot_INT)*y_Psi_dot_INT;
    Coeff_Vy(:,1)      = (transpose(   M_Vy_INT  )*M_Vy_INT     ) \ transpose(M_Vy_INT     )*y_Vy_INT;
    Coeff_Vx(:,1)      = (transpose(   M_Vx_INT  )*M_Vx_INT     ) \ transpose(M_Vx_INT     )*y_Vx_INT;

else        % if lap number == 1 
    coeffCost(:,1)     = zeros(2*( Order+1 ),1);
    coeffConst         = zeros(2*( (Order+1)*(nz-2) ),1); % nz-2 beacuse no coeff for s and filter
%     coeffCost(:,1)     = zeros(( Order+1 ),1);
%     coeffConst         = zeros(( (Order+1)*(nz-1) ),1);
    Coeff_Psi_dot(:,1) = zeros(3,1);
    Coeff_Vy(:,1)      = zeros(4,1);
    Coeff_Vx(:,1)      = zeros(3,1);
end

