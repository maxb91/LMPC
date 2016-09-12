type mpcParams
    Hp
    nz
    OrderCostCons
end

function coeffConstraintCost(InTime, OldTrajectory,OldInput, InTraj, InInput,State_In, Input_In, s_start, s, s_target, v_target, mpcParams,LapNumber,Q,R,OldTrajectory_ID,OldInput_ID, OldInput_ID_1, OldTrajectory_ID_1, OldInput_1, OldTrajectory_1)
    # this computes the coeff for the cost and constraints

    # Outputs: 
    # coeffConst
    # coeffCost
    MultiLapSysID = 1;

    if MultiLapSysID == 0
        ActivateSysIdMultiLap = 10000;
    else
        ActivateSysIdMultiLap = 2;
    end
    # Parameters
    N = mpcParams.Hp;
    nz = mpcParams.nz;
    Order = mpcParams.OrderCostCons;
    coeffCost     = zeros(2*(Order+1),1);
    coeffConst    = zeros(2*(Order+1)*(nz-2),1); #  nz-2 beacuse no coeff for s and filter
    # %coeffCost     = zeros((Order+1),1);
    # %coeffConst    = zeros((Order+1)*(nz-1),1);
    Coeff_Psi_dot = zeros(3,1);
    Coeff_Vy      = zeros(4,1);
    Coeff_Vx      = zeros(3,1);
    RMatrix = diag(R);

    # % Select the old data
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
    # % eY_INT_Curr      = OutTraj(5,:);
    delta_IND_Curr   = OutInput(1,:);
    a_IND_Curr       = OutInput(2,:);
        
    Qfunction = zeros(N_points,1);
    Qfunction_1 = zeros(N_points_1,1);