# This function evaluates and returns the coefficients for constraints and cost which are used later in the MPC formulation
# Inputs are:
# oldTraj   -> contains information about previous trajectories and Inputs
# lapStatus -> contains information about the current lap
# mpcCoeff  -> contains information about previous coefficient results
# posInfo   -> contains information about the car's current position along the track
# mpcParams -> contains information about the MPC formulation (e.g. Q, R)
# stateIn   -> current state of the car
# inputIn   -> last input to the system

function coeffConstraintCost(oldTraj, lapStatus, mpcCoeff, posInfo, mpcParams, stateIn, inputIn)
    # this computes the coefficients for the cost and constraints

    # Outputs: 
    # coeffConst
    # coeffCost
    MultiLapSysID = 1

    if MultiLapSysID == 0
        ActivateSysIdMultiLap = 10000
    else
        ActivateSysIdMultiLap = 2
    end

    # Read Inputs
    oldTrajectory       = oldTraj.oldTraj
    oldTrajectory_1     = oldTraj.oldTraj_1
    oldInput            = oldTraj.oldInput
    oldInput_1          = oldTraj.oldInput_1
    oldTrajectory_ID    = oldTraj.oldTrajectory_ID
    oldInput_ID         = oldTraj.oldInput_ID
    oldInput_ID_1       = oldTraj.oldInput_ID
    oldTrajectory_ID_1  = oldTraj.oldTrajectory_ID_1

    inTime              = mpcCoeff.currentTime
    inTraj              = mpcCoeff.currentTraj
    inInput             = mpcCoeff.currentInput

    lapNumber           = lapStatus.currentLap

    s_start             = posInfo.s_start
    s                   = posInfo.s
    s_target            = posInfo.s_target


    # Parameters
    N               = mpcParams.N
    nz              = mpcParams.nz
    Order           = mpcParams.OrderCostCons         # interpolation order
    R               = mpcParams.R

    coeffCost       = zeros(2*(Order+1),1)            # polynomial coefficients for cost
    coeffConst      = zeros(2*(Order+1)*(nz-2),1)     # nz-2 beacuse no coeff for s and filter

    Coeff_Psi_dot   = zeros(3,1)
    Coeff_Vy        = zeros(4,1)
    Coeff_Vx        = zeros(3,1)
    RMatrix         = diagm(R)

    # Select the old data
    oldxDot         = oldTrajectory[1,:]
    oldyDot         = oldTrajectory[2,:]
    oldpsiDot       = oldTrajectory[3,:]
    oldPsi          = oldTrajectory[4,:]
    oldeY           = oldTrajectory[5,:]
    oldS            = oldTrajectory[6,:]

    oldxDot_1       = oldTrajectory_1[1,:]
    oldyDot_1       = oldTrajectory_1[2,:]
    oldpsiDot_1     = oldTrajectory_1[3,:]
    oldpsi_1        = oldTrajectory_1[4,:]
    oldeY_1         = oldTrajectory_1[5,:]
    oldS_1          = oldTrajectory_1[6,:]

    N_points_1      = length(oldS_1)
    N_points        = length(oldS)

    DelayBeforeStaringID    = 4
    N_CurrentPoints         = 50 + DelayBeforeStaringID
    PointsBefore            = 80
    PointsAfter             = 160
        
    NumberPoints_Curr       = 0


    OutTraj  = [inTraj[1:6, 2:(N_CurrentPoints) ] stateIn[1:6,1]]
    OutInput = [inInput[1:2, 2:(N_CurrentPoints) ] [inputIn[1,1]; stateIn[7,1]]]
    OutTime  = inTime + 1


    # oldxDot_ID      = oldTrajectory_ID[1,:]
    # oldyDot_ID      = oldTrajectory_ID[2,:]
    # oldpsiDot_ID    = oldTrajectory_ID[3,:]
    # oldeY_ID        = oldTrajectory_ID[5,:]

    # N_points_ID     = length(oldxDot_ID)

    # oldxDot_ID_1    = oldTrajectory_ID_1[1,:]
    # oldyDot_ID_1    = oldTrajectory_ID_1[2,:]
    # oldpsiDot_ID_1  = oldTrajectory_ID_1[3,:]
    # oldeY_ID_1      = oldTrajectory_ID_1[5,:]

    # N_points_ID_1   = length(oldxDot_ID_1)

    Vx_INT_Curr      = OutTraj[1,:]
    Vy_INT_Curr      = OutTraj[2,:]
    Psi_dot_INT_Curr = OutTraj[3,:]
    # eY_INT_Curr      = OutTraj[5,:]
    delta_IND_Curr   = OutInput[1,:]
    a_IND_Curr       = OutInput[2,:]
        
    Qfunction   = zeros(N_points,1)
    Qfunction_1 = zeros(N_points_1,1)


    if lapNumber > 1
        deltaOld = oldInput[1,:]
        aOld     = oldInput[2,:]

        deltaOld_1 = oldInput_1[1,:]
        aOld_1     = oldInput_1[2,:]
        
        # deltaOld_ID = oldInput_ID[1,:]
        # deltaOld_ID_1 = oldInput_ID_1[1,:]

        # aOld_ID = oldInput_ID[2,:]
        # aOld_ID_1 = oldInput_ID_1[2,:]
        
        # Compute the total s
        s_total = s_start + s

        # Compute the index
        DistS = ( s_total - oldS ).^2
        idx_s = indmin(DistS)
        IndexBezierS=idx_s[1]
        
        DistS_1 = ( s_total - oldS_1 ).^2
        idx_s_1 = indmin(DistS_1)
        IndexBezierS_1=idx_s_1[1]
        
        # Create the vector used for the interpolation
        bxDot_Vector    = oldxDot[IndexBezierS:IndexBezierS+floor(2*N)]'
        byDot_Vector    = oldyDot[IndexBezierS:IndexBezierS+floor(2*N)]'
        bpsiDot_Vector  = oldpsiDot[IndexBezierS:IndexBezierS+floor(2*N)]'
        bpsi_Vector     = oldPsi[IndexBezierS:IndexBezierS+floor(2*N)]'
        beY_Vector      = oldeY[IndexBezierS:IndexBezierS+floor(2*N)]'

        
        bxDot_Vector_1   = oldxDot_1[IndexBezierS_1:IndexBezierS_1+floor(2*N)]'
        byDot_Vector_1   = oldyDot_1[IndexBezierS_1:IndexBezierS_1+floor(2*N)]'
        bpsiDot_Vector_1 = oldpsiDot_1[IndexBezierS_1:IndexBezierS_1+floor(2*N)]'
        bpsi_Vector_1    = oldpsi_1[IndexBezierS_1:IndexBezierS_1+floor(2*N)]'
        beY_Vector_1     = oldeY_1[IndexBezierS_1:IndexBezierS_1+floor(2*N)]'
        
        # The states are parametrized with resprect to the curvilinear abscissa,
        # so we select the point used for the interpolation. Need to subtract an
        # offset to be coherent with the MPC formulation
        s_forinterpy   = oldS[IndexBezierS:IndexBezierS+floor(2*N)] - s_start*ones(1,floor(2*N)+1)
        s_forinterpy_1 = oldS_1[IndexBezierS_1:IndexBezierS_1+floor(2*N)] - s_start*ones(1,floor(2*N)+1)
        
        # Create the Matrix for the interpolation
        MatrixInterp = zeros(floor(2*N)+1,Order+1)
        for i = 1:(floor(2*N)+1)
            for k = 0:Order
                MatrixInterp[i,Order+1-k] = [s_forinterpy[i]^k]
            end
        end

        MatrixInterp_1 = zeros(floor(2*N)+1,Order+1)
        for i = 1:(floor(2*N)+1)
            for k = 0:Order
                MatrixInterp_1[i,Order+1-k] = [s_forinterpy_1[i]^k]
            end
        end
        
        # Compute the coefficients
        CoefficientsFor_xDot   = MatrixInterp\bxDot_Vector
        CoefficientsFor_yDot   = MatrixInterp\byDot_Vector
        CoefficientsFor_psiDot = MatrixInterp\bpsiDot_Vector
        CoefficientsFor_psi    = MatrixInterp\bpsi_Vector
        CoefficientsFor_ey     = MatrixInterp\beY_Vector

        CoefficientsFor_xDot_1   = MatrixInterp_1\bxDot_Vector_1
        CoefficientsFor_yDot_1   = MatrixInterp_1\byDot_Vector_1
        CoefficientsFor_psiDot_1 = MatrixInterp_1\bpsiDot_Vector_1
        CoefficientsFor_psi_1    = MatrixInterp_1\bpsi_Vector_1
        CoefficientsFor_ey_1     = MatrixInterp_1\beY_Vector_1
        
        # Stack the coefficients
        coeffConst[1:(Order+1),1]                 = CoefficientsFor_xDot
        coeffConst[Order+1+(1:(Order+1)),1]       = CoefficientsFor_yDot
        coeffConst[2*(Order+1)+(1:(Order+1)),1]   = CoefficientsFor_psiDot    
        coeffConst[3*(Order+1)+(1:(Order+1)),1]   = CoefficientsFor_psi    
        coeffConst[4*(Order+1)+(1:(Order+1)),1]   = CoefficientsFor_ey 
        
        coeffConst[5*(Order+1)+(1:(Order+1)),1]   = CoefficientsFor_xDot_1    
        coeffConst[6*(Order+1)+(1:(Order+1)),1]   = CoefficientsFor_yDot_1    
        coeffConst[7*(Order+1)+(1:(Order+1)),1]   = CoefficientsFor_psiDot_1    
        coeffConst[8*(Order+1)+(1:(Order+1)),1]   = CoefficientsFor_psi_1        
        coeffConst[9*(Order+1)+(1:(Order+1)),1]   = CoefficientsFor_ey_1 
        
        # Now compute the final cost coefficients

        PowerIndex = 0
        for k=1:N_points
            # Start from the end ---> reverse the index
            indx = N_points-k+1;

            # Here the if for the cost, minimize the distance to a straight line
            if oldS(indx) <= s_target
                x = [1; 0; 0; 0; 0; 0]
                QMatrix = diagm([100, 0, 0, 0, 0, 0])
            else
                x = [1; 0; 0; 0; 0; 0]
                QMatrix = diagm([0, 0, 0, 0, 0, 0])
            end

            # Here actually computing the cost
            if (indx >=  IndexBezierS)
                if k ==1
                    # If last point --> No Input
                    Qfunction[indx] = (x')*QMatrix*(x)
                else
                    # Stack input
                    u = [deltaOld[indx]; aOld[indx]]
                    u = [0; 0]; # % This is just for numerical issues
                    # Compute the cost Q(N-1) = Q(N) + x'Qx + u'Ru
                    Qfunction[indx] = Qfunction[indx + 1] + x'*QMatrix*x + u'*RMatrix*u
                    if Qfunction[indx] > 10^PowerIndex
                        PowerIndex = PowerIndex + 1
                    end
                end
            end
        end

        # Select the part needed for the interpolation
        bQfunction_Vector    = Qfunction[IndexBezierS:IndexBezierS+floor(2*N)]
        # Just scaling
        bQfunction_Vector_Scaled    = bQfunction_Vector
        # Compute coefficient for the cost
        coeffCost[1:Order+1,1]  = MatrixInterp\bQfunction_Vector_Scaled

        # Now compute the final cost coefficients

        PowerIndex = 0;
        for k=1:N_points_1
            # Start from the end ---> reverse the index
            indx_1 = N_points_1-k+1

            # Here the if for the cost, minimize the distance to a straight line
            if oldS_1[indx_1] <= s_target
                x = [1; 0; 0; 0; 0; 0]
                QMatrix = diagm([100, 0, 0, 0, 0, 0])
            else
                x = [1; 0; 0; 0; 0; 0]
                QMatrix = diagm([0, 0, 0, 0, 0, 0])
            end

            # Here actually computing the cost
            if (indx_1 >=  IndexBezierS_1)
                if k ==1
                    # If last point --> No Input
                    Qfunction_1[indx_1] = (x')*QMatrix*(x)
                else
                    # Stack input
                    u = [deltaOld_1[indx_1]; aOld_1[indx_1]]
                    u = [0; 0] # This is just for numerical issues
                    # Compute the cost Q(N-1) = Q(N) + x'Qx + u'Ru
                    Qfunction_1[indx_1] = Qfunction_1[indx_1 + 1] + x'*QMatrix*x + u'*RMatrix*u
                    if Qfunction_1[indx_1] > 10^PowerIndex
                        PowerIndex = PowerIndex + 1
                    end
                end
            end
        end

        # Select the part needed for the interpolation
        bQfunction_Vector_1    = Qfunction_1[IndexBezierS_1:IndexBezierS_1+floor(2*N)]
        # Just scaling
        bQfunction_Vector_Scaled_1    = bQfunction_Vector_1
        # Compute coefficient for the cost
        coeffCost[Order+1+(1:(Order+1)),1] = MatrixInterp_1\bQfunction_Vector_Scaled_1

        # ---------------
        # HERE STARTS THE SYSTEM ID part (skipped for now)
        # ---------------

    else        # if it is the first lap
        coeffCost[:,1]     = zeros(2*( Order+1 ),1)
        coeffConst         = zeros(2*( (Order+1)*(nz-2) ),1) # nz-2 because no coeff for s and filter
    #     coeffCost(:,1)     = zeros(( Order+1 ),1)
    #     coeffConst         = zeros(( (Order+1)*(nz-1) ),1)
        Coeff_Psi_dot[:,1] = zeros(3,1)
        Coeff_Vy[:,1]      = zeros(4,1)
        Coeff_Vx[:,1]      = zeros(3,1)
    end

    return MpcCoeff(OutTime, OutTraj, OutInput, coeffCost, coeffConst, Coeff_Vx,  Coeff_Vy, Coeff_Psi_dot)
    # Todo: Have to find expressions for Coeff_Vy, Coeff_Vx, Coeff_Psi_dot without System ID
end
